from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode
import networkx as nx
import warnings
from copy import deepcopy

class RiboGraphFlux(RiboGraph):
    """
    A directed graph representing ribosomal flux through phase space, built from a TransitionMap.

    RiboGraphFlux models the movement of ribosomes along an mRNA transcript by propagating
    flux from a bulk cytoplasmic pool through a network of ribosomal phase-space nodes.
    Each edge carries flux values that account for initiation, elongation, termination,
    reinitiation, and decay during scanning and translation.

    The graph is constructed automatically on instantiation via construct(), which:
      - Removes edges below the weight cutoff from the transition map
      - Normalises entry flux across all initiation sites
      - Propagates flux iteratively through the graph, applying decay and branching
      - Collapses redundant intermediate nodes

    Parameters
    ----------
    transition_map : TransitionMap
        The transition map from which the flux graph is built. Encodes the probability
        of moving between ribosomal phase-space nodes.
    half_life_scanning : float or None
        The half-life (in nucleotides) of a scanning 40S subunit. Controls the rate of
        decay applied to flux along scanning edges (phase == 0). If None, no scanning
        decay is applied.
    half_life_translation : float or None
        The half-life (in nucleotides) of an elongating ribosome. Controls
        decay along translation edges (phase > 0). If None, no translational decay is applied.
    weight_cutoff : float
        Edges in the transition map with weight below this value are removed before
        flux propagation. Default is 0.0.
    reinitiation_half_life : float or None
        Controls the decay of reinitiation potential as a function of distance from the
        initiation node. Larger values allow reinitiation over greater distances.
        If None, reinitiation potential does not decay.
    ternary_complex_half_life : float or None
        Half-life governing the replenishment of ternary complex during scanning after
        a termination event. Used to scale initiation probability following reinitiation.
        If None, ternary complex is assumed to be fully available at all times.
    flux_cutoff : float
        Flux values below this threshold are not propagated further, pruning negligible
        branches from the graph. Default is 0.001.
    retention_limit : int or None
        Maximum number of consecutive retention events (ribosome returning to a scanning state on the same
        mRNA after termination) allowed per path. If None, retention is unlimited.
        Default is 1.

    Attributes
    ----------
    transitions : TransitionMap
        The underlying transition map used to build the graph.
    ribopaths : list of list of tuple
        All paths representing continuous 40S association with the mRNA, from the
        first loading node to the bulk node. Each path is a list of edge tuples (u, v).
    translons : list of list of tuple
        All sub-paths within ribopaths representing continuous 60S association
        (i.e., active translation). Each translon is a list of edge tuples (u, v)
        where both nodes have phase > 0.

    Edge Data
    ---------
    Each edge (u, v) in the graph carries the following data:
      flux_start : float
          Flux entering the edge at node u.
      flux_end : float
          Flux remaining at node v after decay along the edge.
      decay : float
          Flux lost along the edge due to ribosome drop-off.


    Methods
    -------
    flux_proportion(u, v) -> float
        Returns the total proportion of bulk flux that travels between nodes u and v,
        summed across all simple paths connecting them.
    flux_proportion_path(path) -> float
        Returns the proportion of bulk flux that travels along a specific path,
        given as an ordered list of RiboNodes.
    edge_weight(u, v) -> float
        Returns the fraction of flux leaving node u that is carried by the edge to v.
    edge_decay(u, v) -> float
        Returns the fraction of flux lost to ribosome drop-off along the edge from u to v.
    rein_decay(u, v) -> float
        Returns the decay factor applied to reinitiation potential between nodes u and v.
    ternary_complex_proportion(u, v) -> float
        Returns the fraction of ternary complex available at v, given prior scanning from u.
    add_transition(source, target, probability)
        Adds a single transition to the underlying map and reconstructs the flux graph.
    add_transitions_from(tbunch)
        Adds multiple transitions from an iterable of (source, target, weight) tuples
        and reconstructs the flux graph.

    Notes
    -----
    - Node phase conventions: phase == -1 represents the bulk cytoplasmic pool (off-mRNA),
      phase == 0 represents a scanning 40S subunit, and phase > 0 represents an elongating
      80S ribosome.
    - Flux is normalised so that the maximum edge flux equals 1.0.
    - Floating-point accumulation errors in complex graphs may cause small discrepancies
      between total inbound and outbound flux; a warning is raised if this exceeds the
      internal tolerance (flux_error = 1e-15).
    """

    def __init__(self, transition_map: TransitionMap, 
                 incoming_graph_data=None, 
                 half_life_scanning: float|None = None, 
                 half_life_translation: float|None = None, 
                 weight_cutoff=0.0, 
                 reinitiation_half_life = None,
                 ternary_complex_half_life = None,
                 flux_cutoff = 0.001,
                 retention_limit: int|None = 1,
                 **attr):
        

        super().__init__(incoming_graph_data, **attr)

        self.transitions = transition_map
        self.begun = False
        if incoming_graph_data is not None:
            raise ValueError('Incoming graph data must be left empty, graph is calculated from transition_map')
        self.half_life_translation = half_life_translation
        self.half_life_scanning = half_life_scanning
        self.weight_cutoff = weight_cutoff
        self.reinitiation_potential=reinitiation_half_life
        self.tc_half_life = ternary_complex_half_life
        self.flux_cutoff = flux_cutoff
        self.flux_error = 0.000000000000001
        self.retention_limit = retention_limit
        if map:
            self.construct()   
    
    def construct(self):
        below_cutoff = []

        for u, v in self.transitions.edges:
            if self.transitions[u][v]['weight'] < self.weight_cutoff:
                below_cutoff.append((u,v))
        self.transitions.remove_edges_from(below_cutoff)

        nodes_to_remove = [node for node, degree in self.transitions.degree() if degree == 0]

        self.transitions.remove_nodes_from(nodes_to_remove)
        total_weight_in = 0
        weights = []
        for u,v in self.transitions.edges:
            if u.phase == -1:
                total_weight_in += self.transitions[u][v]['weight']
                weights.append(self.transitions[u][v]['weight'])
        norm_weight = {}

        for weight in weights:
            norm_weight[weight] = weight/total_weight_in

            
                
        
        for u, v in self.transitions.edges:

            if u.phase == -1:
                self.add_edge(self.bulk_node, u)
                flux = norm_weight[self.transitions[u][v]['weight']]
                self.add_edge(self.bulk_node, u, weight=flux, flux_start=flux, flux_end=flux)
                self.add_edge(u, v, weight=flux, flux_start=flux, flux_end=flux )

                self._iterate_graph(v, flux, initiation_node=v) 
                self._normalize_flux()

        self._collapse_unused_nodes()
        self._is_valid()


    def _downstream_node(self, node: RiboNode):
        if not any([p.position > node.position 
                            and p.phase == node.phase
                            for p in self.transitions.nodes ]):
            return None
            
        return RiboNode(( min([p.position
                            for p in self.transitions.nodes 
                            if p.position > node.position 
                            and p.phase == node.phase]),

                            node.phase))
    

    def _iterate_graph(self, node: RiboNode, flux, retained=0,
                   initiation_node: RiboNode | None = None,
                   retention_node: RiboNode | None = None):
        """
        Builds the graph using an explicit stack
        """

        stack = [(node, flux, retained, initiation_node, retention_node)]

        while stack:
            node, flux, retained, initiation_node, retention_node = stack.pop()

            if node == self.bulk_node:
                continue

            next_node = self._downstream_node(node)

            if next_node is None:
                drop_node = RiboNode(node.position, -1)
                self.add_edge(node, drop_node, flux_start=flux, flux_end=flux)
                self.add_edge(drop_node, self.bulk_node, flux_start=flux, flux_end=flux)
                continue

            # ── Decay ──────────────────────────────────────────────────────────
            
            drop_flux = flux * self.edge_decay(node, next_node)
            remaining_flux = flux - drop_flux

            if drop_flux != 0:
                drop_node = RiboNode(next_node.position, -1)
                self.add_edge(next_node, drop_node, flux_start=drop_flux, flux_end=drop_flux)
                self.add_edge(drop_node, self.bulk_node, flux_start=drop_flux, flux_end=drop_flux)

            self.add_edge(node, next_node,
                        flux_start=flux, flux_end=remaining_flux, decay=drop_flux)


            # ── Pass 1: retention edges ────────────────────────────────────────
            for u, v, w in self.transitions.out_edges(next_node, data='weight'):
                if not self.is_retention(u, v):
                    continue

                if initiation_node is None:
                    continue

                retention_flux = remaining_flux * self.rein_decay(initiation_node, next_node)
                if retention_flux < self.flux_cutoff:
                    continue
                if self.retention_limit is not None and retained >= self.retention_limit:
                    continue

                new_flux = retention_flux * w

                if new_flux < self.flux_cutoff:
                    continue
                remaining_flux -= new_flux


                retention_node    = v

                self.add_edge(u, v, flux_start=new_flux, flux_end=new_flux)
                stack.append((v, new_flux, retained + 1, initiation_node, retention_node))

            # ── Pass 2: initiation / termination / shift edges ─────────────────────────
            for u, v, w in self.transitions.out_edges(next_node, data='weight'):
                if self.is_retention(u, v):
                    continue
    
                if self.is_initiation(u, v):
                    if retained:
                        w *= self.ternary_complex_proportion(retention_node, next_node)
                    initiation_node = v

                new_flux = remaining_flux * w

                if new_flux < self.flux_cutoff and v.phase != -1:
                    continue

                remaining_flux -= new_flux


                self.add_edge(u, v, flux_start=new_flux, flux_end=new_flux)

                if v.phase == -1:
                    self.add_edge(v, self.bulk_node, flux_start=new_flux, flux_end=new_flux)
                    continue

                stack.append((v, new_flux, retained, initiation_node, retention_node))

                

            # ── Continuation on same phase ─────────────────────────────────────
            if remaining_flux < -self.flux_error:
                raise ValueError(
                    f'Weight from node: {next_node} exceeds 1: {abs(remaining_flux - 1)}')

            if remaining_flux >= self.flux_error:
                stack.append((next_node, remaining_flux,
                            retained, initiation_node, retention_node))
                

    def _collapse_unused_nodes(self):
        changed = True
        test_graph = deepcopy(self)
        test_graph.remove_node(self.bulk_node)
    

        while changed:
            changed = False
            for node in list(nx.topological_sort(test_graph)):
                if node.phase == -1:
                    continue

                in_edges = self.in_edges(node)
                #ignore nodes with more than 1 in edge
                if len(list(in_edges)) > 1:
                    continue

                in_u = False
                out_v = False
                for u,_ in in_edges:
                    if u.phase == node.phase:
                        in_u = u
                        break

                out_edges = self.out_edges(node)
                for _,v in out_edges:
                    if v.phase == node.phase:
                        out_v = v
                        break

                if not in_u or not out_v:
                    continue
                
                
                

                in_flux_end = self[in_u][node]['flux_end']
                
                out_flux = self[node][out_v]['flux_start']
                in_flux_start = self[in_u][node]['flux_start']
                out_decay = self[node][out_v]['flux_start'] - self[node][out_v]['flux_end']

                if abs(in_flux_end - out_flux) < self.flux_error:

                    self.remove_node(node)


                    drop_flux = in_flux_start * self.edge_decay(in_u, out_v)
                    endflux = in_flux_start - drop_flux

                    if drop_flux != 0:
                        drop_node = RiboNode(out_v.position, -1)
                        out_drop_flux = self[out_v][drop_node]['flux_start']
                        non_decay_drop = out_drop_flux - out_decay
                        non_decay_drop = non_decay_drop if non_decay_drop > self.flux_error else 0

                        drop_flux += non_decay_drop

                        if self.has_edge(out_v, drop_node):
                            self.remove_edge(out_v, drop_node)


                        #new drop edge
                        self.add_edge(out_v, drop_node,
                            flux_start=drop_flux,
                            flux_end=drop_flux)
                        
                        #new recycling edge
                        self.add_edge(drop_node, self.bulk_node, flux_start=drop_flux, flux_end=drop_flux)

                    #new horizontal edge
                    self.add_edge(in_u, out_v,
                                flux_start=in_flux_start,
                                flux_end=endflux,
                                decay=drop_flux)
                    changed = True

                    break
        
        parentless_nodes = []
        for u, v, data in self.edges(data=True):
            if v == self.bulk_node:
                influx = 0
                for _, _, flux in self.in_edges(u, data='flux_end'):
                    influx += flux
                data['flux_start'] = influx
                data['flux_end'] = influx
            if self.in_degree(u) < 1:
                parentless_nodes.append(u)
        self.remove_nodes_from(parentless_nodes)


    def rein_decay(self, u: RiboNode, v: RiboNode):
        if u.phase == v.phase:
            if u.phase > 0:

                if not self.reinitiation_potential:
                    return 1
                else:
                    half_life = self.reinitiation_potential
            
            else:
                return 0

            
            return 0.5 ** (abs(u.position-v.position) / half_life )
        
        else:
            return 0
        
    def is_retention(self, u, v):
        return u.phase > 0 and v.phase == 0
    
    def is_initiation(self, u, v):
        return u.phase == 0 and v.phase > 0
    
    
    def add_transition(self, source, target, probability):
        """
        Adds new tranistion to graph. 
        """
        self.transitions.add_weighted_edge(source, target, probability)
        self.clear_edges()
        self.construct()

    def add_transitions_from(self, tbunch):
        """
        Adds transitions from an iterable. Must be of form (source, target, weight)
        Flux is recalculated after adding.
        """
        self.transitions.add_weighted_edges_from(tbunch)
        self.clear_edges
        self.construct

    def edge_decay(self, u: RiboNode, v: RiboNode):
        if u.phase == v.phase:
            if u.phase > 0:
                half_life = self.half_life_translation
            elif u.phase == 0:
                half_life = self.half_life_scanning
            else:
                return 0
            if half_life == None:
                return 0
            
            return 1 - (0.5 ** (abs(u.position-v.position) / half_life ))
        
        else:
            return 0
        
    def ternary_complex_proportion(self, u, v) -> float:
        if u.phase == 0 and v.phase == 0 and self.tc_half_life != None:
            half_life = self.tc_half_life
            if half_life:
                return 1-0.5 ** (abs(u.position-v.position) / half_life )
            else: 
                return 1
        else:
            return 1
        

    def _normalize_flux(self):
        flux_keys = ('flux_start', 'flux_end', 'decay')
        fluxes = []
        for u,v, data in self.edges(data=True):
            for key in flux_keys:
                if key in data:
                    fluxes.append(data[key])
        fluxes = set(fluxes)
        max_flux = max(fluxes)
        
        factor = 1/ max_flux if max_flux > 1 else 1
        if factor == 1:
            return
        flux_dict = {}
        for flux in fluxes:
            flux_dict[flux] = flux * factor
        for u,v,data in self.edges(data=True):
            for key in flux_keys:
                if key in data:
                    data[key] = flux_dict[data[key]]
            


    def _is_valid(self):
        self._valid_in_out()

    @property
    def ribopaths(self) -> list:
        """
        Returns a list of the paths with continued 40S association, each as a list of edge tuples.
        """
        paths = []
        for loading in self.successors(self.bulk_node):
            for path in nx.all_simple_edge_paths(self, loading, self.bulk_node):
                paths.append(path)
        return paths
    
    @property
    def translons(self) -> list:
        """
        Returns a list of all translons in the graph (continued 60S association) as a list of edge tuples
        """
        translon_list = []
        for path in self.ribopaths:
            translon = False
            current_translon = []
            for edge in path:

                if translon:
                    if edge[1].phase < 1:
                        translon = False
                        translon_list.append(current_translon)

                    else:
                        current_translon.append(edge)

                elif edge[0].phase > 0:
                    translon = True
                    current_translon.append(edge)

        return translon_list

    def flux_proportion(self, u:RiboNode, v:RiboNode) -> float:
        total_proportion = 0
        for path in nx.all_simple_paths(self, u , v):
            total_proportion += self.flux_proportion_path(path)
        return total_proportion
    

    def flux_proportion_path(self, path: list[RiboNode]) -> float:
        if not path:
            raise ValueError('Cannot compute flux proportion of an empty path')

        # Total flux proportion arriving at path[0] from bulk_node
        entry_proportion = self.node_flux(path[0])

        # Proportion of that flux which travels along the given path
        path_proportion = 1.0
        for i in range(len(path) - 1):
            path_proportion *= self.edge_weight(path[i], path[i + 1])

        return entry_proportion * path_proportion

    def node_flux(self, nbunch: RiboNode) -> float:
        total_flux = 0.0
        for _, _, flux in self.out_edges(nbunch=nbunch, data='flux_start'):
            total_flux += flux
        return total_flux 

    def _valid_in_out(self):
        out_flux = 0
        for _,_, flux in self.in_edges(self.bulk_node, data='flux_end'):
            out_flux += flux

        in_flux = 0
        for node in self.successors(self.bulk_node):
            for _, _, flux in self.out_edges(node, data='flux_start'):
                in_flux += flux

        
        if abs(out_flux - in_flux) > self.flux_error:

            warnings.warn(f'Flux in: {in_flux} does not equal Flux out: {out_flux}, '
                          'this may occur due to accumulated errors in floating point numbers,' 
                          'especially in very complex graphs, ' 
                          'and can be ignored to your deisred level of accuracy.')
            
    def edge_weight(self, u: RiboNode,v: RiboNode):
        return self[u][v]['flux_start'] / self.node_flux(u)