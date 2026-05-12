from ..graph import RiboGraph
from .transitionmap import TransitionMap
from .skeleton import RiboSkeleton
from ..core import RiboNode, State
import networkx as nx
import warnings
from copy import deepcopy
from collections import defaultdict
import heapq

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
        (i.e., active translation). Each translon is a list of RiboNodes
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

    def __init__(self, skeleton: RiboSkeleton, 
                 incoming_graph_data=None, 
                 weight_cutoff=0.0, 
                 flux_cutoff = 0.001,
                 retention_limit: int|None = 1,
                 **attr):
        

        super().__init__(incoming_graph_data, **attr)

        self.transitions = skeleton
        self.begun = False
        if incoming_graph_data is not None:
            raise ValueError('Incoming graph data must be left empty, graph is calculated from transition_map')
        self.weight_cutoff = weight_cutoff
        self.flux_cutoff = flux_cutoff
        self.flux_error = 0.000000000000001
        self.retention_limit = retention_limit
        if self.transitions:
            self.construct()   
    
    def construct(self):
        self._iterate_graph_topo(self.transitions.bulk_node, 1)
        self._normalize_flux()
        self._is_valid()

    def _iterate_graph_topo(
        self,
        start_node: RiboNode,
        start_flux: float,
        start_retained: int = 0,
    ):
        self.transitions.remove_node(self.transitions.bulk_node)
        topo_order = list(nx.topological_sort(self.transitions))

        # pending[node][retained] = accumulated flux
        pending: dict[RiboNode, dict[int, float]] = defaultdict(lambda: defaultdict(float))
        pending[start_node][start_retained] += start_flux

        # Only process a node once all upstream flux has arrived
        in_degree_remaining = {node:self.transitions.in_degree(node) for node in self.transitions.nodes}
        in_degree_remaining[start_node] = 0

        def dispatch(target: RiboNode, flux: float, retained: int):
            """Accumulate flux into target and mark one upstream sender as done."""
            if target == self.bulk_node:
                return
            pending[target][retained] += flux
            in_degree_remaining[target] -= 1

        for node in topo_order:
            if in_degree_remaining.get(node, 0) > 0 or not pending.get(node):
                continue

            for retained, flux in list(pending[node].items()):
                if flux < self.flux_error:
                    continue

                for u, v, w in self.transitions.out_edges(node, data="weight"):
                    new_flux = flux * w

                    # Below cutoff — drop the ribosome if translating
                    if v.phase != -1 and new_flux < self.flux_cutoff:
                        if node.phase >= 1:
                            drop_node = RiboNode(u.position, state=State(-1))
                            self.add_edge(u, drop_node, flux=new_flux)
                            self.add_edge(drop_node, self.bulk_node, flux=new_flux)
                        continue

                    next_retained = retained + 1 if self.is_retention(u, v) else retained
                    self.add_edge(u, v, flux=new_flux)

                    if v.phase == -1:
                        self.add_edge(v, self.bulk_node, flux=new_flux)
                        continue

                    dispatch(v, new_flux, next_retained)

            del pending[node]



        
    def collapse_unused_nodes(self):
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



        
    def is_retention(self, u, v):
        return u.phase > 0 and v.phase == 0
    
    def is_initiation(self, u, v):
        return u.phase == 0 and v.phase > 0
    
    
    def add_transition(self, source, target, weight):
        """
        Adds new tranistion to graph. 
        """
        self.transitions.add_edge(source, target, weight=weight)
        self.clear_edges()
        self.construct()

    def add_transitions_from(self, tbunch):
        """
        Adds transitions from an iterable. Must be of form (source, target, weight)
        Flux is recalculated after adding.
        """
        self.transitions.add_weighted_edges_from(tbunch)
        self.clear_edges()
        self.construct()



    def _normalize_flux(self):
        flux_keys = ('flux')
        fluxes = []
        for u,v, data in self.edges(data=True):
            if 'flux' in data:
                fluxes.append(data['flux'])
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
    def ribopaths(self) -> list[list[RiboNode]]:
        """
        Returns a list of the paths with continued 40S association, each as a list of edge tuples.
        """
        paths = []
        for loading in self.successors(self.bulk_node):
            for path in nx.all_simple_paths(self, loading, self.bulk_node):
                paths.append(path)
        return paths
    
    @property
    def translons(self) -> list[list[RiboNode]]:
        """
        Returns a list of all translons in the graph (continued 60S association) as a list of edge tuples
        """
        translon_list = []
        for path in self.ribopaths:
            translon = False
            current_translon = []
            for node in path:

                if translon:
                    if node.phase < 1:
                        translon = False
                        translon_list.append(current_translon)

                    else:
                        current_translon.append(node)

                elif node.phase > 0:
                    translon = True
                    current_translon.append(node)

        return list(set(translon_list))

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
        for _, _, flux in self.out_edges(nbunch=nbunch, data='flux'):
            total_flux += flux
        return total_flux 

    def _valid_in_out(self):
        out_flux = 0
        for _,_, flux in self.in_edges(self.bulk_node, data='flux'):
            out_flux += flux

        in_flux = 0
        for node in self.successors(self.bulk_node):
            for _, _, flux in self.out_edges(node, data='flux'):
                in_flux += flux

        
        if abs(out_flux - in_flux) > self.flux_error:

            warnings.warn(f'Flux in: {in_flux} does not equal Flux out: {out_flux}, '
                          'this may occur due to accumulated errors in floating point numbers,' 
                          'especially in very complex graphs, ' 
                          'and can be ignored to your deisred level of accuracy.')
            
    def edge_weight(self, u: RiboNode,v: RiboNode):
        return self[u][v]['flux'] / self.node_flux(u)
    
    @property
    def simple(self):
        out = deepcopy(self)
        out.clear()
        out.bulk_node = self.bulk_node.simple
        for u, v, flux in self.edges(data='flux'):

            if u.simple == v.simple:
                continue
            elif v.phase == -1 and u.position != v.position and u.phase != -1:
                out.add_edge(u.simple, RiboNode(v.position, State(u.phase)), flux_start=flux, flux_end=0, decay=flux)
                out.add_edge(RiboNode(v.position, State(u.phase)), v.simple, flux_start=flux, flux_end=flux)
    
            else:
                out.add_edge(u.simple, v.simple,
                            flux_start=flux,
                            flux_end=flux)
        return out
    
    def prune_recycle_edges(self) -> None:
        """
        Remove bulk-to-bulk edges and any nodes they leave isolated.

        Edges where both u and v have phase == -1 are recycling arcs internal
        to the bulk pool.  They carry no mRNA-positional information and
        produce degenerate geometry, so they are stripped before layout runs.
        """
        dead = []

        dead = [(u, v) for u, v in self.edges
                if u.phase == -1 and v.phase == -1]
        self.remove_edges_from(dead)
        isolated = [n for n, d in self.degree() if d < 1]
        self.remove_nodes_from(isolated)
