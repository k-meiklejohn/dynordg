from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode, RiboPath
import networkx as nx
import warnings
from copy import deepcopy

class RiboGraphFlux(RiboGraph):

    """
    A RiboGraph that contains information about the available paths of a ribosome \
    aswell as the flux along each path. It is calculated from a TransitionMap, \
    using the from_transition_map() or upon initialisation by massing a TransitionMap instance. 
    Attributes:
    ribopaths: returns a list of lists where each inner list is the path of a 40S subunit through \
    ribosomal phase space where the 40S maintains continuous association with the mRNA

    translons: returns a list of lists where each each inner list is the path of a ribosome through ribosomal phase space \
    where the 60S maintains continuos association with the mRNA.

    Methods:
    flux_proportion(path): returns a float representing the proportion of flux owned by a path. 

    """
    def __init__(self, transition_map: TransitionMap, 
                 incoming_graph_data=None, 
                 half_life_scanning: float|None = None, 
                 half_life_translation: float|None = None, 
                 weight_cutoff=0.0, 
                 reinitiation_half_life = None,
                 ternary_complex_half_life = None,
                 flux_cutoff = 0.001,
                 retention_limit = 1,
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
        
        for u, v in self.transitions.edges:

            if u.phase == -1:
                self.add_edge(self.bulk_node, u)
                flux = flux=self.transitions[u][v]['weight']
                self.add_edge(self.bulk_node, u, weight=flux, flux_start=flux, flux_end=flux)
                self.add_edge(u, v, weight=flux, flux_start=flux, flux_end=flux )

                self._iterate_graph(v, flux) 
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
                    raise RuntimeError(
                        f'Initiation node not found when retaining at node: {u}')

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

            # ── Pass 2: initiation / termination edges ─────────────────────────
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
        out_flux = 0
        for u,v, flux in self.in_edges(self.bulk_node, data='flux_end'):
            out_flux += flux
        print(out_flux)

        
        changed = True
        test_graph = deepcopy(self)
        test_graph.remove_node(self.bulk_node)
    

        while changed:
            changed = False
            for node in list(nx.topological_sort(test_graph)):
                if node.phase == -1:
                    continue

                in_edges = self.in_edges(node)
                #ignore nodes with  more than 1 in edge
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
                out_decay = self[node][out_v]['decay']

                if abs(in_flux_end - out_flux) < self.flux_error:

                    drop_node = RiboNode(node.position, -1)

                    if node.phase == 0 and self.half_life_scanning:
                        
                        if self.has_node(drop_node) and self.in_degree(drop_node) > 1:
                            self.remove_edge(node, drop_node)

                        else:
                            if self.has_node(drop_node):
                                self.remove_node(drop_node)

                    elif node.phase > 0 and self.half_life_translation:
                        if self.has_node(drop_node) and self.in_degree(drop_node) > 1:
                            self.remove_edge(node, drop_node)
                        else:
                            if self.has_node(drop_node):
                                self.remove_node(drop_node)
                    self.remove_node(node)


                    drop_flux = in_flux_start * self.edge_decay(in_u, out_v)
                    endflux = in_flux_start - drop_flux

                    if drop_flux != 0:
                        drop_node = RiboNode(out_v.position, -1)
                        print(out_v, self[out_v][drop_node])
                        print('drop flux:', drop_flux)
                        print('decay:', out_decay)
                        if self.has_node(drop_node) and self.in_degree(drop_node) > 1:
                            self.remove_edge(out_v, drop_node)
                        else:
                            if self.has_node(drop_node):
                                if  self.out_degree(out_v) > 1:
                                    self.remove_node(drop_node)

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
        for u, v, data in self.edges(data=True):
            if v == self.bulk_node:
                influx = 0
                for _, _, flux in self.in_edges(u, data='flux_end'):
                    influx += flux
                data['flux_start'] = influx
                data['flux_end'] = influx

    

                    



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
                paths.append(RiboPath(path))
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
        if path[0] == self.bulk_node:
            entry_proportion = 1.0
        else:
            paths_to_start = list(nx.all_simple_paths(self, self.bulk_node, path[0]))
            if not paths_to_start:
                raise ValueError(f"No path from {self.bulk_node} to {path[0]}")
            entry_proportion = 0.0
            for _path in paths_to_start:
                p = 1.0
                for i in range(len(_path) - 1):
                    p *= self.edge_weight(_path[i], _path[i + 1])
                entry_proportion += p

        # Proportion of that flux which travels along the given path
        path_proportion = 1.0
        for i in range(len(path) - 1):
            path_proportion *= self.edge_weight(path[i], path[i + 1])

        return entry_proportion * path_proportion


    def _valid_in_out(self):
        out_flux = 0
        for u,v, flux in self.in_edges(self.bulk_node, data='flux_end'):
            out_flux += flux

        in_flux = 0
        for node in self.successors(self.bulk_node):
            for u, v, flux in self.out_edges(node, data='flux_start'):
                in_flux += flux

        
        if abs(out_flux - in_flux) > self.flux_error:

            warnings.warn(f'Flux in: {in_flux} does not equal Flux out: {out_flux}, '
                          'this may occur due to accumulated errors in floating point numbers,' 
                          'especially in very complex graphs, ' 
                          'and can be ignored to your deisred level of accuracy.')
            
    def edge_weight(self, u: RiboNode,v: RiboNode):
        total_flux = 0
        for _,_, flux in self.out_edges(u, data='flux_start'):
            total_flux += flux

        return self[u][v]['flux_start'] / total_flux

    
    