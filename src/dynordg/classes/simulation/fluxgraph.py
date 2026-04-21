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

                self._iterate_graph((v), flux) 
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
    
    def _iterate_graph(self, node: RiboNode,
                       flux,
                       weight=1,
                       retained=0,
                       initiation_node: RiboNode|None = None,
                       retention_node: RiboNode|None = None):

        if node==self.bulk_node:
            return
        
        next_node = self._downstream_node(node)
        if next_node is None:
            return 
    
        #### Calculate decay of ribosomes based on half life ####

        endflux = flux * self.edge_decay(node, next_node)
        drop_flux = flux - endflux
        if drop_flux != 0:
            drop_node = RiboNode(next_node.position, -1)
            self.add_edge(next_node, drop_node, flux_start=drop_flux, flux_end=drop_flux)
            self.add_edge(drop_node, self.bulk_node, flux_start=drop_flux, flux_end=drop_flux)

        self.add_edge(node, next_node, flux_start=flux, flux_end = endflux, weight=weight*endflux/flux, decay=drop_flux) # this is the horizontal edge

        ### Calculate reiniation potential of ribosomes
  


        #### calculate flux for each edge off next node ####

        remaining_weight = 1
        abs_endflux = endflux
        for u, v, w in self.transitions.out_edges(next_node, data='weight'):


            if self.is_retention(u,v):
                if initiation_node != None:
                    retention_flux = endflux * self.rein_decay(initiation_node, next_node)
                else: 
                    raise RuntimeError(f'Initation node not found when retaining at node: {u}')
                if retention_flux < self.flux_cutoff:
                    continue
                if self.retention_limit != None:
                    if retained >= self.retention_limit:
                        continue

                
                new_flux = retention_flux * w
                w = new_flux/abs_endflux
                endflux -= new_flux
            else:
                continue
            if new_flux < self.flux_cutoff:
                continue

            retention_node = v
         
            remaining_weight -= w

            #adds the edge corresponding to the events defined by the transitions from this node
            self.add_edge(u, v, flux_start=new_flux, flux_end=new_flux) 
            self._iterate_graph(v, new_flux, w, retained+1, retention_node=retention_node)

        for u, v, w in self.transitions.out_edges(next_node, data='weight'):
  
            if self.is_retention(u,v):
                continue
            elif self.is_initiation(u,v):
                initiation_node = v
                if retained:
                    w *= self.ternary_complex_proportion(retention_node, next_node)

            new_flux = endflux * w
            if new_flux < self.flux_cutoff and v.phase != -1:
                    continue
            remaining_weight -= new_flux/abs_endflux

            #adds the edge corresponding to the events defined by the transitions from this node
            self.add_edge(u, v, flux_start=new_flux, flux_end=new_flux)


            if v.phase == -1:
                self.add_edge(v, self.bulk_node, flux_start =new_flux, flux_end=new_flux)
                continue

            self._iterate_graph(v, new_flux, w, retained=retained, initiation_node=initiation_node, retention_node=retention_node)

        #### Continue graph on same phase if weight remaining #####
        if remaining_weight < -self.flux_error:
            raise ValueError(f'Weight from node: {next_node} exceeds 1: {abs(remaining_weight - 1 )}')
        
        if remaining_weight < self.flux_error:
            return
        else:
            self._iterate_graph(next_node,
                                endflux*remaining_weight,
                                weight=remaining_weight,
                                retained=retained,
                                initiation_node=initiation_node,
                                retention_node=retention_node
                                )

    def _collapse_unused_nodes(self):
        
        changed = True
        test_graph = deepcopy(self)
        test_graph.remove_node(self.bulk_node)
    

        while changed:
            for node in list(nx.topological_sort(test_graph)):
                changed = False
                if node.phase == -1:
                    continue

                in_edges = self.in_edges(node)
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


                    endflux = in_flux_start * self.edge_decay(in_u, out_v)
                    drop_flux = in_flux_start - endflux

                    if drop_flux != 0:
                        drop_node = RiboNode(out_v.position, -1)
                        if self.has_node(drop_node) and self.in_degree(drop_node) > 1:
                            self.remove_edge(out_v, drop_node)
                        else:
                            if self.has_node(drop_node):
                                if  self.out_degree(out_v) > 1:
                                    self.remove_node(drop_node)

                        #new drop edge
                        self.add_edge(out_v, drop_node,
                            flux_start=drop_flux,
                            flux_end=drop_flux,
                            weight=drop_flux / in_flux_start)
                        
                        #new recycling edge
                        self.add_edge(drop_node, self.bulk_node, flux_start=drop_flux, flux_end=drop_flux, weight=1)

                    #new horizontal edge
                    self.add_edge(in_u, out_v,
                                flux_start=in_flux_start,
                                flux_end=endflux,
                                weight=endflux / in_flux_start, 
                                decay=drop_flux) # this is the horizontal edge
                    changed = True
                    break
        for u, v, data in self.edges(data=True):
            if u.phase == -1 and v == self.bulk_node:
                influx = 0
                for _, _, flux in self.in_edges(u, data='flux_end'):
                    influx += flux
                data['flux_end'] = influx
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
                return 1
            if half_life == None:
                return 1
            
            return 0.5 ** (abs(u.position-v.position) / half_life )
        
        else:
            return 1
        
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

    
    