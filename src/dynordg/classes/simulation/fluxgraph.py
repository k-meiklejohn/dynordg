from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode, RiboPath
import networkx as nx
import warnings

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
    def __init__(self, transition_map: TransitionMap, incoming_graph_data=None, half_life_scanning: float|None = None, half_life_translation: float|None = None, cutoff=0.0, **attr):
        super().__init__(incoming_graph_data, **attr)
        self.transitions = transition_map
        self.begun = False
        if incoming_graph_data is not None:
            raise ValueError('Incoming graph data must be left empty, graph is calculated from transition_map')
        self.half_life_translation = half_life_translation
        self.half_life_scanning = half_life_scanning
        if map:
            self.construct(cutoff=cutoff)   

    @classmethod
    def from_transition_map(cls, transition_map: TransitionMap, half_life_translation:float|None =None, half_life_scanning:float|None =None):
        return cls(transition_map=transition_map, half_life_scanning=half_life_scanning, half_life_translation=half_life_translation)
    
    def construct(self, cutoff=0.0):
        for u, v in self.transitions.edges:
            if u.phase == -1:
                self.add_edge(self.bulk_node, u)
                flux = flux=self.transitions[u][v]['weight']
                self.add_edge(u, v, weight=flux, flux_start=flux, flux_end=flux )

                self._iterate_graph((v), flux, cutoff=cutoff) 

        self._is_valid()


    def _downstream_node(self, node: RiboNode):
        if not any([p.position > node.position 
                            and p.phase == node.phase
                            for p in self.transitions.nodes ]):
            return self.bulk_node
            
        return RiboNode(( min([p.position
                            for p in self.transitions.nodes 
                            if p.position > node.position 
                            and p.phase == node.phase]),

                            node.phase))
    
    def _iterate_graph(self, node: RiboNode, flux, weight=1, cutoff=0.0):

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
            self.add_edge(next_node, drop_node, flux_start=drop_flux, flux_end=drop_flux, weight=weight*drop_flux/flux) # this is the drop edge
            self.add_edge(drop_node, self.bulk_node, flux_start=drop_flux, flux_end=drop_flux, weight=1) #returning to the bulk
        self.add_edge(node, next_node, flux_start=flux, flux_end = endflux, weight=weight*endflux/flux, decay=drop_flux) # this is the horizontal edge

        #### calculate flux for each edge off next node ####

        remaining_weight = 1
        for u, v, w in self.transitions.out_edges(next_node, data='weight'):

            new_flux = endflux * w
            if new_flux < cutoff:
                continue
            remaining_weight -= w

            #adds the edge corresponding to the events defined by the transitions from this node
            self.add_edge(u, v, flux_start=new_flux, flux_end=new_flux, weight=w) 

            if v.phase == -1:
                #Adds the edge from the from the transcript bulk node to the general bulk node
                self.add_edge(v, self.bulk_node, flux_start = new_flux, flux_end=new_flux, weight=w)    
                continue

            self._iterate_graph(v, new_flux, w, cutoff=cutoff)

        #### Continue graph on same phase if weight remaining ####
        if remaining_weight == 0:
            return
        else:
            self._iterate_graph(next_node, endflux*remaining_weight, weight=remaining_weight, cutoff=cutoff)

    
    
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
        if u.phase == v.phase or (u.phase > 0 and v.phase > 0):
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

    def flux_proportion(self, path: list[tuple[RiboNode,RiboNode]]) -> float|None:
        if path is not None:
            weights = nx.get_edge_attributes(self, name='weight')
            paths_to_start = nx.all_simple_edge_paths(self, self.bulk_node, path[0][0])
            proportion = 1.0


            if paths_to_start is None:
                raise ValueError(f"No path exists from {self.bulk_node} to {self[0][0]}")
            
            start_proportion = 0
            for _path in paths_to_start:
                for edge in _path:
                    proportion *= weights[edge]
                start_proportion += proportion


            proportion = start_proportion
            if path[-1][1] != self.bulk_node:
                for finish_edge in self.out_edges(path[-1][1]):
                    if finish_edge[1].phase < 1:
                        new_path = path.copy()
                        for edge in new_path:
                            proportion *= weights[edge]
                
            else:
                for edge in path:
                    proportion *= weights[edge]
        else:
            raise ValueError('No weight associated with empty path')

        return proportion

    def _valid_in_out(self):
        out_flux = 0
        for u,v, flux in self.in_edges(self.bulk_node, data='flux_end'):
            out_flux += flux

        in_flux = 0
        for node in self.successors(self.bulk_node):
            for u, v, flux in self.out_edges(node, data='flux_start'):
                in_flux += flux

        if out_flux != in_flux:
            return
            raise RuntimeError(f'Flux in: {in_flux} does not equal Flux out: {out_flux}')
    
    