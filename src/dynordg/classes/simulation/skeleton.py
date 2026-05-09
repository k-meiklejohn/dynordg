from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode

class RiboSkeleton(TransitionMap):
    def __init__(self, transition_map:TransitionMap,
                 incoming_graph_data=None, 
                half_life_scanning: float|None = None, 
                half_life_translation: float|None = None, 
                reinitiation_half_life = None,
                ternary_complex_half_life = None,
                **attr,):
        super().__init__(incoming_graph_data, **attr) 
        self.half_life_translation = half_life_translation
        self.half_life_scanning = half_life_scanning
        self.reinitiation_potential=reinitiation_half_life
        self.tc_half_life = ternary_complex_half_life
        self.transitions = transition_map
        self._construct()

    def _construct(self):
        to_add = []
        for node in self.transitions.nodes:
            if node.phase == -1:
                continue
            opp_node = RiboNode(node.position, node.phase, not node.factors)
            to_add.append(opp_node)

        self.transitions.add_nodes_from(to_add)


        next_node = {n: self.transitions.downstream_node(n) for n in self.transitions.nodes if n.phase != -1}
        remaining_weights = {}
        for node in self.transitions.nodes:
            remaining_weight= 1
            for _,_,w in self.transitions.out_edges(node, data='weight'):
                remaining_weight -= w
            remaining_weights[node] = remaining_weight

        for u,v,w in self.transitions.edges(data='weight'):
            self.add_weighted_edge(u,v, w)
            if u.phase != -1:
                if next_node[u]:
                    for weight, target in self.half_life_swap(u, next_node[u], remaining_weights[u]):
                        if weight != 0:
                            self.add_weighted_edge(u, target, weight)
            else:
                self.add_weighted_edge(self.bulk_node, u, 1)
            if v.phase != -1:
                if next_node[v] != None:
                    for weight, target in self.half_life_swap(v, next_node[v], remaining_weights[v]):
                        if weight != 0:
                            self.add_weighted_edge(v, target, weight)
            else:
                self.add_weighted_edge(v, self.bulk_node, 1)

        self.remove_all_indegree_zero_nodes()
        self.remove_node(self.bulk_node)

    def half_life_swap(self, u: RiboNode, v: RiboNode, cont_weight:float) -> list[tuple[float, RiboNode]] :
        swaps: list[tuple[float, RiboNode]] = []
        if u.phase != v.phase:
            print(u,v)
            raise ValueError('Decay must be in same phase')
        if u.phase != -1 and v.phase != -1:
            swaps.append((cont_weight*self.edge_decay(u,v), RiboNode(v.position, -1, False)))
            print(u,v, self.edge_decay(u,v))
            swaps.append((cont_weight - cont_weight*self.edge_decay(u,v), v))
        if u.phase == 0 and not u.factors:
            swaps.append((self.ternary_complex_proportion(u, v), RiboNode(v.position, v.phase, True)))
        elif u.phase > 0:
            swaps.append((self.scanning_factor_proportion_lost(u,v), RiboNode(v.position, v.phase, False)))

        return swaps

    def scanning_factor_proportion_lost(self, u: RiboNode, v: RiboNode):
        if u.phase == v.phase:
            if u.phase > 0:

                if not self.reinitiation_potential:
                    return 0
                else:
                    half_life = self.reinitiation_potential
            
            else:
                return 0

            
            return 1 - 0.5 ** (abs(u.position-v.position) / half_life )
        
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
