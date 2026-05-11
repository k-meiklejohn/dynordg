from ..graph import RiboGraph
from .transitionmap import TransitionMap
from ..core import RiboNode, Transition, State

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
        next_node = {n: self.transitions.downstream_node(n) for n in self.transitions.nodes if n.phase != -1}
        print("-----------------\nDICT of NEXT NODES")
        print(next_node)
        print("-----------------\n")
        cont_weight = {}
        for node in self.transitions.nodes:
            remaining_weight= 1
            for _,_,w in self.transitions.out_edges(node, data='weight'):
                remaining_weight -= w
            cont_weight[node] = remaining_weight

        print("continue weights:\n", cont_weight)

        for u,v,w in self.transitions.edges(data='weight'):
            print('==========================')
            print('u:', u,'v:', v)
            self.add_transition(Transition(u,v,w))
            if u.phase != -1:
                tnode = next_node[u]
                print('after u:', tnode)
                if tnode:
                    if tnode.state == u.state:
                        for weight, target in self.half_life_swap(u, tnode, cont_weight[u]):
                            if weight != 0:
                                self.add_transition(Transition(u, target, weight))
                    else:
                        same_state = RiboNode(tnode.position, state=u.state)
                        for weight, target in self.half_life_swap(u, same_state, cont_weight[u]):
                            if weight != 0:
                                self.add_transition(Transition(u, target, weight))
                        if same_state in self:
                            self.add_transition(Transition(same_state, tnode, 1- cont_weight[tnode]))
                            cont_weight[same_state] = cont_weight[tnode]
                else:
                    drop_node = RiboNode(u.position, State(-1))
                    self.add_transition(Transition(u, drop_node, cont_weight[u]))
                    self.add_transition(Transition(drop_node, self.bulk_node, 1))

            else:
                self.add_transition(Transition(self.bulk_node, u, 1))
            if v.phase != -1:
                tnode = next_node[v]
                print('after v:', tnode)
                if tnode:
                    if tnode.state == v.state:
                        for weight, target in self.half_life_swap(v, tnode, cont_weight[v]):
                            if weight != 0:
                                self.add_transition(Transition(v, target, weight))
                    else:
                        same_state = RiboNode(tnode.position, state=v.state)
                        for weight, target in self.half_life_swap(v, tnode, cont_weight[v]):
                            if weight != 0:
                                self.add_transition(Transition(v, target, weight))
                        if same_state in self:
                            self.add_transition(Transition(same_state, tnode, 1- cont_weight[tnode]))
                            cont_weight[same_state] = cont_weight[tnode]
                else:
                    drop_node = RiboNode(v.position, State(-1))
                    self.add_transition(Transition(v, drop_node, cont_weight[v]))
                    self.add_transition(Transition(drop_node, self.bulk_node, 1))
            else:
                self.add_transition(Transition(v, self.bulk_node, 1))

        self.remove_node(self.bulk_node)

    def half_life_swap(self, u: RiboNode, v: RiboNode, cont_weight:float) -> list[tuple[float, RiboNode]] :
        swaps: list[tuple[float, RiboNode]] = []
        if u.phase != v.phase:
            raise ValueError('Decay must be in same phase')
        if u.phase != -1 and v.phase != -1:
            swaps.append((cont_weight*self.edge_decay(u,v), RiboNode(v.position, State(-1))))
            swaps.append((cont_weight - cont_weight*self.edge_decay(u,v), v))
        if u.phase == 0 and 'ternary_complex' not in u.factors:
            swaps.append((self.ternary_complex_proportion(u, v), RiboNode(v.position, v.state + 'ternary_complex')))
        elif u.phase > 0 and 'scanning_factors' in u.factors:
            swaps.append((self.scanning_factor_proportion_lost(u,v), RiboNode(v.position, v.state - 'scanning_factors')))

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
