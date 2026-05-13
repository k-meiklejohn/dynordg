from ..graph.ribograph import RiboGraph
from ..core import Transition, RiboNode

class TransitionMap(RiboGraph):
    """
    A validated directed graph representing the transition probabilities between
    ribosomal phase-space nodes along an mRNA transcript.

    TransitionMap stores the discrete probabilistic steps a ribosome can take —
    scanning, initiating, elongating, terminating, or reinitiating — as weighted
    edges between RiboNodes. It enforces correctness at the edge level and serves
    as the primary input to RiboGraphFlux for flux propagation.

    Edge weights represent the probability of a ribosome at node u transitioning
    to node v. For any non-bulk node with outgoing edges, weights should sum to
    no greater than 1.0; the remainder represents flux lost to drop-off, which
    RiboGraphFlux handles implicitly via its decay parameters.

    Construction
    ------------
    Edges must be added via add_weighted_edge() or add_weighted_edges_from()
    rather than the standard NetworkX add_edge()/add_edges_from(), which are
    disabled. This ensures all edges are coerced through RiboTransition
    validation before being added to the graph.

    Parameters
    ----------
    incoming_graph_data : None
        Reserved for NetworkX compatibility. Must be left as None; the graph
        is populated exclusively through the weighted-edge API.

    Methods
    -------
    add_weighted_edge(source, target, probability)
        Add a single transition. Accepts either three positional arguments
        (source, target, probability), a length-3 tuple, or a RiboTransition
        instance. Validates the edge via RiboTransition before adding.
    add_weighted_edges_from(ebunch_to_add)
        Add multiple transitions from an iterable of RiboTransition instances
        or (source, target, weight) tuples.
    to_fluxgraph(half_life_translation=None, half_life_scanning=None)
        Construct and return a RiboGraphFlux from this TransitionMap, optionally
        specifying ribosome half-lives to control decay during translation and
        scanning respectively.

    Notes
    -----
    - Phase conventions mirror those in RiboGraphFlux: phase == -1 is the bulk
      cytoplasmic pool, phase == 0 is a scanning 40S, phase > 0 is an elongating
      80S ribosome.
    - Validation is run after every edge addition. Nodes with no outgoing edges
      (terminal nodes) and bulk nodes (phase == -1) are exempt from weight checks.
    - add_edge() and add_edges_from() raise NotImplementedError to prevent
      unvalidated edges from being inserted directly.

    Examples
    --------
    >>> tm = TransitionMap()
    >>> tm.add_weighted_edge(RiboNode(0, -1), RiboNode(10, 0), 1.0)
    >>> tm.add_weighted_edges_from([
    ...     (RiboNode(10, 0), RiboNode(10, 1), 0.3),
    ...     (RiboNode(10, 0), RiboNode(20, 0), 0.7),
    ... ])
    >>> flux = tm.to_fluxgraph(half_life_scanning=50.0)
    """

    def __init__(self, incoming_graph_data=None, **attr):
        super().__init__(incoming_graph_data, **attr)
        self._is_valid()
    
    def add_transition(self, transition: Transition):
        if isinstance(transition, Transition):
            super().add_edge(transition.source, transition.target, weight=transition.weight)
        else:
            raise TypeError(f"add_transition requires type 'Transition', got {type(transition).__name__}")

    def add_transitions_from(self, transitions):
        for t in transitions:
            self.add_transition(t)



    def _is_valid_weight(self):

        for node in self.nodes:

            if node.phase == -1:
                continue

            if not any(True for _ in self.successors(node)):
                continue

            total_weight = 0

            for (u, v, w) in self.out_edges(node, data='weight'):
                total_weight += w

            # # if total_weight <= 0 or total_weight > 1:
            #     raise ValueError(f'Total weight of edges from a single node (not from phase=-1) must be greater than 0 and less than or equal to 1\n' \
            #     f'Offending Node: {node}, Weight: {total_weight}')
            
    def _is_valid(self):
        self._is_valid_weight()

    def downstream_node(self, node: RiboNode) -> RiboNode | None:
        """Return the nearest downstream node in the same phase whose factors are a subset of node's."""
        positions = sorted(
            p.position for p in self.nodes
            if p.position > node.position and p.phase == node.phase
        )

        for pos in positions:
            candidates: list[RiboNode] = [n for n in self.nodes if n.position == pos and n.phase == node.phase]
            options = [
                opt for opt in candidates
                if not opt.factors or all([f in node.factors for f in opt.factors])
            ]

            if options:
                return max(options)

        return None

    # def downstream_node(self, node: RiboNode) -> RiboNode|None:
        
    #     same_phase_ahead = sorted([p.position for p in self.nodes 
    #                         if p.position > node.position 
    #                         and p.phase == node.phase])
        
    #     if not same_phase_ahead:
    #         return None
    #     print('==============DOWNSTREAM NODE===============')
    #     print("node:",node)
    #     print("ahead", same_phase_ahead)
    #     next_pos = same_phase_ahead.pop(0)
    #     next_nodes: list[RiboNode] = [n for n in self.nodes
    #                   if n.position == next_pos
    #                   and n.phase == node.phase]
    #     print("NEXT NODES", next_nodes)
    #     next_node = None
    #     options: list[RiboNode] = []
    #     while not next_node:
    #         for opt in next_nodes:
    #             print("----trying-------------")
    #             print('option:',opt)
    #             if len(opt.factors) and len(node.factors):
    #                 if all([f in node.factors for f in opt.factors]):
    #                     options.append(opt)
    #             elif not len(opt.factors):
    #                 options.append(opt)

    #                 continue
    #         if len(options):
    #             next_node = max(options)
    #             break
    #         if not len(same_phase_ahead):
    #             return None
    #         next_pos = same_phase_ahead.pop(0)
    #         next_nodes: list[RiboNode] = [n for n in self.nodes
    #                     if n.position == next_pos
    #                     and n.phase == node.phase]

    #     print('out:', next_node)
    #     return next_node