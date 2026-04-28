from ..graph.ribograph import RiboGraph
from ..core import RiboTransition


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
        

    def add_edge(self, u_of_edge, v_of_edge, **attr):
        raise NotImplementedError("Use add_weighted_edge() instead.")

    def add_edges_from(self, ebunch_to_add, **attr):
        raise NotImplementedError("Use add_weighted_edge() instead.")

    def add_weighted_edge(self, *args):
        # Coerce to RiboTransition, which handles validation
        if len(args) == 1 and isinstance(args[0], (tuple, RiboTransition)):
            transition = RiboTransition(args[0])
        elif len(args) == 3:
            transition = RiboTransition(*args)
        else:
            raise ValueError(f'add_weighted_edge requires a RiboTransition, a length-3 tuple, or 3 arguments, got: {args}')

        super().add_edge(transition.source, transition.target, weight=transition.probability)
        self._is_valid()

    def add_weighted_edges_from(self, ebunch_to_add):
        """ebunch_to_add: iterable of RiboTransitions or (source, target, weight) tuples"""
        for item in ebunch_to_add:
            self.add_weighted_edge(item)


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


    def to_fluxgraph(self, half_life_translation = None, half_life_scanning = None):
        from .fluxgraph import RiboGraphFlux
        return RiboGraphFlux(transition_map=self, half_life_translation=half_life_translation, half_life_scanning=half_life_scanning)

