from typing import Literal
from dataclasses import dataclass, field
from ...core import RiboNode

# ── Type aliases ─────────────────────────────────────────────────────────────

Pt   = tuple[float, float]   # (x, y) coordinate pair
Edge = tuple[RiboNode, RiboNode]


# ─────────────────────────────────────────────────────────────────────────────
# Phase 1 output  –  EdgeSpec
# ─────────────────────────────────────────────────────────────────────────────

EdgeType = Literal[
    'frameshift', '40s_retention', 'drop', 'initiation', 'load',
    '0', '1', '2', '3',
]
"""
Semantic classification of an edge in the flux graph.

    '0', '1', '2', '3'  –  horizontal continuation edges within a phase lane
    'frameshift'              –  phase-change edge where the ribosome moves forward
                            in position (codon displacement given by shift_n)
    'initiation'         –  60S joining; 40S (phase 0) → 80S (phase > 0)
    '40s_retention'      –  80S → 40S retention without forward movement
    'drop'               –  ribosome dissociates from mRNA; enters bulk pool
    'load'               –  ribosome is recruited from the bulk pool onto mRNA
"""


@dataclass
class EdgeSpec:
    """
    Everything known about an edge before any coordinates are assigned (Phase 1).

    EdgeSpec is the primary output of the classification pass. It captures the
    flux carried by the edge and its semantic role in the ribosomal lifecycle,
    but holds no geometric data — coordinate assignment happens in later phases.

    Attributes
    ----------
    u : RiboNode
        Source node.
    v : RiboNode
        Target node.
    flux_start : float
        Flux entering the edge at u (before any decay along the edge).
    flux_end : float
        Flux arriving at v (after decay; equals flux_start for lossless edges).
    etype : EdgeType
        Semantic type of the edge; controls lane routing and visual style.
    direction : int or None
        Resolved travel direction along the mRNA: +1 (5'→3'), -1 (3'→5'),
        0 (vertical / phase-change only), or None for bulk edges where
        direction is determined later.
    shift_n : int
        Codon displacement for 'shift' edges. Zero for all other edge types.
    is_event : bool
        True when u.phase != v.phase, i.e. the edge represents a phase
        transition (initiation, termination, or retention) rather than
        simple progression within a lane.
    """
    u:          RiboNode
    v:          RiboNode
    flux_start: float
    flux_end:   float
    etype:      EdgeType
    direction:  int | None
    shift_n:    int  = 0
    is_event:   bool = False


# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 output  –  NodeLayout
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EdgeSlot:
    """
    A single entry in the ordered slot list on one side of a node (Phase 2).

    Slots define the vertical stacking order of edges at a node face, which
    determines how edge bands are packed without overlap during geometry
    assignment. Each slot records both the edge identity and the direction
    in which flux flows at this particular node.

    Attributes
    ----------
    edge : Edge
        The (u, v) pair identifying the edge this slot belongs to.
    direction : int
        Resolved travel direction at this node face: +1, -1, or 0.
    """
    edge:      Edge
    direction: int


@dataclass
class NodeLayout:
    """
    Ordered in/out edge slots for a single node, ready for geometry (Phase 2).

    NodeLayout is produced by the slot-assignment pass and consumed by the
    geometry pass (Phase 3). It encodes which edges enter and leave the node
    and in what order they should be stacked vertically on each face, along
    with a hint for the direction of any associated drop edge.

    Attributes
    ----------
    node : RiboNode
        The node this layout describes.
    in_slots : list of EdgeSlot
        Ordered slots for edges arriving at this node, bottom-to-top.
    out_slots : list of EdgeSlot
        Ordered slots for edges leaving this node, bottom-to-top.
    drop_direction : int or None
        Direction (+1 or -1) of the drop edge leaving this node, if one
        exists. None when the node has no associated drop edge.
    """
    node:           RiboNode
    in_slots:       list[EdgeSlot] = field(default_factory=list)
    out_slots:      list[EdgeSlot] = field(default_factory=list)
    drop_direction: int | None     = None


# ─────────────────────────────────────────────────────────────────────────────
# Phase 3 output  –  EdgeGeom
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EdgeGeom:
    """
    Full coordinate data for a single edge, progressing through Phases 3 and 4.

    After Phase 3, all point coordinates are in *local, node-relative y space*:
    x is in world units but y is expressed relative to the node's own vertical
    origin. After Phase 4 (world-space alignment), all coordinates are absolute
    and ready for the renderer.

    Primary rectangle
    -----------------
    The four corners (in0, in1, out0, out1) define the main flux band as a
    quadrilateral. in0/in1 are the bottom and top corners at the source node
    face; out0/out1 are the bottom and top corners at the target node face.

    Curve taper hints
    -----------------
    in_quad and out_quad encode which quadrant the band enters/exits at each
    node face, allowing the renderer to apply the correct Bézier taper.

    Bulk arrow geometry
    -------------------
    out_extent, in_extent, and the four *_base points define the arrowhead
    geometry for load and drop edges that connect to the bulk pool.

    Decay triangle
    --------------
    decay0, decay1, decay2 are the three vertices of the triangle drawn
    where flux is lost to ribosome drop-off along a scanning or translation
    edge. Only populated when flux_start > flux_end.

    Alignment anchors
    -----------------
    in_bot and out_bot are temporary anchors used during Phase 4 to align
    this edge's local y coordinates into world space. They are not used by
    the renderer and may be treated as internal after Phase 4.

    Helper rectangles
    -----------------
    in_helper_rects and out_helper_rects store auxiliary rectangles used to
    compute vertical offsets for stacked bulk edges. Accessible together via
    the helper_rects property.

    Attributes
    ----------
    edge : Edge
        The (u, v) pair this geometry belongs to.
    in0, in1 : Pt or None
        Bottom and top corners of the band at the source node face.
    out0, out1 : Pt or None
        Bottom and top corners of the band at the target node face.
    in_quad, out_quad : int or None
        Quadrant hints for Bézier taper at the source and target faces.
    out_extent, in_extent : Pt or None
        Arrow tip positions for bulk (load/drop) edges.
    out_left_base, out_right_base : Pt or None
        Arrowhead base corners at the target end of a bulk edge.
    in_left_base, in_right_base : Pt or None
        Arrowhead base corners at the source end of a bulk edge.
    decay0, decay1, decay2 : Pt or None
        Vertices of the decay-loss triangle. None when edge is lossless.
    in_bot, out_bot : Pt or None
        Phase 4 alignment anchors; not used by the renderer.
    in_helper_rects : list of list of Pt
        Auxiliary stacking rectangles at the source end.
    out_helper_rects : list of list of Pt
        Auxiliary stacking rectangles at the target end.
    etype : EdgeType
        Semantic edge type, carried from EdgeSpec for renderer styling.
    is_event : bool
        True when the edge represents a phase transition; carried from EdgeSpec.
    flux_end : float
        Flux arriving at v; carried from EdgeSpec for renderer use.
    """
    edge: Edge

    # Primary rectangle corners
    in0:  Pt | None = None
    in1:  Pt | None = None
    out0: Pt | None = None
    out1: Pt | None = None

    # Curve taper quadrant hints
    in_quad:  int | None = None
    out_quad: int | None = None

    # Bulk arrow geometry
    out_extent:      Pt | None = None
    in_extent:       Pt | None = None
    out_left_base:   Pt | None = None
    out_right_base:  Pt | None = None
    in_left_base:    Pt | None = None
    in_right_base:   Pt | None = None

    # Decay triangle
    decay0: Pt | None = None
    decay1: Pt | None = None
    decay2: Pt | None = None

    # Phase 4 alignment anchors (internal; stripped after use)
    in_bot:  Pt | None = None
    out_bot: Pt | None = None

    # Stacked bulk offset helpers
    in_helper_rects:  list[list[Pt]] = field(default_factory=list)
    out_helper_rects: list[list[Pt]] = field(default_factory=list)

    @property
    def helper_rects(self):
        return self.in_helper_rects + self.out_helper_rects

    # Carried from EdgeSpec for renderer use
    etype:    EdgeType = '0'
    is_event: bool     = False
    flux_end: float    = 0.0