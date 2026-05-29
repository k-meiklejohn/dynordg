from __future__ import annotations
from .core import RiboNode
from dataclasses import dataclass
from networkx import topological_sort
import math
from typing import Literal
from matplotlib.pyplot import show
from matplotlib.figure import Figure
from dataclasses import field
from .graph import RiboGraph, SimpleFluxGraph, FluxGraph




ARROW_DEPTH_SCALE: float = 0.6
_GEOM_POINT_KEYS = ('out0', 'in0', 'out1', 'in1',
                        'decay0', 'decay1', 'decay2', 
                        'out_extent', 'in_extent', 
                        'out_left_base', 'out_right_base', 'in_left_base', 'in_right_base')


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

    @property
    def point_list(self) -> list[Pt]:
        pts = []
        for k in _GEOM_POINT_KEYS:
            v = getattr(self,k)
            if v is not None:
                pts.append(v)
        return pts


@dataclass
class LayoutResult:
    """
    The complete, renderer-ready output of the layout pipeline.

    LayoutResult is the final product of all four layout phases. It collects
    the fully resolved EdgeGeom for every edge in the flux graph, with all
    coordinates in absolute world space and ready to be consumed directly
    by the renderer.

    Attributes
    ----------
    geoms : dict of Edge -> EdgeGeom
        Mapping from each (u, v) edge pair to its resolved geometry. Every
        edge present in the flux graph has an entry; edges with no drawable
        geometry (e.g. internal bulk connectors) may have unpopulated fields.

    Properties
    ----------
    all_points : list of Pt
        Flat list of every non-None coordinate point across all EdgeGeoms,
        including primary rectangle corners, bulk arrow points, decay triangle
        vertices, and helper rectangle vertices. Useful for computing the
        bounding box of the figure or for applying a global coordinate transform.
    """
    geoms: dict[Edge, EdgeGeom]
    _all_points_cache = None

    @property
    def all_points(self) -> list[Pt]:
        if self._all_points_cache is not None:
            return self._all_points_cache
        pts: list[Pt] = []
        for g in self.geoms.values():
            for attr in _GEOM_POINT_KEYS:
                pt = getattr(g, attr)
                if pt is not None:
                    pts.append(pt)
            for rect in g.helper_rects:
                pts.extend(rect)
        self._all_points_cache = pts
        return pts




_IN_EDGE_ORDER: dict[tuple, int] = {
    ('frameshift',         +1, -1): 0,
    ('initiation',    +1):     1,
    ('frameshift',         +1, +1): 2,
    ('load',          +1):     3,
    # 4 reserved for direction == 0
    ('load',          -1):     5,
    ('frameshift',         -1, +1): 6,
    ('40s_retention', -1):     7,
    ('frameshift',         -1, -1): 8,
}

_OUT_EDGE_ORDER: dict[tuple, int] = {
    ('frameshift',         -1, -1): 0,
    ('40s_retention', -1):     1,
    ('frameshift',         -1, +1): 2,
    ('drop',          -1):     3,
    # 4 reserved for direction == 0
    ('drop',          +1):     5,
    ('frameshift',         +1, +1): 6,
    ('initiation',    +1):     7,
    ('frameshift',         +1, -1): 8,
}


"""
_IN_EDGE_ORDER and _OUT_EDGE_ORDER define the vertical stacking priority of
edges on each face of a node.  Lower integers sit closer to the bottom of
the node face.  Slot 4 is reserved for horizontal (direction == 0) edges.

Keys are (etype, direction) or (etype, direction, shift_n) for 'frameshift' edges.
Any combination not present in the table is rejected by _sort_key() with an
assertion error.
"""



def _sort_key(spec: EdgeSpec, direction: int, order: dict) -> tuple[int, int]:
    if direction == 0:
        return (4, 0)
    if spec.etype == 'frameshift':
        priority = order.get(('frameshift', direction, spec.shift_n), 99)
    else:
        priority = order.get((spec.etype, direction), 99)
    assert priority != 99, \
        f'Unknown edge type/direction combo: {spec.etype!r}, {direction!r}'
    return (priority, spec.shift_n)


# ─────────────────────────────────────────────────────────────────────────────
# Coordinate helpers
# ─────────────────────────────────────────────────────────────────────────────

def _shift_pt(pt: Pt, delta: float, axis: Literal['x', 'y']) -> Pt:
    x, y = pt
    return (x + delta, y) if axis == 'x' else (x, y + delta)


def _shift_geom_y(geom: EdgeGeom, delta: float, keys: list[str]) -> None:
    for k in keys:
        pt = getattr(geom, k)
        if pt is not None:
            setattr(geom, k, (pt[0], pt[1] + delta))
    if 'out0' in keys:
        geom.out_helper_rects = [
            [(x, y + delta) for x, y in rect]
            for rect in geom.out_helper_rects
        ]
    if 'in0' in keys:
        geom.in_helper_rects = [
            [(x, y + delta) for x, y in rect]
            for rect in geom.in_helper_rects
        ]

def _shift_geom_x(geom: EdgeGeom, delta: float, keys: list[str]) -> None:
    for k in keys:
        pt = getattr(geom, k)
        if pt is not None:
            setattr(geom, k, (pt[0] + delta, pt[1]))
    if 'out0' in keys:
        geom.out_helper_rects = [
            [(x + delta, y) for x, y in rect]
            for rect in geom.out_helper_rects
        ]
    if 'in0' in keys:
        geom.in_helper_rects = [
            [(x + delta, y) for x, y in rect]
            for rect in geom.in_helper_rects
        ]
class LayoutEngine:
    """
    Converts a RiboGraph flux graph into a renderer-ready LayoutResult.

    The engine runs four sequential phases, each with a single well-defined
    responsibility.  All phase methods are pure with respect to the graph —
    they return new data structures rather than mutating the graph — and can
    be overridden individually in a subclass.

    Pipeline
    --------
    Pre-pass  _node_x_positions
        Map each node's nucleotide position to a world x coordinate,
        optionally compressing large inter-node gaps with a log scale.

    Phase 1   classify_edges  →  dict[Edge, EdgeSpec]
        Label every edge with its semantic type (load, drop, initiation,
        40s_retention, shift, or phase-lane continuation), travel direction,
        and codon displacement.  Bulk edge directions are resolved in a
        second sub-pass once all non-bulk edges are classified.

    Phase 2   order_nodes  →  dict[RiboNode, NodeLayout]
        Sort each node's incoming and outgoing edges into an ordered slot
        list using the priority tables _IN_EDGE_ORDER / _OUT_EDGE_ORDER.
        The slot order determines vertical stacking on the node face.
        No coordinates are assigned here.

    Phase 3   compute_geometries  →  dict[Edge, EdgeGeom]
        Assign in0/in1/out0/out1 corners to every edge band in local-y
        space (x is already world x; y is relative to each node's own
        vertical origin).  Also computes decay triangles, bulk arrowheads,
        and helper rectangles for stacked bulk offsets.

    Phase 4   align_layout  →  LayoutResult
        Three sequential global passes convert local-y to absolute world
        coordinates:

        4a  _align_horizontal
            Walk nodes in topological order, shifting each node upward until
            the out_bot anchor of every outgoing horizontal edge agrees with
            the in_bot anchor of the same edge at its target node.  Repeated
            until convergence (or a maximum-pass safety limit).

        4b  _stack_phases
            Shift each phase lane (0, 1, 2, 3) upward so that the lanes do
            not overlap, separated by a configurable buffer.

        4c  _centre_events
            For non-shift event edges (initiation, 40s_retention), split the
            x gap between source and target nodes equally so the edge is
            centred between its two node faces.

    Parameters
    ----------
    log_scale : float
        Base for the logarithmic x-axis compression applied to gaps between
        node positions.  A value ≤ 1 disables compression and uses raw
        nucleotide positions directly.  Values > 1 reduce the visual impact
        of long inter-node gaps (e.g. log_scale=10 compresses a 1000 nt gap
        to log₁₀(1000)+1 = 4 display units).

    Methods
    -------
    run(graph)
        Execute the full pipeline and return a LayoutResult.
    classify_edges(graph)
        Phase 1 only — useful for unit testing edge classification in isolation.
    order_nodes(graph, specs)
        Phase 2 only.
    compute_geometries(graph, specs, layouts, node_x)
        Phase 3 only.
    align_layout(graph, geoms, layouts, node_x)
        Phase 4 only — returns a LayoutResult.

    Notes
    -----
    - Coordinate conventions: x increases 5'→3' along the mRNA; y increases
      upward; phase lanes are stacked bottom (phase 0) to top (phase 3).
    - The _OWNED class variable maps (is_event, side) pairs to the point
      attributes that belong to that node face, ensuring that _shift_node_geoms
      never moves a coordinate owned by the opposite end of an edge.
    - Phase 4a issues a UserWarning if it fails to converge within
      (4 * node_count + 10) passes; the resulting layout is approximate but
      still renderable.
    """

    def __init__(self, log_scale: float = 1):
        self.log_scale = log_scale



    def _build_edge_index(self, graph):
        out_idx = {n: [] for n in graph.nodes}
        in_idx  = {n: [] for n in graph.nodes}
        # bulk edges anchored at each non-bulk node
        bulk_out_idx = {n: [] for n in graph.nodes}  # node → drop edges (u=node, v=bulk)
        bulk_in_idx  = {n: [] for n in graph.nodes}  # node → load edges (u=bulk, v=node)
        for u, v in graph.edges:
            out_idx[u].append((u, v))
            in_idx[v].append((u, v))
            if v.phase == -1:
                bulk_out_idx[u].append((u, v))
            if u.phase == -1:
                bulk_in_idx[v].append((u, v))
        return out_idx, in_idx, bulk_out_idx, bulk_in_idx

    # ── Entry point ──────────────────────────────────────────────────────────

    def run(self, graph: RiboGraph) -> LayoutResult:
        node_x = self._node_x_positions(graph)
        self._out_idx, self._in_idx, self._bulk_out_idx, self._bulk_in_idx= self._build_edge_index(graph)
        self._topo_order = list(topological_sort(graph))
        specs     = self.classify_edges(graph)
        layouts   = self.order_nodes(graph, specs)
        geoms     = self.compute_geometries(graph, specs, layouts, node_x)
        result    = self.align_layout(graph, geoms, layouts, node_x)
        return result

    # ── Pre-pass: node x positions ───────────────────────────────────────────

    def _node_x_positions(self, graph: RiboGraph) -> dict[RiboNode, float]:
        raw = {n: n.position for n in graph.nodes}
        if self.log_scale <= 1:
            return raw

        xs      = sorted(set(raw.values()))
        gaps    = [xs[i + 1] - xs[i] for i in range(len(xs) - 1)]
        log_g   = [math.log(g, self.log_scale) + 1 for g in gaps]
        log_xs  = [xs[0]]
        for g in log_g:
            log_xs.append(log_xs[-1] + g)
        log_map = dict(zip(xs, log_xs))
        return {n: log_map[raw[n]] for n in graph.nodes}

    # ── Phase 1: classify ────────────────────────────────────────────────────

    def classify_edges(self, graph: RiboGraph) -> dict[Edge, EdgeSpec]:
        """
        Assign etype, direction, and shift_n to every edge.
        Returns a plain dict — the graph is not mutated.
        """
        specs: dict[Edge, EdgeSpec] = {}

        for u, v, data in graph.edges(data=True):
            specs[(u, v)] = self._classify_one(u, v, data)

        # Resolve bulk (load/drop) directions now that all non-bulk are known
        self._resolve_bulk_directions(graph, specs)

        return specs

    def _classify_one(self, u: RiboNode, v: RiboNode,
                      data: dict) -> EdgeSpec:
        is_event = u.phase != v.phase

        if not is_event:
            return EdgeSpec(
                u=u, v=v,
                flux_start=data['flux_start'], flux_end=data['flux_end'],
                etype=str(u.phase), direction=0,
                is_event=False,
            )

        if u.phase == -1:
            etype = 'load'
        elif v.phase == -1:
            etype = 'drop'
        elif u.phase == 0 and v.phase > 0:
            etype = 'initiation'
        elif u.phase > 0 and v.phase == 0:
            etype = '40s_retention'
        else:
            etype = 'frameshift'

        shift_n   = (v.position - u.position) if etype == 'frameshift' else 0
        direction = (
            None if (u.phase == -1 or v.phase == -1)
            else int((v.phase - u.phase) / abs(v.phase - u.phase))
        )

        return EdgeSpec(
            u=u, v=v,
            flux_start=data['flux_start'], flux_end=data['flux_end'],
            etype=etype, direction=direction,
            shift_n=shift_n, is_event=True,
        )

    def _resolve_bulk_directions(self, graph, specs):
        """
        Fill in direction=None on load/drop EdgeSpecs.

        Bulk edge direction is chosen to oppose the net direction of the
        non-bulk edges at the same node: if most traffic goes up  (+1),
        the load arrow points right (-1) so it does not cross
        the main stream.  Called at the end of classify_edges once all
        non-bulk specs are populated.
        """

    # Build per-node lookup once instead of scanning all specs each time
        in_specs_by_node  = {}   # node → list of specs where v == node, u.phase != -1
        out_specs_by_node = {}   # node → list of specs where u == node, v.phase != -1
        for (u, v), s in specs.items():
            if u.phase != -1 and v.phase != -1:
                out_specs_by_node.setdefault(u, []).append(s)
                in_specs_by_node.setdefault(v, []).append(s)

        for node in graph.nodes:
            if node.phase == -1:
                continue
            bulk_edges = self._bulk_edges_at(node, graph)
            for bulk_node, direction in bulk_edges.items():
                key = (bulk_node, node)
                if key in specs and specs[key].direction is None:
                    total = sum(s.direction or 0 for s in in_specs_by_node.get(node, []))
                    specs[key].direction = -1 if total >= 0 else 1

                key = (node, bulk_node)
                if key in specs and specs[key].direction is None:
                    total_out = sum(s.direction or 0 for s in out_specs_by_node.get(node, []))
                    total_in  = sum(s.direction or 0 for s in in_specs_by_node.get(node, []))
                    specs[key].direction = -1 if (total_out - total_in) > 0 else 1

    def _bulk_edges_at(self, node, graph):
        return {
            **{v: 1  for u, v in self._bulk_out_idx[node]},
            **{u: -1 for u, v in self._bulk_in_idx[node]},
        }

    # ── Phase 2: order ───────────────────────────────────────────────────────

    def order_nodes(self, graph: RiboGraph,
                    specs: dict[Edge, EdgeSpec]) -> dict[RiboNode, NodeLayout]:
        """
        Sort edges into in/out slots for each non-bulk node.
        Returns a NodeLayout per node — no coordinates involved.
        """
        layouts: dict[RiboNode, NodeLayout] = {}

        for node in graph.nodes:
            if node.phase == -1:
                continue

            nl = NodeLayout(node=node)

            # in-slots
            in_specs = [
                (specs[(u, v)], specs[(u, v)].direction)
                for u, v in graph.in_edges(node)
                if (u, v) in specs
            ]
            nl.in_slots = [
                EdgeSlot(edge=(s.u, s.v), direction=d)
                for s, d in sorted(
                    in_specs,
                    key=lambda sd: _sort_key(sd[0], sd[1], _IN_EDGE_ORDER))
            ]

            # out-slots
            out_specs = [
                (specs[(u, v)], specs[(u, v)].direction)
                for u, v in graph.out_edges(node)
                if (u, v) in specs
            ]
            nl.out_slots = [
                EdgeSlot(edge=(s.u, s.v), direction=d)
                for s, d in sorted(
                    out_specs,
                    key=lambda sd: _sort_key(sd[0], sd[1], _OUT_EDGE_ORDER))
            ]

            # record drop direction so geometry phase can use it
            for s, _ in out_specs:
                if s.etype == 'drop':
                    nl.drop_direction = s.direction

            layouts[node] = nl

        return layouts

    # ── Phase 3: local geometry ───────────────────────────────────────────────

    def compute_geometries(
        self,
        graph:   RiboGraph,
        specs:   dict[Edge, EdgeSpec],
        layouts: dict[RiboNode, NodeLayout],
        node_x:  dict[RiboNode, float],
    ) -> dict[Edge, EdgeGeom]:
        """
        Assign in0/in1/out0/out1 coordinates to every edge.
        Coordinates are *local-y* at this stage (x is already world x).
        """
        geoms: dict[Edge, EdgeGeom] = {
            (s.u, s.v): EdgeGeom(
                edge=(s.u, s.v),
                etype=s.etype,
                is_event=s.is_event,
                flux_end=s.flux_end,
            )
            for s in specs.values()
        }

        for node, nl in layouts.items():
            x = node_x[node]

            self._fill_side(
                node, nl, specs, geoms, x,
                slots=nl.in_slots, side='left', flux_key='flux_end',
            )
            self._fill_side(
                node, nl, specs, geoms, x,
                slots=nl.out_slots, side='right', flux_key='flux_start',
            )

        self._resolve_decay_anchors(layouts, geoms)

        self._fill_bulk_geoms(graph, specs, geoms, node_x, bulk_length_factor=1.5)

        return geoms

    def _fill_side(

        self,
        node:     RiboNode,
        nl:       NodeLayout,
        specs:    dict[Edge, EdgeSpec],
        geoms:    dict[Edge, EdgeGeom],
        x:        float,
        slots:    list[EdgeSlot],
        side:     Literal['left', 'right'],
        flux_key: Literal['flux_start', 'flux_end'],
    ) -> None:

        """
        Assign coordinates to all edges on one face (left=in / right=out)
        of *node* in local-y space.

        Edges are processed in three passes — bottom diagonals, horizontal,
        top diagonals — matching the slot ordering from Phase 2.  Helper
        rectangles are appended for any diagonal edge that has an x offset
        (i.e. sits behind another diagonal at the same face).
        """

        sign     = -1 if side == 'left' else 1
        inout    = 'in' if side == 'left' else 'out'
        bottom_d = -1   if side == 'right' else 1
        top_d    = -bottom_d

        x_offset  = 0.0
        current_y = 0.0
        decay     = 0.0

        def flux(slot: EdgeSlot) -> float:
            s = specs[slot.edge]
            return s.flux_end if flux_key == 'flux_end' else s.flux_start

        def geom(slot: EdgeSlot) -> EdgeGeom:
            return geoms[slot.edge]

        def spec(slot: EdgeSlot) -> EdgeSpec:
            return specs[slot.edge]

        # ── bottom diagonals ─────────────────────────────────────────────────
        for slot in [s for s in slots if s.direction == bottom_d]:
            f = flux(slot)
            g = geom(slot)
            setattr(g, inout + '0', (x + x_offset,              current_y))
            setattr(g, inout + '1', (x + x_offset + sign * f,   current_y))
            setattr(g, inout + '_quad', 2 if side == 'left' else 1)
            if x_offset:
                self._add_helper_rect(g, x, current_y, x_offset, f, side)
            x_offset  += sign * f
            current_y += f

        # ── horizontal ───────────────────────────────────────────────────────
        h_slots = [s for s in slots if s.direction == 0]
        assert len(h_slots) <= 1

        for slot in h_slots:
            g  = geom(slot)
            s  = spec(slot)
            f  = s.flux_end if flux_key == 'flux_end' else s.flux_start
            decay = s.flux_start - s.flux_end
            setattr(g, inout + '_bot', (x, current_y))

            if side == 'left':
                if nl.drop_direction is not None:
                    if nl.drop_direction == 1:
                        g.in1, g.in0     = (x, current_y),         (x, current_y + f)
                        g.decay1, g.decay2 = (x, current_y + f),   (x, current_y + f + decay)
                    else:
                        g.in1, g.in0     = (x, current_y + decay), (x, current_y + decay + f)
                        g.decay1, g.decay2 = (x, current_y),       (x, current_y + decay)
                    current_y += f + decay
                    continue

                setattr(g, inout + '1', (x, current_y))
                setattr(g, inout + '0', (x, current_y + f))

            else:
                g.out0 = (x, current_y)
                g.out1 = (x, current_y + f)
                # decay anchor is set by the downstream node's drop_direction
                v_node = slot.edge[1]
                if v_node in {sl.edge[0] for sl in slots}:
                    pass   # resolved in a second pass if needed

            current_y += f + decay

        # ── top diagonals ────────────────────────────────────────────────────
        x_offset  = 0.0
        add        = decay if side == 'left' else 0.0
        current_y  = sum(flux(s) for s in slots) + add

        for slot in reversed([s for s in slots if s.direction == top_d]):
            f = flux(slot)
            g = geom(slot)
            setattr(g, inout + '0', (x + x_offset,             current_y))
            setattr(g, inout + '1', (x + x_offset + sign * f,  current_y))
            setattr(g, inout + '_quad', 3 if side == 'left' else 4)
            if x_offset:
                self._add_helper_rect(g, x, current_y, x_offset, -f, side)
            x_offset  += sign * f
            current_y -= f

    def _fill_bulk_geoms(
        self,
        graph:              RiboGraph,
        specs:              dict[Edge, EdgeSpec],
        geoms:              dict[Edge, EdgeGeom],
        node_x:             dict[RiboNode, float],
        bulk_length_factor: float,
    ) -> None:
        
        """
        Add arrowhead geometry (extent and base points) to load and drop edges.

        Arrow length scales with flux_start * bulk_length_factor, with a
        minimum of 0.2 display units to keep zero-flux arrows visible.
        Arrow depth scales with flux via ARROW_DEPTH_SCALE.
        """

        DELTA_SIGN = {
            (-1, 'load'): +1, (-1, 'drop'): -1,
            (+1, 'load'): -1, (+1, 'drop'): +1,
        }
        for (u, v), s in specs.items():
            if u.phase != -1 and v.phase != -1:
                continue
            g      = geoms[(u, v)]
            length = max(s.flux_start * bulk_length_factor, 0.2)
            delta  = length * DELTA_SIGN[(s.direction, s.etype)]

            change_side = 'in' if v.phase == -1 else 'out'
            get_side    = 'out' if change_side == 'in' else 'in'

            for pos, opp in (('0', '1'), ('1', '0')):
                pt = getattr(g, f'{get_side}{pos}')
                if pt is not None:
                    setattr(g, f'{change_side}{opp}', _shift_pt(pt, delta, 'y'))
    
            flux = s.flux_end
            direction = DELTA_SIGN[(s.direction, s.etype)]
            depth = flux * ARROW_DEPTH_SCALE
            if s.etype == 'load' and g.out0 and g.out1:
                # anchor at midpoint of OUT edge
                x0, y0 = g.out0
                x1, y1 = g.out1

                g.out_extent = ((x0 + x1) / 2, y0 + ARROW_DEPTH_SCALE * flux * direction*2)
                g.out_left_base = (g.out0[0] - depth, g.out_extent[1])
                g.out_right_base =(g.out1[0] + depth, g.out_extent[1])

            elif s.etype == 'drop' and g.in0 and g.in1:
                # anchor at midpoint of IN edge
                x0, y0 = g.in0
                x1, y1 = g.in1

                g.in_extent = ((x0 + x1) / 2, y0 + ARROW_DEPTH_SCALE * flux * direction*2)
                g.in_left_base = (g.in0[0] + depth, y0)
                g.in_right_base =(g.in1[0] - depth, y0)



    def _resolve_decay_anchors(
        self,
        layouts: dict[RiboNode, NodeLayout],
        geoms:   dict[Edge, EdgeGeom],
    ) -> None:
        """
        Set decay0 on every horizontal (non-event) edge.

        decay0 is the y-anchor where the decay wedge meets the horizontal
        stream at the source node.  Its position depends on the downstream
        node's drop_direction, which is only known after Phase 2 has run
        for every node, so this is a deferred second pass within Phase 3.

        Rule:
          drop_direction == -1  →  decay0 = out0  (bottom of stream)
          drop_direction == +1  →  decay0 = out1  (top of stream)
          no drop at v          →  decay0 left unset
        """
        for (u, v), g in geoms.items():
            if g.is_event:
                continue                          # only horizontal edges
            if g.out0 is None or g.out1 is None:
                continue                          # geometry not filled yet
 
            v_layout = layouts.get(v)
            if v_layout is None or v_layout.drop_direction is None:
                continue                          # downstream node has no drop
 
            g.decay0 = g.out0 if v_layout.drop_direction == -1 else g.out1
 

    @staticmethod
    def _add_helper_rect(g: EdgeGeom, x: float, y: float,
                        width: float, height: float, side:str ) -> None:
        if side == 'left':
            g.in_helper_rects.append([
                (x,         y),
                (x + width, y),
                (x + width, y + height),
                (x,         y + height),
            ])
        elif side == 'right':
            g.out_helper_rects.append([
                (x,         y),
                (x + width, y),
                (x + width, y + height),
                (x,         y + height),
            ])


    # ── Phase 4: global alignment ─────────────────────────────────────────────

    def align_layout(
        self,
        graph:   RiboGraph,
        geoms:   dict[Edge, EdgeGeom],
        layouts: dict[RiboNode, NodeLayout],
        node_x:  dict[RiboNode, float],
    ) -> LayoutResult:
        """
        Three sequential global alignment passes:
          4a. Align horizontal edges (make out_bot == in_bot)
          4b. Stack phases vertically with a buffer
          4c. Centre non-shift event edges horizontally
        Each pass calls _shift_node_geoms() — no direct dict mutation.
        """
        self._align_horizontal(graph, geoms, layouts)
        self._diagnose_horizontal(geoms)   # ← add this
        self._stack_phases(graph, geoms, layouts, buffer=1.5)
        self._centre_events(graph, geoms, layouts)
        return LayoutResult(geoms=geoms)
    

    def _diagnose_horizontal(self, geoms):
        for (u, v), g in geoms.items():
            if g.is_event or g.out_bot is None or g.in_bot is None:
                continue
            diff = abs(g.out_bot[1] - g.in_bot[1])
            if diff > 1e-6:
                print(f"MISALIGNED {u} → {v}: out_bot={g.out_bot[1]:.4f} in_bot={g.in_bot[1]:.4f} diff={diff:.4f}")
    
    def _align_horizontal(self, graph, geoms, layouts):
        MAX_PASSES = len(list(graph.nodes)) * 4 + 10
        watch = {  # positions from your misaligned edges
            (150, 3), (248, 2)
        }
        
        for pass_num in range(MAX_PASSES):
            changed = False
            for node in self._topo_order:
                if node.phase == -1:
                    continue
                
                max_out_delta = 0.0
                max_in_delta  = 0.0

                for slot in layouts[node].out_slots:
                    g = geoms[slot.edge]
                    if g.is_event or g.out_bot is None or g.in_bot is None:
                        continue
                    agreed = max(g.out_bot[1], g.in_bot[1])
                    if agreed > g.out_bot[1]:
                        max_out_delta = max(max_out_delta, agreed - g.out_bot[1])

                for slot in layouts[node].in_slots:
                    g = geoms[slot.edge]
                    if g.is_event or g.out_bot is None or g.in_bot is None:
                        continue
                    agreed = max(g.out_bot[1], g.in_bot[1])
                    if agreed > g.in_bot[1]:
                        max_in_delta = max(max_in_delta, agreed - g.in_bot[1])

                delta = max(max_out_delta, max_in_delta)
                
                
                if delta > 1e-9:
                    self._shift_node_geoms(node, delta, 'y', geoms, layouts, graph)
                    changed = True

            if not changed:
                break
                
    # ── 4b: phase stacking ────────────────────────────────────────────────────

    def _stack_phases(
        self,
        graph:   RiboGraph,
        geoms:   dict[Edge, EdgeGeom],
        layouts: dict[RiboNode, NodeLayout],
        buffer:  float,
    ) -> None:
        def phase_nodes(phase: int) -> list[RiboNode]:
            core = [n for n in graph.nodes if n.phase == phase]

            adj_m1 = [
                n for n in graph.nodes
                if n.phase == -1 and any(
                    graph.has_edge(n, nb) and nb.phase == phase
                    for nb in graph.successors(n)
                )
            ]

            adj_m2 = [
                n for n in graph.nodes
                if n.phase == -1 and any(
                    u.phase == phase for u, _ in graph.in_edges(n)
                )
            ]

            return core + adj_m1 + adj_m2

        def min_y(nodes: list[RiboNode]) -> float:
            pts = self._points_exclusive(nodes, geoms, graph)
            return min(p[1] for p in pts) if pts else 0.0

        def max_y(nodes: list[RiboNode]) -> float:
            pts = self._points_exclusive(nodes, geoms, graph)
            return max(p[1] for p in pts) if pts else 0.0

        core0 = [n for n in graph.nodes if n.phase == 0]
        p0    = phase_nodes(0)
        if p0:
            shift = -min_y(p0)
            for n in core0:
                self._shift_node_geoms(n, shift, 'y', geoms, layouts, graph)

        skip       = 0 if core0 else 1
        prev_nodes = p0

        for phase in range(1, 4):
            core_p = [n for n in graph.nodes if n.phase == phase]
            if not core_p:
                skip += 1
                continue
            p_nodes = phase_nodes(phase)
            shift   = (max_y(prev_nodes) + buffer + skip) - min_y(p_nodes)
            for n in core_p:
                self._shift_node_geoms(n, shift, 'y', geoms, layouts, graph)
            prev_nodes = p_nodes
            skip       = 0

    # ── 4c: event x-centring ─────────────────────────────────────────────────

    def _centre_events(
        self,
        graph:   RiboGraph,
        geoms:   dict[Edge, EdgeGeom],
        layouts: dict[RiboNode, NodeLayout],
    ) -> None:
        for (u, v), g in geoms.items():
            if u.phase == -1 or v.phase == -1:
                continue
            if not g.is_event or g.etype == 'frameshift':
                continue
            # Determine x-gap based on edge direction
            # direction +1 → in1/out0, direction -1 → in0/out1
            if g.in1 and g.out0:
                diff = g.in1[0] - g.out0[0]
            elif g.in0 and g.out1:
                diff = g.in0[0] - g.out1[0]
            else:
                continue
            if diff == 0:
                continue
            self._shift_node_geoms(u,  diff / 2, 'x', geoms, layouts, graph)
            self._shift_node_geoms(v, -diff / 2, 'x', geoms, layouts, graph)

    # ── Shared shift helper ───────────────────────────────────────────────────

    # Keys owned by each (is_event, side) combination
    _OWNED: dict[tuple[bool, str], list[str]] = {
        (False, 'out'): ['out0', 'out1', 'decay0', 'out_bot'],
        (False, 'in'):  ['in0',  'in1',  'decay1', 'decay2', 'in_bot'],
        (True,  'out'): ['out0', 'out1', 'decay0', 'out_extent',  'out_left_base', 'out_right_base'],
        (True,  'in'):  ['in0',  'in1',  'decay1', 'decay2', 'in_extent', 'in_left_base', 'in_right_base'],
    }

    def _shift_node_geoms(
        self,
        node:    RiboNode,
        delta:   float,
        axis:    Literal['x', 'y'],
        geoms:   dict[Edge, EdgeGeom],
        layouts: dict[RiboNode, NodeLayout],
        graph:   RiboGraph,
    ) -> None:
        """
        Shift all geometry coordinates *owned by node* by delta along axis.

        Ownership is defined by _OWNED: each (is_event, side) pair maps to
        the specific point attributes that belong to that node face.  This
        ensures that shifting node A never moves the far-end coordinates of
        an edge that are owned by node B.  Bulk edge endpoints anchored at
        this node are handled separately because bulk nodes have no NodeLayout.
        """
        for edges, side in (
                (self._out_idx[node], 'out'),
                (self._in_idx[node],  'in'),
            ):
                for u, v in edges:
                    g = geoms.get((u, v))
                    if g is None:
                        continue
                    keys = self._OWNED[(g.is_event, side)]
                    (_shift_geom_y if axis == 'y' else _shift_geom_x)(g, delta, keys)

        # bulk — now O(degree) not O(E)
        for u, v in self._bulk_out_idx[node]:
            g = geoms.get((u, v))
            if g:
                (_shift_geom_y if axis == 'y' else _shift_geom_x)(g, delta, ['in1', 'in0', 'in_extent', 'in_left_base', 'in_right_base'])
        for u, v in self._bulk_in_idx[node]:
            g = geoms.get((u, v))
            if g:
                (_shift_geom_y if axis == 'y' else _shift_geom_x)(g, delta, ['out1', 'out0', 'out_extent', 'out_left_base', 'out_right_base'])

    # ── Point-set helpers ─────────────────────────────────────────────────────



    def _points_exclusive(
        self,
        nodes: list[RiboNode],
        geoms: dict[Edge, EdgeGeom],
        graph: RiboGraph,
    ) -> list[Pt]:
        node_set = set(nodes)
        return [
            pt
            for (u, v), g in geoms.items()
            if u in node_set and v in node_set
            for key in _GEOM_POINT_KEYS
            if (pt := getattr(g, key)) is not None
        ]

"""
ribo_graph_vis.py  –  Visualisation layer for RiboGraphFlux.

Architecture
------------
RiboRenderer          – stateless factory: LayoutResult → Figure
  └─ EdgePainter      – per-edge patch builder (body / taper / decay / arrow)

To change how any edge type looks, edit EdgePainter or override
COLOR_DICT / STYLE_OVERRIDES / edge_style() on RiboRenderer.
The layout mathematics (LayoutEngine) is entirely separate.
"""

from dataclasses import dataclass
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.path import Path

# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class EdgeStyle:
    """Visual properties for a single rendered edge segment."""
    facecolor: str
    edgecolor: str
    alpha: float = 1.0
    linewidth: float = 0.5
    zorder: int = 2


@dataclass
class RenderPrimitive:
    """A single matplotlib patch to be added to the axes."""
    patch: mpatches.Patch
    zorder: int = 2


# ─────────────────────────────────────────────────────────────────────────────
# EdgePainter  –  builds patches for one edge
# ─────────────────────────────────────────────────────────────────────────────

class EdgePainter:
    """
    Converts the pre-computed geometry for a single edge into renderable patches.

    Instantiate with an EdgeGeom and an EdgeStyle, then call primitives() to
    get the complete list of RenderPrimitives for that edge.  Each primitive
    maps to one PathPatch on the axes.

    This is the only class that needs to change to alter edge appearance.

    Patch types produced
    --------------------
    Body rect       The main flux band quadrilateral (in0/in1 → out0/out1).
    Helper rects    Auxiliary fill rectangles for stacked bulk-edge offsets.
    Decay wedge     Triangle marking flux lost to ribosome drop-off along a
                    scanning or translation edge (purple, zorder - 1).
    Curved tapers   Bézier caps at the node faces of event edges (initiation,
                    40s_retention, shift), smoothing the band entry/exit.
    Bulk arrow      Isoceles triangle arrowhead at the free end of load/drop
                    edges pointing toward (load) or away from (drop) the mRNA.

    Class attributes
    ----------------
    TAPER_DEPTH_SCALE : float
        Controls the x-extent of the Bézier taper cap relative to flux width.
    """

    TAPER_DEPTH_SCALE: float = 0.4

    def __init__(self, geom: EdgeGeom, style: EdgeStyle):
        self.geom  = geom
        self.u, self.v = geom.edge
        self.style = style

    # ── Public entry point ───────────────────────────────────────────────────

    def primitives(self) -> list[RenderPrimitive]:
        """Return all patches needed to draw this edge, in paint order."""
        
        out: list[RenderPrimitive] = []

        out.extend(self._body_primitives())
        out.extend(self._helper_rect_primitives())
        out.extend(self._decay_primitives())
        out.extend(self._taper_primitives())
        out.extend(self._bulk_arrow_primitives())

        return out

    # ── Body rect ────────────────────────────────────────────────────────────

    def _body_primitives(self) -> list[RenderPrimitive]:
        """Main flux band as a filled quadrilateral (in0/in1/out0/out1)."""
        g = self.geom
        corners = [g.in0, g.in1, g.out0, g.out1]
        
        patch = self._rect_patch(corners)
        return [RenderPrimitive(patch, zorder=self.style.zorder)]

    # ── Helper rects (stacked bulk offsets) ──────────────────────────────────

    def _helper_rect_primitives(self) -> list[RenderPrimitive]:
        """Fill rectangles covering x-offset gaps behind stacked bulk edges."""
        return [
            RenderPrimitive(self._rect_patch(r), zorder=self.style.zorder)
            for r in self.geom.helper_rects
        ]
    # ── Decay wedge ──────────────────────────────────────────────────────────

    def _decay_primitives(self) -> list[RenderPrimitive]:
        """
        Triangular wedge (decay0/decay1/decay2) representing flux lost to
        ribosome drop-off.  Drawn in purple at zorder - 1 so it sits behind
        the main band.  Returns nothing when all three vertices are unset or
        when decay1 == decay2 (zero decay).
        """
        g = self.geom

        if g.decay0 is None or g.decay1 is None or g.decay2 is None:

            return []

        if g.decay1 == g.decay2:
            return []
        patch = self._rect_patch(
            [g.decay0, g.decay1, g.decay2],
            facecolor='purple',
            edgecolor='purple',
        )
        return [RenderPrimitive(patch, zorder=self.style.zorder - 1)]

    # ── Curved tapers (event edges only) ─────────────────────────────────────

    def _taper_primitives(self) -> list[RenderPrimitive]:
        """
        Bézier taper caps for event edges (is_event == True).

        One cap is drawn at the source node face (out0, out_quad) and one at
        the target node face (in0, in_quad), provided the respective node is
        not a bulk node.  Caps are rendered at zorder + 1 to sit on top of
        the body rect.  Returns nothing for non-event (horizontal) edges.
        """
        g = self.geom

        if not g.is_event:
            return []

        out: list[RenderPrimitive] = []
        flux = g.flux_end

        if self.v.phase != -1 and g.in0 is not None and g.in_quad is not None:
            patch = self._taper_patch(
                centre=g.in0,
                flux=flux,
                quad=g.in_quad,
                inout='in',
            )
            out.append(RenderPrimitive(patch, zorder=self.style.zorder + 1))

        if self.u.phase != -1 and g.out0 is not None and g.out_quad is not None:
            patch = self._taper_patch(
                centre=g.out0,
                flux=flux,
                quad=g.out_quad,
                inout='out',
            )
            out.append(RenderPrimitive(patch, zorder=self.style.zorder + 1))

        return out
    

    # ── Bulk arrow ───────────────────────────────────────────────────────────
 
    def _bulk_arrow_primitives(self) -> list[RenderPrimitive]:
        """
        Arrowhead at the free (bulk) end of a load or drop edge.

        Load edge (u.phase == -1): triangle points toward the mRNA, with its
        tip between out0 and out1 and its base extending outward.

        Drop edge (v.phase == -1): triangle points away from the mRNA, with
        its tip at in_extent and its base at in_left_base/in_right_base.

        Returns nothing for edges that connect two non-bulk nodes.
        """
        g = self.geom

        is_load = self.u.phase == -1
        is_drop = self.v.phase == -1



        if not (is_load or is_drop):
            return []

        if is_load:     
            verts = [g.out0, g.out_left_base, g.out_right_base, g.out1, g.out0]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        else:
            verts = [g.in_left_base, g.in_right_base, g.in_extent, g.in_left_base]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        
        patch  = mpatches.PathPatch(
            Path(verts, codes),
            **self._patch_kwargs(facecolor='purple', edgecolor='purple'),
        )
        return [RenderPrimitive(patch, zorder=self.style.zorder + 1)]
 


    # ── Patch constructors ───────────────────────────────────────────────────

    def _patch_kwargs(self, **overrides) -> dict:
        kw = dict(
            facecolor=self.style.facecolor,
            edgecolor=self.style.edgecolor,
            alpha=self.style.alpha,
            linewidth=0,
            snap=False
        )
        kw.update(overrides)
        return kw

    def _rect_patch(self, points: list[Pt|None], **overrides) -> mpatches.PathPatch:
        pts = [p for p in points if p is not None]
        if len(pts) < 3:
            raise ValueError('rect_patch needs ≥ 3 valid points')

        verts = pts + [pts[0]]
        
        codes = [Path.MOVETO] + [Path.LINETO] * (len(pts) - 1) + [Path.CLOSEPOLY]

        return mpatches.PathPatch(Path(verts, codes), **self._patch_kwargs(**overrides))

    def _taper_patch(
        self,
        centre: tuple[float, float],
        flux: float,
        quad: int,
        inout: str,
    ) -> mpatches.PathPatch:
        """
        Single Bézier taper cap anchored at *centre*.

        The cap is a CURVE3 quadratic Bézier that rounds the corner where a
        diagonal event edge meets a node face.  *quad* (1–4) encodes which
        quadrant the band approaches from, controlling the x/y sign of the
        control point offset.  *inout* selects whether the cap faces the
        incoming or outgoing side.
        """
        x, y   = centre
        x_sign = -1 if quad in (1, 4) else 1
        y_sign =  1 if quad in (1, 2) else -1

        if inout == 'out':
            verts = [
                (x, y),
                (x,                  y + y_sign * flux),
                (x - x_sign * flux,  y + y_sign * flux),
                (x - x_sign * flux,  y),
            ]
        else:  # 'in'
            verts = [
                (x,                  y),
                (x,                  y + y_sign * flux),
                (x - x_sign * flux,  y + y_sign * flux),
                (x - x_sign * flux,  y),
            ]

        codes = [
            Path.MOVETO, Path.LINETO,
            Path.CURVE3, Path.CURVE3,
        ]
        return mpatches.PathPatch(Path(verts, codes), **self._patch_kwargs())


# ─────────────────────────────────────────────────────────────────────────────
# RiboRenderer  –  assembles all patches for a fully-laid-out graph
# ─────────────────────────────────────────────────────────────────────────────

class RiboRenderer:
    """
    Stateless renderer: converts a LayoutResult into a matplotlib Figure.

    Usage
    -----
    >>> renderer = RiboRenderer()
    >>> fig = renderer.render(layout, node_x=node_x_dict)

    Customisation
    -------------
    Override COLOR_DICT to change per-type colours.  Override STYLE_OVERRIDES
    to merge additional EdgeStyle kwargs (alpha, linewidth, zorder) on top of
    the colour for specific edge types.  For full control, override edge_style()
    to replace style resolution entirely.  To change patch shapes or geometry,
    subclass EdgePainter and override _collect_primitives() to use it.

    Class attributes
    ----------------
    COLOR_DICT : dict[str, str]
        Maps EdgeType strings to matplotlib colour strings.
    STYLE_OVERRIDES : dict[str, dict]
        Optional per-type keyword overrides merged into EdgeStyle on top of
        the colour from COLOR_DICT.  Empty by default.

    Parameters
    ----------
    fig_size : tuple of float
        Matplotlib figure size in inches (width, height). Default (12, 6).
    dpi : int
        Figure resolution. Default 150.
    """

    # ── Colour palette ───────────────────────────────────────────────────────

    COLOR_DICT: dict[str, str] = {
        'frameshift':           'red',
        '40s_retention':   'orange',
        'drop':            'purple',
        'initiation':      'green',
        'load':            'purple',
        'decay':           'purple',
        '0':               'gray',
        '1':               'lightblue',
        '2':               'steelblue',
        '3':               'darkblue',
    }

    # ── Per-type style overrides (merged on top of the default) ──────────────

    STYLE_OVERRIDES: dict[str, dict] = {
        # e.g. 'frameshift': {'alpha': 0.8, 'linewidth': 1.0},
    }

    # ── Construction / entry point ───────────────────────────────────────────

    def __init__(self, fig_size: tuple = (12, 6), dpi: int = 150):
        self.fig_size = fig_size
        self.dpi      = dpi

    def _draw_position_labels(
        self,
        ax,
        layout: LayoutResult,
        node_x: dict,          # RiboNode → world-x float
    ) -> None:
        positions = layout.all_points
        if not positions:
            return
        y_min = min(p[1] for p in positions)
        label_y = y_min - 1.5          # sits just below the lowest geometry

        seen_x: set[float] = set()
        for node, wx in node_x.items():
            if node.phase == -1:       # skip bulk/decay pseudo-nodes
                continue
            if wx in seen_x:           # multiple phases share an x — label once
                continue
            seen_x.add(wx)
            ax.text(
                wx, label_y,
                str(node.position),    # nucleotide coordinate
                ha='center', va='top',
                fontsize=6, color='#444444',
                clip_on=False,
            )

    def render(self, layout: LayoutResult, node_x: dict | None = None) -> Figure:
        fig = plt.figure(figsize=self.fig_size, dpi=self.dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_aspect('equal')
        ax.axis('off')
        self._set_axis_limits(ax, layout)
        primitives = self._collect_primitives(layout)
        for prim in sorted(primitives, key=lambda p: p.zorder):
            prim.patch.set_transform(ax.transData)
            ax.add_patch(prim.patch)
        if node_x is not None:
            self._draw_position_labels(ax, layout, node_x)
        return fig

    # ── Primitive collection ─────────────────────────────────────────────────

    def _collect_primitives(self, layout: LayoutResult) -> list[RenderPrimitive]:
        out: list[RenderPrimitive] = []

        for geom in layout.geoms.values():
            style   = self.edge_style(geom)
            painter = EdgePainter(geom, style)
            out.extend(painter.primitives())

        return out
    
    def _repr_mimebundle_(self, **kwargs):
        return {}
    
    # ── Style resolution ─────────────────────────────────────────────────────

    def edge_style(self, geom: EdgeGeom) -> EdgeStyle:
        etype     = geom.etype
        color     = self.COLOR_DICT.get(etype, 'black')
        overrides = self.STYLE_OVERRIDES.get(etype, {})

        return EdgeStyle(
            facecolor=color,
            edgecolor=color,
            **overrides,
        )

    # ── Helpers ──────────────────────────────────────────────────────────────

    @staticmethod
    def _set_axis_limits(ax, layout: 'LayoutResult') -> None:
        positions = layout.all_points
        if not positions:
            return
        xs = [p[0] for p in positions]
        ys = [p[1] for p in positions]
        ax.set(xlim=(min(xs), max(xs)), ylim=(min(ys), max(ys)))

"""
ribo_graph_vis_layout.py  –  Orchestration layer for RiboGraphFlux visualisation.

Wires together LayoutEngine (geometry) and RiboRenderer (patches → Figure).
See ribo_graph_vis_layout_engine.py for the four-phase layout pipeline and
ribo_graph_vis.py for the rendering layer.
"""





class RiboGraphVis(RiboGraph):
    """
    A RiboGraph subclass that owns the layout and rendered figure for a
    RiboGraphFlux instance.

    On construction, RiboGraphVis copies the flux graph's nodes and edges,
    strips bulk-to-bulk recycling edges, runs the full layout pipeline, and
    renders the result into a matplotlib Figure — all in one step.  The graph
    itself holds no mutable layout state; all geometry lives in layout_result.

    Parameters
    ----------
    incoming_graph_data : RiboGraphFlux
        The flux graph to visualise.  Nodes, edges, and edge data are copied;
        the original graph is not mutated.
    fig_size : tuple of float
        Matplotlib figure size in inches (width, height).  Passed to
        RiboRenderer if no custom renderer is supplied.  Default (12, 6).
    dpi : int
        Figure resolution.  Passed to RiboRenderer if no custom renderer is
        supplied.  Default 150.
    log_scale : float
        Log-compression base for the x axis.  Passed to LayoutEngine if no
        custom engine is supplied.  Values ≤ 1 use raw nucleotide positions.
    engine : LayoutEngine or None
        Custom layout engine.  If None, a default LayoutEngine is constructed
        from log_scale.  Supply a subclass to override any layout phase.
    renderer : RiboRenderer or None
        Custom renderer.  If None, a default RiboRenderer is constructed from
        fig_size and dpi.  Supply a subclass to override colours or patch shapes.

    Attributes
    ----------
    layout_result : LayoutResult
        World-space geometry for every edge, produced by the layout engine.
        Recomputed on each call to compute_layout().
    fig : matplotlib.figure.Figure
        The rendered figure.  Available immediately after construction and
        after any subsequent compute_layout() call.

    Methods
    -------
    compute_layout()
        Re-run the layout pipeline and re-render the figure.  Call this after
        manually modifying graph edges if you need a fresh figure.
    get_figure()
        Return the current Figure object.
    show()
        Display the figure.  Uses IPython display() inside a notebook;
        falls back to plt.show(block=True) in scripts.
    save(filename, dpi, format, **kwargs)
        Save the figure to disk via Figure.savefig().
    positions
        Flat list of all world-space (x, y) points in the layout, useful for
        computing bounding boxes or applying coordinate transforms.

    Notes
    -----
    - Bulk-to-bulk recycling edges (phase == -1 on both ends) are removed
      before layout because they carry no positional information and would
      produce degenerate geometry.
    - Isolated nodes left after pruning are also removed.
    - LayoutEngine._node_x_positions() is called a second time after
      engine.run() to pass the same x mapping to the renderer for axis labels.
      This is a known redundancy; avoid overriding that method with side effects.
    """

    def __init__(
        self,
        graph_to_render: FluxGraph|SimpleFluxGraph,
        fig_size:  tuple = (12, 6),
        dpi:       int   = 150,
        log_scale: float = 1,
        engine:    LayoutEngine  | None = None,
        renderer:   RiboRenderer | None = None,   # RiboRenderer from previous module
        **attr,
    ):
        self.engine   = engine   or LayoutEngine(log_scale=log_scale)
        self.renderer = renderer or RiboRenderer(fig_size=fig_size, dpi=dpi)

        super().__init__(**attr)
        if isinstance(graph_to_render, SimpleFluxGraph):
            source = graph_to_render
        elif isinstance(graph_to_render, FluxGraph):
            source = graph_to_render.simple
        source.prune_recycle_edges()
        for node in source.nodes:
            self.add_node(node)
        for u, v, data in source.edges(data=True):
            self.add_edge(u, v, **data)
        self.compute_layout()

    # ── Public API ───────────────────────────────────────────────────────────


    def compute_layout(self) -> None:
        """Run the layout pipeline and re-render the figure.

        Stores the result in self.layout_result and self.fig.  Safe to call
        more than once, e.g. after programmatic edge changes.
        """
        self.layout_result = self.engine.run(self)
        node_x = self.engine._node_x_positions(self)   # same call the engine made
        self.fig: Figure = self.renderer.render(self.layout_result, node_x=node_x)

    def get_figure(self) -> Figure:
        """Return the matplotlib Figure object"""
        return self.fig

    def show(self) -> None:
        if self.fig is None:
            raise RuntimeError("Call compute_layout() first")
        try:
            from IPython.core.getipython import get_ipython
            if get_ipython() is not None:
                from IPython.display import display
                display(self.fig)
                return
        except Exception:
            pass
        import matplotlib.pyplot as plt
        plt.figure(self.fig)
        plt.show(block=True)

    def save(self, filename='output.png', dpi=150, format=None, **kwargs):
        """
        Save the figure to a file.

        Parameters
        ----------
        filename : str
            Output path including extension.  Default 'output.png'.
        dpi : int
            Override the figure's own dpi for the saved file.  Default 150.
        format : str or None
            File format string (e.g. 'pdf', 'svg').  Inferred from the
            filename extension when None.
        **kwargs
            Passed through to Figure.savefig().
        """
        if self.fig is None:
            raise RuntimeError("Call compute_layout() first")
        self.fig.savefig(filename, dpi=dpi, format=format, **kwargs)


    @property
    def positions(self) -> list[Pt]:
        """All world-space (x, y) points in the current layout."""
        return self.layout_result.all_points
