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

from __future__ import annotations
from dataclasses import dataclass
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.path import Path
from ..data import EdgeGeom, LayoutResult

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

    def _rect_patch(self, points: list[tuple], **overrides) -> mpatches.PathPatch:
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

    def render(self, layout: LayoutResult,
            node_x: dict | None = None) -> Figure:
        fig, ax = plt.subplots(figsize=self.fig_size, dpi=self.dpi)
        ax.set_aspect('equal')
        ax.axis('off')

        self._set_axis_limits(ax, layout)

        primitives = self._collect_primitives(layout)
        for prim in sorted(primitives, key=lambda p: p.zorder):
            prim.patch.set_transform(ax.transData)
            ax.add_patch(prim.patch)

        if node_x is not None:
            self._draw_position_labels(ax, layout, node_x)

        fig.tight_layout()
        return fig

    # ── Primitive collection ─────────────────────────────────────────────────

    def _collect_primitives(self, layout: LayoutResult) -> list[RenderPrimitive]:
        out: list[RenderPrimitive] = []

        for geom in layout.geoms.values():
            style   = self.edge_style(geom)
            painter = EdgePainter(geom, style)
            out.extend(painter.primitives())

        return out
    
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

