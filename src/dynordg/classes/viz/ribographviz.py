"""
ribo_graph_vis_layout.py  –  Layout engine for RiboGraphFlux.

Pipeline
--------
RiboGraphVis.compute_layout() runs four explicit phases, each of which
consumes typed dataclasses and produces typed dataclasses:

  Phase 1  classify_edges()      RiboGraphFlux  →  dict[edge, EdgeSpec]
  Phase 2  order_nodes()         EdgeSpec dict  →  dict[node, NodeLayout]
  Phase 3  compute_geometries()  NodeLayout     →  dict[edge, EdgeGeom]
  Phase 4  align_layout()        EdgeGeom dict  →  LayoutResult

"""

from __future__ import annotations
from ..graph import RiboGraph
from ..simulation import RiboGraphFlux
from .data import LayoutResult, Pt
from .render import LayoutEngine, RiboRenderer
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


class RiboGraphVis(RiboGraph):
    """
    Orchestrates layout (via LayoutEngine) and rendering (via RiboRenderer).
    The graph itself stores no layout state — all geometry lives in
    self.layout_result.
    """

    def __init__(
        self,
        incoming_graph_data: RiboGraphFlux,
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
        for node in incoming_graph_data.nodes:
            self.add_node(node)
        for u, v, data in incoming_graph_data.edges(data=True):
            self.add_edge(u, v, **data)
        self.graph.update(incoming_graph_data.graph)

        self._prune_recycle_edges()
        self.compute_layout()

    # ── Public API ───────────────────────────────────────────────────────────


    def compute_layout(self) -> None:
        self.layout_result = self.engine.run(self)
        self.fig: Figure = self.renderer.render(self.layout_result)

    def get_figure(self) -> Figure:
        """Return the matplotlib Figure object"""
        return self.fig

    def show(self) -> None:
        """
        Render the graph, works for scripting and jupyter notebooks
        """
        if self.fig is None:
            raise RuntimeError("Call compute_layout() first")

        try:
            # Works nicely in notebooks
            from IPython import get_ipython
            if get_ipython() is not None:
                display(self.fig)
                return
        except Exception:
            pass

        # Fallback for scripts / desktop
        plt.show(block=True)

    def save(self, filename='output.png', dpi=150, format=None, **kwargs):
        """Save figure to file"""
        if self.fig is None:
            raise RuntimeError("Call compute_layout() first")
        self.fig.savefig(filename, dpi=dpi, format=format, **kwargs)


    @property
    def positions(self) -> list[Pt]:
        return self.layout_result.all_points

    # ── Graph prep ───────────────────────────────────────────────────────────

    def _prune_recycle_edges(self) -> None:
        dead = [(u, v) for u, v in self.edges
                if u.phase == -1 and v.phase == -1]
        self.remove_edges_from(dead)
        isolated = [n for n, d in self.degree() if d < 1]
        self.remove_nodes_from(isolated)
