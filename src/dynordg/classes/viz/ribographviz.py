"""
ribo_graph_vis_layout.py  –  Orchestration layer for RiboGraphFlux visualisation.

Wires together LayoutEngine (geometry) and RiboRenderer (patches → Figure).
See ribo_graph_vis_layout_engine.py for the four-phase layout pipeline and
ribo_graph_vis.py for the rendering layer.
"""

from __future__ import annotations
from ..graph import RiboGraph
from ..simulation import RiboGraphFlux
from .data import Pt
from .render import LayoutEngine, RiboRenderer
import matplotlib.pyplot as plt
from matplotlib.figure import Figure


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
        """
        Display the figure in the current environment.

        Uses IPython.display() when running inside a Jupyter notebook so the
        figure appears inline.  Falls back to plt.show(block=True) for
        scripts and desktop environments.  Raises RuntimeError if
        compute_layout() has not been called.
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

    # ── Graph prep ───────────────────────────────────────────────────────────

    def _prune_recycle_edges(self) -> None:
        """
        Remove bulk-to-bulk edges and any nodes they leave isolated.

        Edges where both u and v have phase == -1 are recycling arcs internal
        to the bulk pool.  They carry no mRNA-positional information and
        produce degenerate geometry, so they are stripped before layout runs.
        """
        dead = [(u, v) for u, v in self.edges
                if u.phase == -1 and v.phase == -1]
        self.remove_edges_from(dead)
        isolated = [n for n, d in self.degree() if d < 1]
        self.remove_nodes_from(isolated)
