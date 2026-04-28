from dataclasses import dataclass
from .edges import Edge, EdgeGeom, Pt, EdgeSpec
from ....constants import _GEOM_POINT_KEYS
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

    @property
    def all_points(self) -> list[Pt]:
        pts: list[Pt] = []
        for g in self.geoms.values():
            for attr in _GEOM_POINT_KEYS:
                pt = getattr(g, attr)
                if pt is not None:
                    pts.append(pt)
            for rect in g.helper_rects:
                pts.extend(rect)
        return pts

