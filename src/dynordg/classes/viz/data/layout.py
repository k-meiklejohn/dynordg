from dataclasses import dataclass
from .edges import Edge, EdgeGeom, Pt, EdgeSpec
from ....constants import _GEOM_POINT_KEYS
@dataclass
class LayoutResult:
    """
    The complete, renderer-ready description of the figure.
    This is the *only* object the renderer needs to read.
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

