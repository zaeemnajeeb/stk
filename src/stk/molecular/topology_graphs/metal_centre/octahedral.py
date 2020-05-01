"""
Octahedral
==========

"""


from .metal_centre import MetalCentre
from .vertices import _MetalVertex, _BinderVertex
from ..topology_graph import Edge


class Octahedral(MetalCentre):
    """
    Represents an octahedral metal centre topology graph.

    See :class:`.MetalCentre` for more details and examples.

    """

    _vertex_prototypes = (
        _MetalVertex(0, [0, 0, 0]),
        _BinderVertex(1, [2, 0, 0]),
        _BinderVertex(2, [0, 2, 0]),
        _BinderVertex(3, [0, 0, 2]),
        _BinderVertex(4, [-2, 0, 0]),
        _BinderVertex(5, [0, -2, 0]),
        _BinderVertex(6, [0, 0, -2])
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(4, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[0], _vertex_prototypes[6]),
    )
