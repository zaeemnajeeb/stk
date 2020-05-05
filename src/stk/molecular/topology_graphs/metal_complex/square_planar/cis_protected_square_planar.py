"""
Cis Protected Square Planar
===========================

"""

from ..metal_complex import MetalComplex
from ..vertices import _MetalVertex, _BiDentateLigandVertex
from ...topology_graph import Edge


class CisProtectedSquarePlanar(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    Metal building blocks with at least two functional groups are
    required for this topology graph.

    Linker building blocks with two functional groups are required
    for this topology graph.

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (
        _MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        _BiDentateLigandVertex(1, [1, 1, 0]),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[1, 0, 0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0, 1, 0],
        ),
    )
