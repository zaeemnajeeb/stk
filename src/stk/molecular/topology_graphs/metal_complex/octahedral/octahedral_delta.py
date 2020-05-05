"""
Octahedral Delta
=================

"""

from ..metal_complex import MetalComplex
from ..vertices import _MetalVertex, _BiDentateLigandVertex
from ...topology_graph import Edge


class OctahedralDelta(MetalComplex):
    """
    Represents a metal complex topology graph.

    Metal building blocks with at least 6 functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    See :class:`.MetalComplex` for more details and examples.

    """

    # Use fg_assignment to ensure the desired stereochemistry is
    # achieved.
    _metal_vertex_prototypes = (
        _MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        _BiDentateLigandVertex(1, [1, 1, 0]),
        _BiDentateLigandVertex(2, [0, -1, 1]),
        _BiDentateLigandVertex(3, [-1, 0, -1]),
    )

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
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
        Edge(
            id=4,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, -1, 0],
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, 0, 1],
        ),
        Edge(
            id=5,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=[0, 0, -1],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=[-1, 0, 0],
        ),
    )
