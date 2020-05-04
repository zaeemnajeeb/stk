"""
Octahedral Delta
=================

"""

from ..complex import MetalComplex
from ..vertices import _MetalVertex, _MonoDentateLigandVertex
from ...topology_graph import Edge


class Octahedral(MetalComplex):
    """
    Represents a metal complex topology graph.

    Ligand building blocks with one functional group are required for
    this topology.

    See :class:`.MetalComplex` for more details and examples.

    """
    #
    # def map_functional_groups_to_edges(self, building_block, edges):
    #     return {
    #         fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
    #     }

    # Use fg_assignment to ensure the desired stereochemistry is
    # achieved.
    _metal_vertex_prototypes = (
        _MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        _MonoDentateLigandVertex(1, [1, 0, 0]),
        _MonoDentateLigandVertex(2, [0, 1, 0]),
        _MonoDentateLigandVertex(3, [0, 0, 1]),
        _MonoDentateLigandVertex(4, [-1, 0, 0]),
        _MonoDentateLigandVertex(5, [0, -1, 0]),
        _MonoDentateLigandVertex(6, [0, 0, -1]),
    )

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[3],
        ),
        Edge(
            id=4,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[4],
        ),
        Edge(
            id=5,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[5],
        ),
    )
