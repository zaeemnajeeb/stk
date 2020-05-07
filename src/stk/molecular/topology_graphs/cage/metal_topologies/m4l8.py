"""
M4L8
====

"""

from ..cage import Cage
from ..vertices import _MetalVertex, _LinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M4L8(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: (0, 1, 2, 3)
        ligands: (4, 5, 6, 7, 8, 9, 10, 12)

    See :class:`.Cage` for more details and examples.

    """
    #
    # def __init__(
    #     metal_sites,
    #     linker_sites,
    #     vertex_alignments=None,
    #     reaction_factory=GenericReactionFactory(),
    #     num_processes=1,
    # ):
    #     """
    #     Initialize a :class:`.M4L4Square`.
    #
    #     """
    #
    #     print(metal_sites, linker_sites)
    #
    #     building_blocks = {**metal_sites, **linker_sites}
    #
    #     super().__init__(
    #         building_blocks,
    #         vertex_alignments=vertex_alignments,
    #         reaction_factory=reaction_factory,
    #         num_processes=num_processes,
    #     )

    _vertex_prototypes = (
        _MetalVertex(0, [1, 0, 0]),
        _MetalVertex(1, [0, 1, 0]),
        _MetalVertex(2, [-1, 0, 0]),
        _MetalVertex(3, [0, -1, 0]),

        _LinearCageVertex(4, [1, 1, 1], False),
        _LinearCageVertex(5, [1, 1, -1], False),

        _LinearCageVertex(6, [1, -1, 1], False),
        _LinearCageVertex(7, [1, -1, -1], False),

        _LinearCageVertex(8, [-1, -1, 1], False),
        _LinearCageVertex(9, [-1, -1, -1], False),

        _LinearCageVertex(10, [-1, 1, 1], False),
        _LinearCageVertex(11, [-1, 1, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[7]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[10]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[11]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[9]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[9]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[7]),
    )

    _num_windows = 2
    _num_window_types = 1
