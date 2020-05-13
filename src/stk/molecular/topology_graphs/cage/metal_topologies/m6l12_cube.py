"""
M6L12 Cube
==========

"""

import numpy as np

from ..cage import Cage
from ..vertices import _MetalVertex, _LinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M6L12Cube(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: (0, 1, 2, 3, 4, 5)
        ligands: (6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17)

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

    _x = np.sqrt(2)
    _vertex_prototypes = (
        _MetalVertex(0, [_x, 0, 0]),
        _MetalVertex(1, [0, _x, 0]),
        _MetalVertex(2, [-_x, 0, 0]),
        _MetalVertex(3, [0, -_x, 0]),
        _MetalVertex(4, [0, 0, _x]),
        _MetalVertex(5, [0, 0, -_x]),

        _LinearCageVertex(6, [1, 1, 0], False),
        _LinearCageVertex(7, [1, -1, 0], False),
        _LinearCageVertex(8, [1, 0, 1], False),
        _LinearCageVertex(9, [1, 0, -1], False),
        _LinearCageVertex(10, [-1, 1, 0], False),
        _LinearCageVertex(11, [-1, -1, 0], False),
        _LinearCageVertex(12, [-1, 0, 1], False),
        _LinearCageVertex(13, [-1, 0, -1], False),
        _LinearCageVertex(14, [0, 1, 1], False),
        _LinearCageVertex(15, [0, 1, -1], False),
        _LinearCageVertex(16, [0, -1, 1], False),
        _LinearCageVertex(17, [0, -1, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[7]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[9]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[10]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[14]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[15]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[12]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[13]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[11]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[16]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[17]),

        Edge(16, _vertex_prototypes[4], _vertex_prototypes[8]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(18, _vertex_prototypes[4], _vertex_prototypes[14]),
        Edge(19, _vertex_prototypes[4], _vertex_prototypes[16]),

        Edge(20, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[13]),
        Edge(22, _vertex_prototypes[5], _vertex_prototypes[15]),
        Edge(23, _vertex_prototypes[5], _vertex_prototypes[17]),
    )

    _num_windows = 8
    _num_window_types = 1