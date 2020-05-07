"""
M3L3 Triangle
===========

"""

import numpy as np

from ..cage import Cage
from ..vertices import _BentMetalComplexCageVertex, _LinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M3L3Triangle(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with two functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: (0, 1, 2)
        ligands: (3, 4, 5)

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

    _x = np.sqrt(3)/4
    _y = 1
    _vertex_prototypes = (
        _BentMetalComplexCageVertex(0, [0, _x, 0]),
        _BentMetalComplexCageVertex(1, [_y/2, -_x, 0]),
        _BentMetalComplexCageVertex(2, [-_y/2, -_x, 0]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _LinearCageVertex.init_at_center(
            id=3,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=4,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=5,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[0]),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[3]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[4]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[0], _vertex_prototypes[5]),
    )

    _num_windows = 1
    _num_window_types = 1
