"""
M4L4 Square
===========

"""

from ..cage import Cage
from ..vertices import _BentMetalComplexCageVertex, _LinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M4L4Square(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with two functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: (0, 1, 2, 3)
        ligands: (4, 5, 6, 7)

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
        _BentMetalComplexCageVertex(0, [1, 1, 0]),
        _BentMetalComplexCageVertex(1, [1, -1, 0]),
        _BentMetalComplexCageVertex(2, [-1, -1, 0]),
        _BentMetalComplexCageVertex(3, [-1, 1, 0]),

        _LinearCageVertex(4, [1, 0, 0], False),
        _LinearCageVertex(5, [0, -1, 0], False),
        _LinearCageVertex(6, [-1, 0, 0], False),
        _LinearCageVertex(7, [0, 1, 0], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[4]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[5]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[3], _vertex_prototypes[6]),

        Edge(6, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[7]),
    )

    _num_windows = 1
    _num_window_types = 1
