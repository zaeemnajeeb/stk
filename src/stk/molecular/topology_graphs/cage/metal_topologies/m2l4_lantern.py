"""
M2L4 Lantern
============

"""

from ..cage import Cage
from ..vertices import _MetalVertex, _LinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M2L4Lantern(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: (0, 1)
        ligands: (2, 3, 4, 5)

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
        _MetalVertex(0, [0, 0.5, 0]),
        _MetalVertex(1, [0, -0.5, 0]),

        _LinearCageVertex(2, [1, 0, 0], False),
        _LinearCageVertex(3, [0, 0, 1], False),
        _LinearCageVertex(4, [-1, 0, 0], False),
        _LinearCageVertex(5, [0, 0, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[5]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[5]),
    )

    _num_windows = 4
    _num_window_types = 1
