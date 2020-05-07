"""
M8L6 Cube
=========

"""

from ..cage import Cage
from ..vertices import _NonLinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M8L6Cube(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with four functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: range(8)
        ligands: range(8, 14)

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
        _NonLinearCageVertex(0, [1, 1, 1]),
        _NonLinearCageVertex(1, [1, -1, 1]),
        _NonLinearCageVertex(2, [-1, -1, 1]),
        _NonLinearCageVertex(3, [-1, 1, 1]),
        _NonLinearCageVertex(4, [1, 1, -1]),
        _NonLinearCageVertex(5, [1, -1, -1]),
        _NonLinearCageVertex(6, [-1, -1, -1]),
        _NonLinearCageVertex(7, [-1, 1, -1]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _NonLinearCageVertex.init_at_center(
            id=8,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[2],
                _vertex_prototypes[3],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=9,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[4],
                _vertex_prototypes[5],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=10,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[3],
                _vertex_prototypes[4],
                _vertex_prototypes[7],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=11,
            vertices=(
                _vertex_prototypes[2],
                _vertex_prototypes[3],
                _vertex_prototypes[6],
                _vertex_prototypes[7],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=12,
            vertices=(
                _vertex_prototypes[4],
                _vertex_prototypes[5],
                _vertex_prototypes[6],
                _vertex_prototypes[7],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=13,
            vertices=(
                _vertex_prototypes[1],
                _vertex_prototypes[2],
                _vertex_prototypes[5],
                _vertex_prototypes[6],
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(2, _vertex_prototypes[2], _vertex_prototypes[8]),
        Edge(3, _vertex_prototypes[3], _vertex_prototypes[8]),

        Edge(4, _vertex_prototypes[4], _vertex_prototypes[9]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[9]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[9]),

        Edge(8, _vertex_prototypes[4], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[0], _vertex_prototypes[10]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[10]),
        Edge(11, _vertex_prototypes[7], _vertex_prototypes[10]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[11]),
        Edge(13, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(14, _vertex_prototypes[6], _vertex_prototypes[11]),
        Edge(15, _vertex_prototypes[7], _vertex_prototypes[11]),

        Edge(16, _vertex_prototypes[5], _vertex_prototypes[12]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(18, _vertex_prototypes[7], _vertex_prototypes[12]),
        Edge(19, _vertex_prototypes[6], _vertex_prototypes[12]),

        Edge(20, _vertex_prototypes[1], _vertex_prototypes[13]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[13]),
        Edge(22, _vertex_prototypes[6], _vertex_prototypes[13]),
        Edge(23, _vertex_prototypes[2], _vertex_prototypes[13]),
    )

    _num_windows = 4
    _num_window_types = 1
