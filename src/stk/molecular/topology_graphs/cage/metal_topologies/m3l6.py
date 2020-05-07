"""
M3L6
====

"""

import numpy as np

from ..cage import Cage
from ..vertices import _MetalVertex, _LinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M3L6(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with two functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    :class:`.BuildingBlock` placements:
        metals: (0, 1, 2)
        ligands: (3, 4, 5, 6, 7, 8)

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

    _R, _theta = 1, 0

    _vertex_prototypes = (
        _MetalVertex(
            id=0,
            position=[_R*np.cos(_theta), _R*np.sin(_theta), 0]
        ),
        _MetalVertex(
            id=1,
            position=[
                _R*np.cos(_theta+(4*np.pi/3)),
                _R*np.sin(_theta+(4*np.pi/3)),
                0
            ]
        ),
        _MetalVertex(
            id=2,
            position=[
                _R*np.cos(_theta+(2*np.pi/3)),
                _R*np.sin(_theta+(2*np.pi/3)),
                0
            ]
        ),

        _LinearCageVertex(
            id=3,
            position=[
                _R*np.cos((_theta+np.pi/4)),
                _R*np.sin((_theta+np.pi/4)),
                1
            ],
            use_neighbor_placement=False
        ),
        _LinearCageVertex(
            id=4,
            position=[
                _R*np.cos((_theta+1*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)),
                -1
            ],
            use_neighbor_placement=False
        ),

        _LinearCageVertex(
            id=5,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(4*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(4*np.pi/3)),
                1
            ],
            use_neighbor_placement=False
        ),
        _LinearCageVertex(
            id=6,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(4*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(4*np.pi/3)),
                -1
            ],
            use_neighbor_placement=False
        ),

        _LinearCageVertex(
            id=7,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(2*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(2*np.pi/3)),
                1
            ],
            use_neighbor_placement=False
        ),
        _LinearCageVertex(
            id=8,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(2*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(2*np.pi/3)),
                -1
            ],
            use_neighbor_placement=False
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[6]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[7]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[8]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[3]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[4]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[7]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[8]),
    )

    _num_windows = 2
    _num_window_types = 1
