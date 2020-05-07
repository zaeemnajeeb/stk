"""
M4L6 Tetrahedron
================

"""

import numpy as np

from ..cage import Cage
from ..vertices import _NonLinearCageVertex
from ...topology_graph import Edge
# from ...reactions import GenericReactionFactory


class M4L6Tetrahedron(Cage):
    """
    Represents a cage topology graph.

    This topology directly connects the vertices of the tetrahedron.

    Building blocks with three functional groups are required for this
    topology.

    :class:`.BuildingBlock` placements:
        metals: range(4)

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
        _NonLinearCageVertex(0, [0, 0, np.sqrt(6)/2]),
        _NonLinearCageVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(5, _vertex_prototypes[2], _vertex_prototypes[3]),
    )

    _num_windows = 4
    _num_window_types = 1
