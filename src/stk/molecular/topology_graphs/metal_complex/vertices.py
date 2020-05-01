import numpy as np

from ..topology_graph import Vertex
from scipy.spatial.distance import euclidean


class _MetalVertex(Vertex):
    """
    Places the host in a :class:`.Complex`.

    """

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class _MonoDentateLigandVertex(Vertex):
    """
    Places monodentate ligand in a :class:`.Complex`.

    """


class _BiDentateLigandVertex(Vertex):
    """
    Places bidentate ligand in a :class:`.Complex`.

    """


class _UnsaturatedLigandVertex(Vertex):
    """
    Handles unsaturated metal site in a :class:`.Complex`.

    """
