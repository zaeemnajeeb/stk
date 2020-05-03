import numpy as np
from scipy.spatial.distance import euclidean

from stk.utilities import (
    get_acute_vector,
    get_plane_normal,
)
from ..utilities import _FunctionalGroupSorter, _EdgeSorter
from ..topology_graph import Vertex


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
        edges = sorted(edges, key=lambda i: i.get_id())
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

    def __init__(
        self,
        id,
        position,
        aligner_edge=0,
    ):
        """
        Initialize a :class:`._CageVertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        """

        self._aligner_edge = aligner_edge
        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        return clone

    def _with_aligner_edge(self, aligner_edge):
        """
        Modify the instance.

        """

        self._aligner_edge = aligner_edge
        return self

    def with_aligner_edge(self, aligner_edge):
        """
        Return a clone with a different `aligner_edge`.

        Parameters
        ----------
        aligner_edge : :class:`int`
            The aligner edge of the clone.

        Returns
        -------
        :class:`._CageVertex`
            The clone. Has the same type as the original instance.

        """

        return self.clone()._with_aligner_edge(aligner_edge)

    def place_building_block(self, building_block, edges):
        # Translate building block to vertex position.
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )

        # Align vector between 2 edges with vector between centroid of
        # placers in 2 FGs.
        fg0, fg1 = building_block.get_functional_groups()
        fg0_position = building_block.get_centroid(
            atom_ids=fg0.get_placer_ids(),
        )
        fg1_position = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        start = fg1_position - fg0_position
        # Vector between connected edges.
        c_edge_positions = [
            i.get_position() for i in edges
        ]
        target = c_edge_positions[1] - c_edge_positions[0]
        building_block = building_block.with_rotation_between_vectors(
            start=start,
            target=target,
            origin=building_block.get_centroid(),
        )

        # Align vector between edge-self.position with vector between
        # placer centroid and deleter centroid.
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        deleter_ids = list((
            j.get_id()
            for i in building_block.get_functional_groups()
            for j in i.get_deleters()
        ))
        deleter_centroid = building_block.get_centroid(
            atom_ids=deleter_ids
        )
        start = placer_centroid - deleter_centroid
        target = self._position - edge_centroid
        building_block = building_block.with_rotation_between_vectors(
            start=start,
            target=target,
            origin=building_block.get_centroid(),
        )

        # Translate building block to vertex position.
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        return building_block.get_position_matrix()
        #
        # print('aaaaaaaaa22222222aaa')
        # # return building_block.get_position_matrix()
        # # Align normal of placer plane and edge plane.
        # edge_centroid = (
        #     sum(edge.get_position() for edge in edges) / len(edges)
        # )
        # core_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_core_atom_ids(),
        # )
        # placer_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_placer_ids(),
        # )
        # deleter_ids = list((
        #     j.get_id()
        #     for i in building_block.get_functional_groups()
        #     for j in i.get_deleters()
        # ))
        # deleter_centroid = building_block.get_centroid(
        #     atom_ids=deleter_ids
        # )
        # print(list(edge.get_position() for edge in edges))
        # print(edge_centroid)
        # target_points = [edge.get_position() for edge in edges]
        # target_points.append(self._position)
        # target_points = np.array(target_points)
        # print(target_points)
        # target_normal = get_acute_vector(
        #     reference=self._position - edge_centroid,
        #     vector=get_plane_normal(points=target_points),
        # )
        # print(list(
        #     i
        #     for i in building_block.get_core_atom_ids()
        # ))
        # print(list(
        #     i
        #     for i in building_block.get_placer_ids()
        # ))
        #
        # print('en', target_normal)
        # print('c', core_centroid)
        # print('p', placer_centroid)
        # print('d', deleter_centroid)
        # print('dp', deleter_centroid - placer_centroid)
        # start = get_acute_vector(
        #     reference=deleter_centroid - placer_centroid,
        #     vector=building_block.get_plane_normal(
        #         atom_ids=building_block.get_placer_ids(),
        #     ),
        # )
        # print('s', start)
        # # return building_block.get_position_matrix()
        # building_block = building_block.with_rotation_between_vectors(
        #     start=start,
        #     target=target_normal,
        #     origin=self._position,
        # )
        # core_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_core_atom_ids(),
        # )
        # placer_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_placer_ids(),
        # )
        # start = get_acute_vector(
        #     reference=placer_centroid - core_centroid,
        #     vector=building_block.get_plane_normal(
        #         atom_ids=building_block.get_placer_ids(),
        #     ),
        # )
        # print('ns', start)
        # print('============')
        # return building_block.get_position_matrix()

        #
        # Align vector 1 and 2.
        # Vector 1: vector between deleter centroid and bonder centroid
        # Vector 2: vector between edge centrod and vertex position.
        # deleter_ids = list((
        #     j.get_id()
        #     for i in building_block.get_functional_groups()
        #     for j in i.get_deleters()
        # ))
        # deleter_centroid = building_block.get_centroid(
        #     atom_ids=deleter_ids
        # )
        # placer_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_placer_ids(),
        # )
        # start = deleter_centroid - placer_centroid
        # edge_centroid = (
        #     sum(edge.get_position() for edge in edges) / len(edges)
        # )
        # target = edge_centroid - self._position
        # building_block = building_block.with_rotation_between_vectors(
        #     start=start,
        #     target=target,
        #     origin=building_block.get_centroid(),
        # )
        # return building_block.get_position_matrix()
        #
        # # Rotate vector1 onto vector2:
        # print('--')
        # # Vector1 : bonder centroid - edge centroid
        # # Vector2 : self._position - edge centroid
        # placer_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_placer_ids(),
        # )
        # start = placer_centroid - edge_centroid
        # print('init', placer_centroid, edge_centroid)
        # edge_normal = get_acute_vector(
        #     reference=edge_centroid,
        #     vector=get_plane_normal(
        #         points=np.array([
        #             edge.get_position() for edge in edges
        #         ]),
        #     ),
        # )
        # print(edge_normal)
        # target = self._position - edge_centroid
        # print(self._position, edge_centroid)
        # print('ts', target, start)
        # building_block = building_block.with_rotation_between_vectors(
        #     start=start,
        #     target=target,
        #     origin=building_block.get_centroid(),
        # )
        # placer_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_placer_ids(),
        # )
        # print('rot', placer_centroid, edge_centroid, placer_centroid - edge_centroid)
        # print('=======')
        # return building_block.get_position_matrix()

        # Rotate building block about edge centroid to minimize the
        # angles between
        # FG_ binder atoms - edge_i position - binder centroid.

        #
        # fg_bonder_centroid = building_block.get_centroid(
        #     atom_ids=next(
        #         building_block.get_functional_groups()
        #     ).get_placer_ids(),
        # )
        # start = fg_bonder_centroid - building_block.get_centroid()
        # edge_centroid = (
        #     sum(edge.get_position() for edge in edges) / len(edges)
        # )
        # target = edge_centroid - self._position
        # edge_centroid = (
        #     sum(edge.get_position() for edge in edges) / len(edges)
        # )
        # edge_normal = get_acute_vector(
        #     reference=edge_centroid,
        #     vector=get_plane_normal(
        #         points=np.array([
        #             edge.get_position() for edge in edges
        #         ]),
        #     ),
        # )
        # core_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_core_atom_ids(),
        # )
        # placer_centroid = building_block.get_centroid(
        #     atom_ids=building_block.get_placer_ids(),
        # )
        # building_block = building_block.with_rotation_between_vectors(
        #     start=get_acute_vector(
        #         reference=core_centroid - placer_centroid,
        #         vector=building_block.get_plane_normal(
        #             atom_ids=building_block.get_placer_ids(),
        #         ),
        #     ),
        #     target=edge_normal,
        #     origin=self._position,
        # )
        # fg_bonder_centroid = building_block.get_centroid(
        #     atom_ids=next(
        #         building_block.get_functional_groups()
        #     ).get_placer_ids(),
        # )
        # edge_position = edges[self._aligner_edge].get_position()
        # building_block = building_block.with_rotation_to_minimize_angle(
        #     start=fg_bonder_centroid - self._position,
        #     target=edge_position - edge_centroid,
        #     axis=edge_normal,
        #     origin=self._position,
        # )
        # building_block = building_block.with_rotation_between_vectors(
        #     start=start,
        #     target=target,
        #     origin=building_block.get_centroid(
        #         atom_ids=building_block.get_placer_ids()
        #     ),
        # )
        # return building_block.with_rotation_between_vectors(
        #     start=start,
        #     target=target,
        #     origin=building_block.get_centroid(),
        # ).get_position_matrix()
        # return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }
