"""
Metal Cage
==========

#. :class:`.M2L4_Lantern`
#. :class:`.M4L8_sqpl`

"""


import logging
import numpy as np

from .topology_graph import TopologyGraph, Vertex, VertexData, EdgeData
from ..reactor import Reactor
from ...utilities import vector_angle

logger = logging.getLogger(__name__)


class _MetalCageVertexData(VertexData):
    """
    Holds the data of a metal cage vertex

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. Must match the index in
        :attr:`TopologyGraph.vertices`.

    position : :class:`numpy.ndarray`
        The position of the vertex.

    edges : :class:`list` of :class:`.EdgeData`
        The edges connected to the vertex.

    cell : :class:`numpy.ndarray`
        The unit cell in which the vertex is found.

    aligner_edge : :class:`int`
        The edge which is used to align the :class:`.BuildingBlock`
        placed on the vertex. The first :class:`.FunctionalGroup`
        in :attr:`.BuildingBlock.func_groups` is rotated such that
        it lies exactly on this :class:`.Edge`. Must be between
        ``0`` and the number of edges the vertex is connected to.

    """

    def __init__(self, x, y, z):
        """
        Initialize a :class:`_MetalCageVertexData`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        """

        self.aligner_edge = None
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            ``True`` if the clone should not be connected to any edges.

        Returns
        -------
        :class:`_CageVertexData`
            The clone.

        """

        clone = super().clone(clear_edges)
        clone.aligner_edge = self.aligner_edge
        return clone

    def get_vertex(self):
        return _MetalCageVertex(self)

    @classmethod
    def init_at_center(cls, *vertex_data):
        obj = super().init_at_center(*vertex_data)
        obj.aligner_edge = None
        return obj


class _MetalCageVertex(Vertex):
    """
    Represents a vertex of a :class:`.MetalCage`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        self._aligner_edge = data.aligner_edge
        super().__init__(data)

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            ``True`` if the clone should not be connected to any edges.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = super().clone(clear_edges)
        clone._aligner_edge = self._aligner_edge
        return clone

    def get_aligner_edge(self):
        return self._aligner_edge

    def place_building_block(self, building_block, vertices, edges):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        bb_fg_names = set((
            i.fg_type.name for i in building_block.func_groups
        ))

        if 'metal' in bb_fg_names:
            return self._place_metal_atom(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 1:
            return self._place_cap_building_block(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 2:
            return self._place_linear_building_block(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        return self._place_nonlinear_building_block(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )

    def _place_metal_atom(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """
        building_block.set_centroid(position=self._position)
        return building_block.get_position_matrix()

    def _place_cap_building_block(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """
        building_block.set_centroid(position=self._position)
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._position
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position()
        target = edge_coord - self._position
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=building_block.get_centroid()
        )

        return building_block.get_position_matrix()

    def _place_linear_building_block(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._position,
            atom_ids=building_block.get_bonder_ids()
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._position
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position()
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        target = edge_coord - self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._position
        )
        start = building_block.get_centroid_centroid_direction_vector()
        e0_coord = edges[self._edge_ids[0]].get_position()
        e1_coord = edges[self._edge_ids[1]].get_position()
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=self._position,
            axis=e0_coord-e1_coord,
            origin=self._position,
        )
        return building_block.get_position_matrix()

    def _place_nonlinear_building_block(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._position,
            atom_ids=building_block.get_bonder_ids()
        )
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        edge_normal = self._get_edge_plane_normal(
            reference=self._get_edge_centroid(
                centroid_edges=connected_edges,
                vertices=vertices
            ),
            plane_edges=connected_edges,
            vertices=vertices
        )
        building_block.apply_rotation_between_vectors(
            start=building_block.get_bonder_plane_normal(),
            target=edge_normal,
            origin=self._position
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_bonder_centroid - self._position
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=edge_normal,
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in :attr:`edges` it
            is assigned to.

        """

        bb_fg_names = set((
            i.fg_type.name for i in building_block.func_groups
        ))
        if 'metal' in bb_fg_names:
            return self._assign_func_groups_to_metal_atom(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 1:
            return self._assign_func_groups_to_cap_edges(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 2:
            return self._assign_func_groups_to_linear_edges(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        return self._assign_func_groups_to_nonlinear_edges(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )

    def after_assign_func_groups_to_edges(
        self,
        building_block,
        func_groups,
        vertices,
        edges
    ):
        """
        Perform operations after functional groups have been assigned.

        This method is always executed serially. It is often useful
        when data needs to be transferred between vertices, which
        have been processed independently, in parallel.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        func_groups : :class:`tuple` of :class:`.FunctionalGroup`
            The functional group clones added to the constructed
            molecule.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        None : :class:`NoneType`

        """

        return super().after_assign_func_groups_to_edges(
            building_block=building_block,
            func_groups=func_groups,
            vertices=vertices,
            edges=edges
        )

    def _assign_func_groups_to_metal_atom(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(building_block, edges)
            ))
        }

    def _assign_func_groups_to_cap_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(building_block, edges)
            ))
        }

    def _assign_func_groups_to_linear_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(building_block, edges)
            ))
        }

    def _get_fg0_distance(self, building_block, edges):
        fg_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )

        def distance(edge_id):
            displacement = edges[edge_id].get_position() - fg_coord
            return np.linalg.norm(displacement)

        return distance

    def _assign_func_groups_to_nonlinear_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        # The idea is to order the functional groups in building_block
        # by their angle from func_groups[0] and the bonder centroid,
        #  going in the clockwise direction.
        #
        # The edges are also ordered by their angle from aligner_edge
        # and the edge centroid going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg0_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.get_bonder_ids()
        )
        fg0_direction = fg0_coord-bonder_centroid
        axis = np.cross(
            fg0_direction,
            building_block.get_bonder_plane_normal()
        )
        func_groups = sorted(
            range(len(building_block.func_groups)),
            key=self._get_func_group_angle(
                building_block=building_block,
                fg0_direction=fg0_direction,
                bonder_centroid=bonder_centroid,
                axis=axis
            )
        )
        assignments = {}
        edge_ids = sorted(
            self._edge_ids,
            key=self._get_edge_angle(axis, vertices, edges)
        )
        for edge_id, fg_id in zip(edge_ids, func_groups):
            assignments[fg_id] = edge_id
        return assignments

    @staticmethod
    def _get_func_group_angle(
        building_block,
        fg0_direction,
        bonder_centroid,
        axis
    ):

        def angle(fg_id):
            func_group = building_block.func_groups[fg_id]
            coord = building_block.get_centroid(
                atom_ids=func_group.get_bonder_ids()
            )
            fg_direction = coord-bonder_centroid
            theta = vector_angle(fg0_direction, fg_direction)

            projection = fg_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def _get_edge_angle(self, axis, vertices, edges):
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        aligner_edge_coord = aligner_edge.get_position()
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        edge_centroid = self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        # This axis is used to figure out the clockwise direction.
        aligner_edge_direction = aligner_edge_coord - edge_centroid

        def angle(edge_id):
            coord = edges[edge_id].get_position()
            edge_direction = coord - edge_centroid
            theta = vector_angle(
                vector1=edge_direction,
                vector2=aligner_edge_direction
            )

            projection = edge_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class MetalCage(TopologyGraph):
    """
    Represents metal cage topology graphs.

    MetalCage topologies are added by creating a subclass which
    defines the :attr:`vertex_data` and :attr:`edge_data` class
    attributes.

    This class is modelled after :class:`Cage` and its subclasses.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------

    Construction of topologies that contain metals requires some extra
    manipulation compared to purely organic systems because SMILES and
    RDKit are not properly defined to handle metal chemistry. To define
    a :class:`.BuildingBlock` containing metal atoms and functional
    groups in stk requires something like below to avoid RDKit
    problems upon initialization. Note that here the functional groups
    do not define the metal centre coordination geometry, but does
    define the number of coordination sites per metal.

    .. code-block:: python

        import stk
        from rdkit.Chem import AllChem as rdkit
        import numpy as np

        # Define an stk BuildingBlock with no functional groups and a
        # single metal (Pd with 2+ charge) atom.
        m = rdkit.MolFromSmiles('[Pd+2]')
        m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
        metal = stk.BuildingBlock.init_from_rdkit_mol(
            m,
            functional_groups=None,
        )

        # Manually set atom position if required. RDKit will default
        # to [0, 0, 0].
        metal.set_position_matrix(np.array([[0.0, 0.0, 0.0]]))

        # Manually set functional group information for the metal
        # centre. The key of this dictionary is the fg.id.
        # Pd2+ is 4 coordinate, so we write 4 metal functiongal groups.
        metal_coord_info = {
            0: {
                'atom_ids': [0],
                'bonder_ids': [0],
                'deleter_ids': [None]
            },
            1: {
                'atom_ids': [0],
                'bonder_ids': [0],
                'deleter_ids': [None]
            },
            2: {
                'atom_ids': [0],
                'bonder_ids': [0],
                'deleter_ids': [None]
            },
            3: {
                'atom_ids': [0],
                'bonder_ids': [0],
                'deleter_ids': [None]
            },
        }
        metal = stk.assign_metal_fgs(
            building_block=metal,
            coordination_info=metal_coord_info
        )

    The metal centre coordination geometry is defined by the topology
    graph used for construction; in this case the metal centre geometry
    is defined by the cage topology used. Ligands that bind to
    metal centres are given functional groups that are designated for
    metal interaction with 'metal' in their name.

    .. code-block:: python

        ligand = stk.BuildingBlock(
            'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
            functional_groups=['pyridine_N_metal']
        )

        m2l4_lantern = stk.metal_cage.M2L4_Lantern()
        lantern = stk.ConstructedMolecule(
            building_blocks=[metal, ligand],
            topology_graph=m2l4_lantern,
            building_block_vertices={
                metal: m2l4_lantern.vertices[0:2],
                ligand: m2l4_lantern.vertices[2:]
            }
        )

    It is crucial to save metal-containing molecules using the
    :meth:`.dump` option to a :class:`dict` because RDKit may fail to
    load the molecule back from a .mol file. Note that the metal
    functional groups assigned to the building blocks will be lost.

    .. code-block:: python

        # Dump ConstructedMolecule.
        lantern.dump('m2l4_lantern.json')
        # Load in ConstructedMolecule.
        loaded = stk.ConstructedMolecule.load('m2l4_lantern.json')

    As with other cage toplogies, different mixtures of ligands can be
    used in any topology using the `building_block_vertices` variable.

    .. code-block:: python

        ligand1 = stk.BuildingBlock(
            'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
            functional_groups=['pyridine_N_metal']
        )
        ligand2 = stk.BuildingBlock(
            'COc1c(OC)c2ccc(-c3ccncc3)cc2c2cc(-c3ccncc3)ccc12',
            functional_groups=['pyridine_N_metal']
        )
        m2l4_lantern = stk.metal_cage.M2L4_Lantern()
        hetero_lantern = stk.ConstructedMolecule(
            building_blocks=[metal, ligand1, ligand2],
            topology_graph=m2l4_lantern,
            building_block_vertices={
                metal: m2l4_lantern.vertices[0:2],
                ligand1: m2l4_lantern.vertices[2:4],
                ligand2: m2l4_lantern.vertices[4:]
            }
        )

    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertex_data):
            vertex.id = i
        for i, edge in enumerate(cls.edge_data):
            edge.id = i
        return super().__init_subclass__(**kwargs)

    def __init__(self, vertex_alignments=None, num_processes=1):
        """
        Initialize a :class:`.MetalCage`.

        Parameters
        ----------
        vertex_alignments : :class:`dict`, optional
            A mapping from the :attr:`.Vertex.id` of a :class:`.Vertex`
            :attr:`vertices` to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is refered to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if vertex_alignments is None:
            vertex_alignments = {}

        vertex_data = {
            data: data.clone(True) for data in self.vertex_data
        }
        for vertex in vertex_data.values():
            vertex.aligner_edge = vertex_alignments.get(vertex.id, 0)
        edge_data = tuple(
            edge.clone(vertex_data)
            for edge in self.edge_data
        )
        vertex_types = sorted(
            {len(v.edges) for v in vertex_data},
            reverse=True
        )
        super().__init__(
            vertex_data=vertex_data.values(),
            edge_data=edge_data,
            construction_stages=tuple(
                lambda vertex, vertex_type=vt:
                    vertex.get_num_edges() == vertex_type
                for vt in vertex_types
            ),
            num_processes=num_processes
        )

    def construct(self, mol):
        """
        Construct a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance which needs to
            be constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        scale = self._get_scale(mol)
        vertices = tuple(self._get_vertex_clones(mol, scale))
        edges = tuple(self._get_edge_clones(scale))

        self._prepare(mol)
        self._place_building_blocks(mol, vertices, edges)

        vertices, edges = self._before_react(mol, vertices, edges)
        reactor = Reactor(mol)
        for edge in edges:
            reactor.add_reaction(
                func_groups=edge.get_func_groups(),
                periodicity=tuple(edge.get_periodicity())
            )
        reactor.finalize()

        # Assign the information of Edge and Vertex classes and their
        # assignments.
        self.vertex_edge_assignments = {}
        for vertex in vertices:
            v_id = vertex.id
            vertex_edges = []
            for edge in edges:
                if edge in vertex_edges:
                    continue
                if v_id in edge._vertex_ids:
                    vertex_edges.append(edge)
            vertex.edges = vertex_edges
            self.vertex_edge_assignments[vertex.id] = vertex

        self._clean_up(mol)

    def _get_scale(self, mol):
        """
        Get the scale used for the positions of :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        """
        organic_bbs = [
            i for i in mol.building_block_vertices
            if 'metal' not in set((
                j.fg_type.name for j in i.func_groups
            ))
        ]
        if organic_bbs:
            return max(bb.get_maximum_diameter() for bb in organic_bbs)
        else:
            # No organic building blocks.
            return 1

    def __repr__(self):
        vertex_alignments = ', '.join(
            f'{v.id}: {v.get_aligner_edge()}'
            for v in self.vertices
        )
        return (
            f'metal_cage.{self.__class__.__name__}('
            f'vertex_alignments={{{vertex_alignments}}})'
        )


class M2L4_Lantern(MetalCage):
    """
    Represents a M2L4 lantern topology graph with square planar metal.

    See :class:`.MetalCage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, -1, 0),
        _MetalCageVertexData(1, 0, 0),
        _MetalCageVertexData(0, 0, 1),
        _MetalCageVertexData(-1, 0, 0),
        _MetalCageVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.2, 1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[0, 1, 0.2]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[4],
            position=[-0.2, 1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[5],
            position=[0, 1, -0.2]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[2],
            position=[0.2, -1, 0]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[3],
            position=[0, -1, 0.2]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[4],
            position=[-0.2, -1, 0]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[5],
            position=[0, -1, -0.2]
        ),
    )

    num_windows = 4
    num_window_types = 1


class M4L8_sqpl(MetalCage):
    """
    Represents a M4L8 topology graph with square planar metal.

    See :class:`.MetalCage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalCageVertexData(1, 0, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(-1, 0, 0),
        _MetalCageVertexData(0, -1, 0),

        _MetalCageVertexData(1, 1, 1),
        _MetalCageVertexData(1, 1, -1),

        _MetalCageVertexData(1, -1, 1),
        _MetalCageVertexData(1, -1, -1),

        _MetalCageVertexData(-1, -1, 1),
        _MetalCageVertexData(-1, -1, -1),

        _MetalCageVertexData(-1, 1, 1),
        _MetalCageVertexData(-1, 1, -1),

    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[4],
            position=[1, 0.2, 0.2]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[5],
            position=[1, 0.2, -0.2]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[6],
            position=[1, -0.2, 0.2]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[7],
            position=[1, -0.2, -0.2]
        ),

        EdgeData(
            vertex_data[1],
            vertex_data[4],
            position=[0.2, 1, 0.2]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[5],
            position=[0.2, 1, -0.2]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[10],
            position=[-0.2, 1, 0.2]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[11],
            position=[-0.2, 1, -0.2]
        ),

        EdgeData(
            vertex_data[2],
            vertex_data[10],
            position=[-1, 0.2, 0.2]
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[11],
            position=[-1, 0.2, -0.2]
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[8],
            position=[-1, -0.2, 0.2]
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[9],
            position=[-1, -0.2, -0.2]
        ),

        EdgeData(
            vertex_data[3],
            vertex_data[8],
            position=[-0.2, -1, 0.2]
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[9],
            position=[-0.2, -1, -0.2]
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[6],
            position=[0.2, -1, 0.2]
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[7],
            position=[0.2, -1, -0.2]
        ),
    )

    num_windows = 6
    num_window_types = 2


class M12L24_sqpl(MetalCage):
    """
    Represents a M12L24 topology graph with square planar metal.

    See :class:`.MetalCage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, 1, 0),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.2, 1, 0]
        ),
    )

    num_windows = 14
    num_window_types = 2


class M6L12_cube(MetalCage):
    """
    Represents a M6L12 cube topology graph with square planar metal.

    See :class:`.MetalCage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalCageVertexData(1, 0, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(-1, 0, 0),
        _MetalCageVertexData(0, -1, 0),
        _MetalCageVertexData(0, 0, 1),
        _MetalCageVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.2, 1, 0]
        ),
    )

    num_windows = 8
    num_window_types = 1


class M24L48_sqpl(MetalCage):
    """
    Represents a M24L48 topology graph with square planar metal.

    See :class:`.MetalCage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalCageVertexData(0.7071, 0.7071, 0.0),
        _MetalCageVertexData(0.7071, 0.7071, 1.0),
        _MetalCageVertexData(0.2929, 0.7071, 0.0),
        _MetalCageVertexData(0.2929, 0.7071, 1.0),
        _MetalCageVertexData(0.2929, 0.2929, 0.0),
        _MetalCageVertexData(0.2929, 0.2929, 1.0),
        _MetalCageVertexData(0.7071, 0.2929, 0.0),
        _MetalCageVertexData(0.7071, 0.2929, 1.0),
        _MetalCageVertexData(0.0, 0.7071, 0.7071),
        _MetalCageVertexData(1.0, 0.7071, 0.7071),
        _MetalCageVertexData(0.2929, 0.0, 0.7071),
        _MetalCageVertexData(0.2929, 1.0, 0.7071),
        _MetalCageVertexData(0.0, 0.2929, 0.7071),
        _MetalCageVertexData(1.0, 0.2929, 0.7071),
        _MetalCageVertexData(0.7071, 0.0, 0.7071),
        _MetalCageVertexData(0.7071, 1.0, 0.7071),
        _MetalCageVertexData(0.0, 0.2929, 0.2929),
        _MetalCageVertexData(1.0, 0.2929, 0.2929),
        _MetalCageVertexData(0.7071, 0.0, 0.2929),
        _MetalCageVertexData(0.7071, 1.0, 0.2929),
        _MetalCageVertexData(0.0, 0.7071, 0.2929),
        _MetalCageVertexData(1.0, 0.7071, 0.2929),
        _MetalCageVertexData(0.2929, 0.0, 0.2929),
        _MetalCageVertexData(0.2929, 1.0, 0.2929),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.2, 1, 0]
        ),
    )

    num_windows = 26
    num_window_types = 2
