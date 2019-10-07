"""
Metal Complex
=============

#. :class:`.SquarePlanar`

"""


import logging
import numpy as np

from .topology_graph import TopologyGraph, VertexData, Vertex, EdgeData
from ..reactor import Reactor

logger = logging.getLogger(__name__)


class _MetalCentreVertexData(VertexData):
    """
    Holds the data of a metal complex vertex.

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
        Initialize a :class:`_MetalCentreVertexData` instance.

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
        return _MetalCentreVertex(self)


class _MetalCentreVertex(Vertex):
    """
    Represents a vertex of a :class:`.MetalCentre`.

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

        return self._place_single_atom(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )

    def _place_single_atom(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        This method assumes the single atom should be at the position
        of the edge defined by the complex geometry.

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

        return self._assign_func_groups_to_single_atoms(
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

    def _assign_func_groups_to_single_atoms(
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

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class MetalCentre(TopologyGraph):
    """
    Represents a metal centre topology graphs.

    MetalCentre topologies are added by creating a subclass which
    defines the :attr:`vertex_data` and :attr:`edge_data` class
    attributes.

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
        # Pd2+ is 4 coordinate, so we write 4 metal functional groups.
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

    To define the metal coordination geometry, stk places functional
    groups at defined positions around the metal atoms. Here, we
    build a :class:`SquarePlanar` metal centre with nitrogen functional
    groups.

    .. code-block:: python

        # Define singular N atoms.
        m = rdkit.MolFromSmiles('N')
        m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
        n_atom = stk.BuildingBlock.init_from_rdkit_mol(
            m,
            functional_groups=['metal_bound_N'],
        )

        # Build the metal centre.
        sqpl = stk.metal_centre.SquarePlanar()
        sqpl_complex = stk.ConstructedMolecule(
            building_blocks=[metal, n_atom],
            topology_graph=sqpl,
            building_block_vertices={
                metal: tuple([sqpl.vertices[0]]),
                n_atom: sqpl.vertices[1:]
            }
        )

    This metal centre should be used as a building block for metal
    containing archiectures by converting it to a
    :class:`BuildingBlock`.

    .. code-block:: python

        metal_centre = stk.BuildingBlock.init_from_molecule(
            sqpl_complex,
            functional_groups=['metal_bound_N']
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
        Initialize a :class:`.MetalCentre`.

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
        return 1

    def _clean_up(self, mol):
        """
        Clean up of metal-centres removes the "metal" func groups.
        """
        fgs = [i for i in mol.func_groups if i.fg_type.name != 'metal']
        mol.func_groups = fgs
        super()._clean_up(mol)

    def __repr__(self):
        vertex_alignments = ', '.join(
            f'{v.id}: {v.get_aligner_edge()}'
            for v in self.vertices
        )
        return (
            f'metal_centre.{self.__class__.__name__}('
            f'vertex_alignments={{{vertex_alignments}}})'
        )


class SquarePlanar(MetalCentre):
    """
    Represents a square planar metal complex topology graph.

    See :class:`.MetalCentre` for more details and examples.

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
        _MetalCentreVertexData(0, 0, 0),
        _MetalCentreVertexData(0, 1, 0),
        _MetalCentreVertexData(0, 0, 1),
        _MetalCentreVertexData(0, -1, 0),
        _MetalCentreVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 2.0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0, 0, 2.0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[0, -2.0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[4],
            position=[0, 0, -2.0]
        ),
    )


class Paddlewheel(MetalCentre):
    """
    Represents a Paddlewheel metal complex topology graph.

    See :class:`.MetalCentre` for more details and examples.

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
        _MetalCentreVertexData(0, 0, 0),
        _MetalCentreVertexData(0, 1, 0),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 2.0, 0]
        ),
    )


class Octahedral(MetalCentre):
    """
    Represents an octahedral metal complex topology graph.

    See :class:`.MetalCentre` for more details and examples.

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
        _MetalCentreVertexData(0, 0, 0),
        _MetalCentreVertexData(0, 1, 0),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 2.0, 0]
        ),
    )
