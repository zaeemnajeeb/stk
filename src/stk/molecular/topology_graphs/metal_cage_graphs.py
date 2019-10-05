"""
Defines metal-cage topology graphs.

"""


import logging
import numpy as np

from .topology_graph import EdgeData
from .metal_cage import MetalCage, _MetalCageVertexData

logger = logging.getLogger(__name__)


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

        _MetalCageVertexData(1, 1, 0),
        _MetalCageVertexData(1, -1, 0),
        _MetalCageVertexData(1, 0, 1),
        _MetalCageVertexData(1, 0, -1),
        _MetalCageVertexData(-1, 1, 0),
        _MetalCageVertexData(-1, -1, 0),
        _MetalCageVertexData(-1, 0, 1),
        _MetalCageVertexData(-1, 0, -1),
        _MetalCageVertexData(0, 1, 1),
        _MetalCageVertexData(0, 1, -1),
        _MetalCageVertexData(0, -1, 1),
        _MetalCageVertexData(0, -1, -1),


    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[6],
            position=[1, 0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[7],
            position=[1, -0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[8],
            position=[1, 0, 0.2]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[9],
            position=[1, 0, -0.2]
        ),

        EdgeData(
            vertex_data[1],
            vertex_data[6],
            position=[0.2, 1, 0]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[10],
            position=[-0.2, 1, 0]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[14],
            position=[0, 1, 0.2]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[15],
            position=[0, 1, -0.2]
        ),

        EdgeData(
            vertex_data[2],
            vertex_data[10],
            position=[-1, 0.2, 0]
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[11],
            position=[-1, -0.2, 0]
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[12],
            position=[-1, 0, 0.2]
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[13],
            position=[-1, 0, -0.2]
        ),

        EdgeData(
            vertex_data[3],
            vertex_data[7],
            position=[0.2, -1, 0]
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[11],
            position=[-0.2, -1, 0]
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[16],
            position=[0, -1, 0.2]
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[17],
            position=[0, -1, -0.2]
        ),

        EdgeData(
            vertex_data[4],
            vertex_data[8],
            position=[0.2, 0, 1]
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[12],
            position=[-0.2, 0, 1]
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[14],
            position=[0, 0.2, 1]
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[16],
            position=[0, -0.2, 1]
        ),

        EdgeData(
            vertex_data[5],
            vertex_data[9],
            position=[0.2, 0, -1]
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[13],
            position=[-0.2, 0, -1]
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[15],
            position=[0, 0.2, -1]
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[17],
            position=[0, -0.2, -1]
        ),
    )

    num_windows = 8
    num_window_types = 1


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

    def metal_bonder_vector(metal_v, lig_v):
        """
        Define metal bonding vector from vertex positions.

        Parameters
        ----------
        metal_v : :class:`._MetalCageVertexData`
            A metal vertex.

        lig_v : :class:`._MetalCageVertexData`
            A ligand vertex bonded to the metal.

        Returns
        -------
        bv : :class:`float`
            The position of the edge between the metal and ligand
            vertex, defined to be 0.2 angstrom from metal_v in the
            direction of the ligand vertex.

        """
        vector = lig_v.position - metal_v.position
        bv = (vector / np.linalg.norm(vector)) * 0.2
        bv = bv + metal_v.position
        return bv

    vertex_data = (
        _MetalCageVertexData(1, 0, 0),
        _MetalCageVertexData(-1, 0, 0),
        _MetalCageVertexData(0, 1, 0),
        _MetalCageVertexData(0, -1, 0),
        _MetalCageVertexData(0.5, 0.5, 0.707),
        _MetalCageVertexData(0.5, -0.5, 0.707),
        _MetalCageVertexData(-0.5, 0.5, 0.707),
        _MetalCageVertexData(-0.5, -0.5, 0.707),
        _MetalCageVertexData(0.5, 0.5, -0.707),
        _MetalCageVertexData(0.5, -0.5, -0.707),
        _MetalCageVertexData(-0.5, 0.5, -0.707),
        _MetalCageVertexData(-0.5, -0.5, -0.707),

        _MetalCageVertexData(1, 0.35, 0.35),
        _MetalCageVertexData(1, 0.35, -0.35),
        _MetalCageVertexData(1, -0.35, 0.35),
        _MetalCageVertexData(1, -0.35, -0.35),

        _MetalCageVertexData(-1, 0.35, 0.35),
        _MetalCageVertexData(-1, 0.35, -0.35),
        _MetalCageVertexData(-1, -0.35, 0.35),
        _MetalCageVertexData(-1, -0.35, -0.35),

        _MetalCageVertexData(0.35, 1, 0.35),
        _MetalCageVertexData(0.35, 1, -0.35),
        _MetalCageVertexData(-0.35, 1, 0.35),
        _MetalCageVertexData(-0.35, 1, -0.35),

        _MetalCageVertexData(0.35, -1, 0.35),
        _MetalCageVertexData(0.35, -1, -0.35),
        _MetalCageVertexData(-0.35, -1, 0.35),
        _MetalCageVertexData(-0.35, -1, -0.35),

        _MetalCageVertexData(0.5, 0, 0.707),
        _MetalCageVertexData(-0.5, 0, 0.707),
        _MetalCageVertexData(0, 0.5, 0.707),
        _MetalCageVertexData(0, -0.5, 0.707),
        _MetalCageVertexData(0.5, 0, -0.707),
        _MetalCageVertexData(-0.5, 0, -0.707),
        _MetalCageVertexData(0, 0.5, -0.707),
        _MetalCageVertexData(0, -0.5, -0.707),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[12],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[12]
            )
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[13],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[13]
            )
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[14],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[14]
            )
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[15],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[15]
            )
        ),

        EdgeData(
            vertex_data[1],
            vertex_data[16],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[16]
            )
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[17],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[17]
            )
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[18],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[18]
            )
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[19],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[19]
            )
        ),

        EdgeData(
            vertex_data[2],
            vertex_data[20],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[20]
            )
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[21],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[21]
            )
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[22],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[22]
            )
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[23],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[23]
            )
        ),

        EdgeData(
            vertex_data[3],
            vertex_data[24],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[24]
            )
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[25],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[25]
            )
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[26],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[26]
            )
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[27],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[27]
            )
        ),

        EdgeData(
            vertex_data[4],
            vertex_data[28],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[28]
            )
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[30],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[30]
            )
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[12],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[12]
            )
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[20],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[20]
            )
        ),

        EdgeData(
            vertex_data[5],
            vertex_data[14],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[14]
            )
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[24],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[24]
            )
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[28],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[28]
            )
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[31],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[31]
            )
        ),

        EdgeData(
            vertex_data[6],
            vertex_data[16],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[16]
            )
        ),
        EdgeData(
            vertex_data[6],
            vertex_data[29],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[29]
            )
        ),
        EdgeData(
            vertex_data[6],
            vertex_data[30],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[30]
            )
        ),
        EdgeData(
            vertex_data[6],
            vertex_data[22],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[22]
            )
        ),

        EdgeData(
            vertex_data[7],
            vertex_data[18],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[18]
            )
        ),
        EdgeData(
            vertex_data[7],
            vertex_data[26],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[26]
            )
        ),
        EdgeData(
            vertex_data[7],
            vertex_data[31],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[31]
            )
        ),
        EdgeData(
            vertex_data[7],
            vertex_data[29],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[29]
            )
        ),

        EdgeData(
            vertex_data[8],
            vertex_data[13],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[13]
            )
        ),
        EdgeData(
            vertex_data[8],
            vertex_data[32],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[32]
            )
        ),
        EdgeData(
            vertex_data[8],
            vertex_data[34],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[34]
            )
        ),
        EdgeData(
            vertex_data[8],
            vertex_data[21],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[21]
            )
        ),

        EdgeData(
            vertex_data[9],
            vertex_data[15],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[15]
            )
        ),
        EdgeData(
            vertex_data[9],
            vertex_data[32],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[32]
            )
        ),
        EdgeData(
            vertex_data[9],
            vertex_data[35],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[35]
            )
        ),
        EdgeData(
            vertex_data[9],
            vertex_data[25],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[25]
            )
        ),

        EdgeData(
            vertex_data[10],
            vertex_data[17],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[17]
            )
        ),
        EdgeData(
            vertex_data[10],
            vertex_data[23],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[23]
            )
        ),
        EdgeData(
            vertex_data[10],
            vertex_data[34],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[34]
            )
        ),
        EdgeData(
            vertex_data[10],
            vertex_data[33],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[33]
            )
        ),

        EdgeData(
            vertex_data[11],
            vertex_data[19],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[19]
            )
        ),
        EdgeData(
            vertex_data[11],
            vertex_data[33],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[33]
            )
        ),
        EdgeData(
            vertex_data[11],
            vertex_data[27],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[27]
            )
        ),
        EdgeData(
            vertex_data[11],
            vertex_data[35],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[35]
            )
        ),
    )

    num_windows = 14
    num_window_types = 2


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

    def metal_bonder_vector(metal_v, lig_v):
        """
        Define metal bonding vector from vertex positions.

        Parameters
        ----------
        metal_v : :class:`._MetalCageVertexData`
            A metal vertex.

        lig_v : :class:`._MetalCageVertexData`
            A ligand vertex bonded to the metal.

        Returns
        -------
        bv : :class:`float`
            The position of the edge between the metal and ligand
            vertex, defined to be 0.2 angstrom from metal_v in the
            direction of the ligand vertex.

        """
        vector = lig_v.position - metal_v.position
        bv = (vector / np.linalg.norm(vector)) * 0.2
        bv = bv + metal_v.position
        return bv

    coord1 = 0.414
    coord2 = -0.414
    coord3 = coord1 - (coord1-coord2)/2
    coord4 = -1
    coord5 = 1
    coord6 = coord5 - (coord5-coord1)/2

    vertex_data = (
        _MetalCageVertexData(coord1, coord1, coord4),
        _MetalCageVertexData(coord1, coord1, coord5),
        _MetalCageVertexData(coord1, coord4, coord2),
        _MetalCageVertexData(coord1, coord5, coord2),
        _MetalCageVertexData(coord1, coord2, coord4),
        _MetalCageVertexData(coord1, coord2, coord5),
        _MetalCageVertexData(coord1, coord4, coord1),
        _MetalCageVertexData(coord1, coord5, coord1),
        _MetalCageVertexData(coord2, coord1, coord4),
        _MetalCageVertexData(coord2, coord1, coord5),
        _MetalCageVertexData(coord2, coord2, coord4),
        _MetalCageVertexData(coord2, coord2, coord5),
        _MetalCageVertexData(coord2, coord4, coord1),
        _MetalCageVertexData(coord2, coord5, coord1),
        _MetalCageVertexData(coord2, coord4, coord2),
        _MetalCageVertexData(coord2, coord5, coord2),
        _MetalCageVertexData(coord4, coord1, coord1),
        _MetalCageVertexData(coord4, coord2, coord1),
        _MetalCageVertexData(coord4, coord2, coord2),
        _MetalCageVertexData(coord4, coord1, coord2),
        _MetalCageVertexData(coord5, coord1, coord1),
        _MetalCageVertexData(coord5, coord2, coord1),
        _MetalCageVertexData(coord5, coord2, coord2),
        _MetalCageVertexData(coord5, coord1, coord2),

        _MetalCageVertexData(coord1, coord3, coord4),
        _MetalCageVertexData(coord1, coord4, coord3),
        _MetalCageVertexData(coord1, coord3, coord5),
        _MetalCageVertexData(coord1, coord5, coord3),
        _MetalCageVertexData(coord2, coord3, coord4),
        _MetalCageVertexData(coord2, coord4, coord3),
        _MetalCageVertexData(coord2, coord3, coord5),
        _MetalCageVertexData(coord2, coord5, coord3),

        _MetalCageVertexData(coord3, coord1, coord4),
        _MetalCageVertexData(coord4, coord1, coord3),
        _MetalCageVertexData(coord3, coord1, coord5),
        _MetalCageVertexData(coord5, coord1, coord3),
        _MetalCageVertexData(coord3, coord2, coord4),
        _MetalCageVertexData(coord4, coord2, coord3),
        _MetalCageVertexData(coord3, coord2, coord5),
        _MetalCageVertexData(coord5, coord2, coord3),

        _MetalCageVertexData(coord3, coord4, coord1),
        _MetalCageVertexData(coord4, coord3, coord1),
        _MetalCageVertexData(coord3, coord5, coord1),
        _MetalCageVertexData(coord5, coord3, coord1),
        _MetalCageVertexData(coord3, coord4, coord2),
        _MetalCageVertexData(coord4, coord3, coord2),
        _MetalCageVertexData(coord3, coord5, coord2),
        _MetalCageVertexData(coord5, coord3, coord2),

        _MetalCageVertexData(coord1, coord6, coord6),
        _MetalCageVertexData(coord1, coord6, -coord6),
        _MetalCageVertexData(coord1, -coord6, coord6),
        _MetalCageVertexData(coord1, -coord6, -coord6),
        _MetalCageVertexData(coord2, coord6, coord6),
        _MetalCageVertexData(coord2, coord6, -coord6),
        _MetalCageVertexData(coord2, -coord6, coord6),
        _MetalCageVertexData(coord2, -coord6, -coord6),

        _MetalCageVertexData(coord6, coord1, coord6),
        _MetalCageVertexData(coord6, coord1, -coord6),
        _MetalCageVertexData(-coord6, coord1, coord6),
        _MetalCageVertexData(-coord6, coord1, -coord6),
        _MetalCageVertexData(coord6, coord2, coord6),
        _MetalCageVertexData(coord6, coord2, -coord6),
        _MetalCageVertexData(-coord6, coord2, coord6),
        _MetalCageVertexData(-coord6, coord2, -coord6),

        _MetalCageVertexData(coord6, coord6, coord1),
        _MetalCageVertexData(coord6, -coord6, coord1),
        _MetalCageVertexData(-coord6, coord6, coord1),
        _MetalCageVertexData(-coord6, -coord6, coord1),
        _MetalCageVertexData(coord6, coord6, coord2),
        _MetalCageVertexData(coord6, -coord6, coord2),
        _MetalCageVertexData(-coord6, coord6, coord2),
        _MetalCageVertexData(-coord6, -coord6, coord2),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[57],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[57]
            )
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[32],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[32]
            )
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[49],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[49]
            )
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[24],
            position=metal_bonder_vector(
                vertex_data[0],
                vertex_data[24]
            )
        ),

        EdgeData(
            vertex_data[1],
            vertex_data[26],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[26]
            )
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[56],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[56]
            )
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[48],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[48]
            )
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[34],
            position=metal_bonder_vector(
                vertex_data[1],
                vertex_data[34]
            )
        ),

        EdgeData(
            vertex_data[2],
            vertex_data[44],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[44]
            )
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[25],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[25]
            )
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[69],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[69]
            )
        ),
        EdgeData(
            vertex_data[2],
            vertex_data[51],
            position=metal_bonder_vector(
                vertex_data[2],
                vertex_data[51]
            )
        ),

        EdgeData(
            vertex_data[3],
            vertex_data[68],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[68]
            )
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[49],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[49]
            )
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[27],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[27]
            )
        ),
        EdgeData(
            vertex_data[3],
            vertex_data[46],
            position=metal_bonder_vector(
                vertex_data[3],
                vertex_data[46]
            )
        ),

        EdgeData(
            vertex_data[4],
            vertex_data[51],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[51]
            )
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[36],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[36]
            )
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[24],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[24]
            )
        ),
        EdgeData(
            vertex_data[4],
            vertex_data[61],
            position=metal_bonder_vector(
                vertex_data[4],
                vertex_data[61]
            )
        ),

        EdgeData(
            vertex_data[5],
            vertex_data[50],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[50]
            )
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[60],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[60]
            )
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[26],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[26]
            )
        ),
        EdgeData(
            vertex_data[5],
            vertex_data[38],
            position=metal_bonder_vector(
                vertex_data[5],
                vertex_data[38]
            )
        ),

        EdgeData(
            vertex_data[6],
            vertex_data[40],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[40]
            )
        ),
        EdgeData(
            vertex_data[6],
            vertex_data[25],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[25]
            )
        ),
        EdgeData(
            vertex_data[6],
            vertex_data[65],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[65]
            )
        ),
        EdgeData(
            vertex_data[6],
            vertex_data[50],
            position=metal_bonder_vector(
                vertex_data[6],
                vertex_data[50]
            )
        ),

        EdgeData(
            vertex_data[7],
            vertex_data[64],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[64]
            )
        ),
        EdgeData(
            vertex_data[7],
            vertex_data[48],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[48]
            )
        ),
        EdgeData(
            vertex_data[7],
            vertex_data[27],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[27]
            )
        ),
        EdgeData(
            vertex_data[7],
            vertex_data[42],
            position=metal_bonder_vector(
                vertex_data[7],
                vertex_data[42]
            )
        ),

        EdgeData(
            vertex_data[8],
            vertex_data[53],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[53]
            )
        ),
        EdgeData(
            vertex_data[8],
            vertex_data[32],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[32]
            )
        ),
        EdgeData(
            vertex_data[8],
            vertex_data[28],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[28]
            )
        ),
        EdgeData(
            vertex_data[8],
            vertex_data[59],
            position=metal_bonder_vector(
                vertex_data[8],
                vertex_data[59]
            )
        ),

        EdgeData(
            vertex_data[9],
            vertex_data[34],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[34]
            )
        ),
        EdgeData(
            vertex_data[9],
            vertex_data[52],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[52]
            )
        ),
        EdgeData(
            vertex_data[9],
            vertex_data[30],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[30]
            )
        ),
        EdgeData(
            vertex_data[9],
            vertex_data[58],
            position=metal_bonder_vector(
                vertex_data[9],
                vertex_data[58]
            )
        ),

        EdgeData(
            vertex_data[10],
            vertex_data[63],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[63]
            )
        ),
        EdgeData(
            vertex_data[10],
            vertex_data[28],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[28]
            )
        ),
        EdgeData(
            vertex_data[10],
            vertex_data[36],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[36]
            )
        ),
        EdgeData(
            vertex_data[10],
            vertex_data[55],
            position=metal_bonder_vector(
                vertex_data[10],
                vertex_data[55]
            )
        ),

        EdgeData(
            vertex_data[11],
            vertex_data[38],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[38]
            )
        ),
        EdgeData(
            vertex_data[11],
            vertex_data[54],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[54]
            )
        ),
        EdgeData(
            vertex_data[11],
            vertex_data[62],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[62]
            )
        ),
        EdgeData(
            vertex_data[11],
            vertex_data[30],
            position=metal_bonder_vector(
                vertex_data[11],
                vertex_data[30]
            )
        ),

        EdgeData(
            vertex_data[12],
            vertex_data[67],
            position=metal_bonder_vector(
                vertex_data[12],
                vertex_data[67]
            )
        ),
        EdgeData(
            vertex_data[12],
            vertex_data[54],
            position=metal_bonder_vector(
                vertex_data[12],
                vertex_data[54]
            )
        ),
        EdgeData(
            vertex_data[12],
            vertex_data[29],
            position=metal_bonder_vector(
                vertex_data[12],
                vertex_data[29]
            )
        ),
        EdgeData(
            vertex_data[12],
            vertex_data[40],
            position=metal_bonder_vector(
                vertex_data[12],
                vertex_data[40]
            )
        ),

        EdgeData(
            vertex_data[13],
            vertex_data[42],
            position=metal_bonder_vector(
                vertex_data[13],
                vertex_data[42]
            )
        ),
        EdgeData(
            vertex_data[13],
            vertex_data[31],
            position=metal_bonder_vector(
                vertex_data[13],
                vertex_data[31]
            )
        ),
        EdgeData(
            vertex_data[13],
            vertex_data[66],
            position=metal_bonder_vector(
                vertex_data[13],
                vertex_data[66]
            )
        ),
        EdgeData(
            vertex_data[13],
            vertex_data[52],
            position=metal_bonder_vector(
                vertex_data[13],
                vertex_data[52]
            )
        ),

        EdgeData(
            vertex_data[14],
            vertex_data[71],
            position=metal_bonder_vector(
                vertex_data[14],
                vertex_data[71]
            )
        ),
        EdgeData(
            vertex_data[14],
            vertex_data[55],
            position=metal_bonder_vector(
                vertex_data[14],
                vertex_data[55]
            )
        ),
        EdgeData(
            vertex_data[14],
            vertex_data[44],
            position=metal_bonder_vector(
                vertex_data[14],
                vertex_data[44]
            )
        ),
        EdgeData(
            vertex_data[14],
            vertex_data[29],
            position=metal_bonder_vector(
                vertex_data[14],
                vertex_data[29]
            )
        ),

        EdgeData(
            vertex_data[15],
            vertex_data[46],
            position=metal_bonder_vector(
                vertex_data[15],
                vertex_data[46]
            )
        ),
        EdgeData(
            vertex_data[15],
            vertex_data[31],
            position=metal_bonder_vector(
                vertex_data[15],
                vertex_data[31]
            )
        ),
        EdgeData(
            vertex_data[15],
            vertex_data[70],
            position=metal_bonder_vector(
                vertex_data[15],
                vertex_data[70]
            )
        ),
        EdgeData(
            vertex_data[15],
            vertex_data[53],
            position=metal_bonder_vector(
                vertex_data[15],
                vertex_data[53]
            )
        ),

        EdgeData(
            vertex_data[16],
            vertex_data[66],
            position=metal_bonder_vector(
                vertex_data[16],
                vertex_data[66]
            )
        ),
        EdgeData(
            vertex_data[16],
            vertex_data[58],
            position=metal_bonder_vector(
                vertex_data[16],
                vertex_data[58]
            )
        ),
        EdgeData(
            vertex_data[16],
            vertex_data[41],
            position=metal_bonder_vector(
                vertex_data[16],
                vertex_data[41]
            )
        ),
        EdgeData(
            vertex_data[16],
            vertex_data[33],
            position=metal_bonder_vector(
                vertex_data[16],
                vertex_data[33]
            )
        ),

        EdgeData(
            vertex_data[17],
            vertex_data[41],
            position=metal_bonder_vector(
                vertex_data[17],
                vertex_data[41]
            )
        ),
        EdgeData(
            vertex_data[17],
            vertex_data[37],
            position=metal_bonder_vector(
                vertex_data[17],
                vertex_data[37]
            )
        ),
        EdgeData(
            vertex_data[17],
            vertex_data[67],
            position=metal_bonder_vector(
                vertex_data[17],
                vertex_data[67]
            )
        ),
        EdgeData(
            vertex_data[17],
            vertex_data[62],
            position=metal_bonder_vector(
                vertex_data[17],
                vertex_data[62]
            )
        ),

        EdgeData(
            vertex_data[18],
            vertex_data[45],
            position=metal_bonder_vector(
                vertex_data[18],
                vertex_data[45]
            )
        ),
        EdgeData(
            vertex_data[18],
            vertex_data[37],
            position=metal_bonder_vector(
                vertex_data[18],
                vertex_data[37]
            )
        ),
        EdgeData(
            vertex_data[18],
            vertex_data[71],
            position=metal_bonder_vector(
                vertex_data[18],
                vertex_data[71]
            )
        ),
        EdgeData(
            vertex_data[18],
            vertex_data[63],
            position=metal_bonder_vector(
                vertex_data[18],
                vertex_data[63]
            )
        ),

        EdgeData(
            vertex_data[19],
            vertex_data[70],
            position=metal_bonder_vector(
                vertex_data[19],
                vertex_data[70]
            )
        ),
        EdgeData(
            vertex_data[19],
            vertex_data[59],
            position=metal_bonder_vector(
                vertex_data[19],
                vertex_data[59]
            )
        ),
        EdgeData(
            vertex_data[19],
            vertex_data[45],
            position=metal_bonder_vector(
                vertex_data[19],
                vertex_data[45]
            )
        ),
        EdgeData(
            vertex_data[19],
            vertex_data[33],
            position=metal_bonder_vector(
                vertex_data[19],
                vertex_data[33]
            )
        ),

        EdgeData(
            vertex_data[20],
            vertex_data[43],
            position=metal_bonder_vector(
                vertex_data[20],
                vertex_data[43]
            )
        ),
        EdgeData(
            vertex_data[20],
            vertex_data[35],
            position=metal_bonder_vector(
                vertex_data[20],
                vertex_data[35]
            )
        ),
        EdgeData(
            vertex_data[20],
            vertex_data[56],
            position=metal_bonder_vector(
                vertex_data[20],
                vertex_data[56]
            )
        ),
        EdgeData(
            vertex_data[20],
            vertex_data[64],
            position=metal_bonder_vector(
                vertex_data[20],
                vertex_data[64]
            )
        ),

        EdgeData(
            vertex_data[21],
            vertex_data[43],
            position=metal_bonder_vector(
                vertex_data[21],
                vertex_data[43]
            )
        ),
        EdgeData(
            vertex_data[21],
            vertex_data[39],
            position=metal_bonder_vector(
                vertex_data[21],
                vertex_data[39]
            )
        ),
        EdgeData(
            vertex_data[21],
            vertex_data[65],
            position=metal_bonder_vector(
                vertex_data[21],
                vertex_data[65]
            )
        ),
        EdgeData(
            vertex_data[21],
            vertex_data[60],
            position=metal_bonder_vector(
                vertex_data[21],
                vertex_data[60]
            )
        ),

        EdgeData(
            vertex_data[22],
            vertex_data[69],
            position=metal_bonder_vector(
                vertex_data[22],
                vertex_data[69]
            )
        ),
        EdgeData(
            vertex_data[22],
            vertex_data[61],
            position=metal_bonder_vector(
                vertex_data[22],
                vertex_data[61]
            )
        ),
        EdgeData(
            vertex_data[22],
            vertex_data[47],
            position=metal_bonder_vector(
                vertex_data[22],
                vertex_data[47]
            )
        ),
        EdgeData(
            vertex_data[22],
            vertex_data[39],
            position=metal_bonder_vector(
                vertex_data[22],
                vertex_data[39]
            )
        ),

        EdgeData(
            vertex_data[23],
            vertex_data[47],
            position=metal_bonder_vector(
                vertex_data[23],
                vertex_data[47]
            )
        ),
        EdgeData(
            vertex_data[23],
            vertex_data[57],
            position=metal_bonder_vector(
                vertex_data[23],
                vertex_data[57]
            )
        ),
        EdgeData(
            vertex_data[23],
            vertex_data[68],
            position=metal_bonder_vector(
                vertex_data[23],
                vertex_data[68]
            )
        ),
        EdgeData(
            vertex_data[23],
            vertex_data[35],
            position=metal_bonder_vector(
                vertex_data[23],
                vertex_data[35]
            )
        ),

    )

    num_windows = 26
    num_window_types = 2
