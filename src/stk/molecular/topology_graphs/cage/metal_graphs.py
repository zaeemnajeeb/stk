"""
Defines metal-cage topology graphs.

"""

import logging
import numpy as np

from .metal_base import MetalCage, _MetalVertexData
from ..topology_graph import EdgeData

logger = logging.getLogger(__name__)


class SquarePlanarMonodentate(MetalCage):
    """
    Represents a square planar metal complex topology graph.

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
        _MetalVertexData(0, 0, 0),
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(0, 0, 1),
        _MetalVertexData(0, -1, 0),
        _MetalVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[1]),
        EdgeData(vertex_data[0], vertex_data[2]),
        EdgeData(vertex_data[0], vertex_data[3]),
        EdgeData(vertex_data[0], vertex_data[4]),
    )


class SquarePlanarBidentate(MetalCage):
    """
    Represents a square planar metal complex topology graph.

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
        _MetalVertexData(0, 0, 0),
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(0, -1, 0)
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0.2, 0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[-0.2, 0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.2, -0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[-0.2, -0.2, 0]
        )
    )


class Octahedral(MetalCage):
    """
    Represents an octahedral metal complex topology graph.

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
        _MetalVertexData(0, 0, 0),
        _MetalVertexData(1, 0, 0),
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(0, 0, 1),
        _MetalVertexData(-1, 0, 0),
        _MetalVertexData(0, -1, 0),
        _MetalVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[1]),
        EdgeData(vertex_data[0], vertex_data[2]),
        EdgeData(vertex_data[0], vertex_data[3]),
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[0], vertex_data[5]),
        EdgeData(vertex_data[0], vertex_data[6]),
    )


class Octahedral_Lambda(MetalCage):
    """
    Represents an octahedral metal complex topology graph.

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
        return 5

    # Use fg_assignment to ensure the desired stereochemistry is
    # achieved.
    vertex_data = (
        _MetalVertexData(
            0, 0, 0,
            fg_assignment={0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
        ),
        _MetalVertexData(
            1, 1, 0,
            fg_assignment={0: 0, 1: 1}
        ),
        _MetalVertexData(
            0, -1, -1,
            fg_assignment={0: 4, 1: 5}
        ),
        _MetalVertexData(
            -1, 0, 1,
            fg_assignment={0: 2, 1: 3}
        ),
    )

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0.1, 0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 0.1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[0, 0, 0.1]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[-0.1, 0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0, -0.1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0, 0, -0.1]
        ),
    )


class Octahedral_Delta(MetalCage):
    """
    Represents an octahedral metal complex topology graph.

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
        return 5

    # Use fg_assignment to ensure the desired stereochemistry is
    # achieved.
    vertex_data = (
        _MetalVertexData(
            0, 0, 0,
            fg_assignment={0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
        ),
        _MetalVertexData(
            1, 1, 0,
            fg_assignment={0: 0, 1: 1}
        ),
        _MetalVertexData(
            0, -1, 1,
            fg_assignment={0: 4, 1: 2}
        ),
        _MetalVertexData(
            -1, 0, -1,
            fg_assignment={0: 5, 1: 3}
        ),
    )

    # The ordering here matters for the stereochemistry.
    # The first edge to appear between two vertices determines the
    # directionality of the binding ligand.
    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0.1, 0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 0.1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0, -0.1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0, 0, 0.1]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[0, 0, -0.1]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[-0.1, 0, 0]
        ),
    )


class Paddlewheel(MetalCage):
    """
    Represents a Paddlewheel metal complex topology graph.

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
        _MetalVertexData(0, 0.5, 0),
        _MetalVertexData(0, -0.5, 0),
        _MetalVertexData(1, 0, 0),
        _MetalVertexData(0, 0, 1),
        _MetalVertexData(-1, 0, 0),
        _MetalVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.1, 0.5, 0]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[2],
            position=[0.1, -0.5, 0]
        ),

        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[0, 0.5, 0.1]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[3],
            position=[0, -0.5, 0.1]
        ),

        EdgeData(
            vertex_data[0],
            vertex_data[4],
            position=[-0.1, 0.5, 0]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[4],
            position=[-0.1, -0.5, 0]
        ),

        EdgeData(
            vertex_data[0],
            vertex_data[5],
            position=[0, 0.5, -0.1]
        ),
        EdgeData(
            vertex_data[1],
            vertex_data[5],
            position=[0, -0.5, -0.1]
        ),
    )


class Porphyrin(MetalCage):
    """
    Represents a porphyrin metal complex topology graph.

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
        _MetalVertexData(0, 0, 0),
        _MetalVertexData(0, 0, 0)
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0.1, 0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 0.1, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[-0.1, 0, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, -0.1, 0]
        ),
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
        _MetalVertexData(0, 0.5, 0),
        _MetalVertexData(0, -0.5, 0),

        _MetalVertexData(1, 0, 0),
        _MetalVertexData(0, 0, 1),
        _MetalVertexData(-1, 0, 0),
        _MetalVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[2]),
        EdgeData(vertex_data[0], vertex_data[3]),
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[0], vertex_data[5]),
        EdgeData(vertex_data[1], vertex_data[2]),
        EdgeData(vertex_data[1], vertex_data[3]),
        EdgeData(vertex_data[1], vertex_data[4]),
        EdgeData(vertex_data[1], vertex_data[5]),
    )

    num_windows = 4
    num_window_types = 1


class M3L6_sqpl(MetalCage):
    """
    Represents a M3L6 topology graph with square planar metal.

    See :class:`.MetalCage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    R, theta = 1, 0

    vertex_data = (
        _MetalVertexData(
            R*np.cos(theta),
            R*np.sin(theta),
            0
        ),
        _MetalVertexData(
            R*np.cos(theta+(4*np.pi/3)),
            R*np.sin(theta+(4*np.pi/3)),
            0
        ),
        _MetalVertexData(
            R*np.cos(theta+(2*np.pi/3)),
            R*np.sin(theta+(2*np.pi/3)),
            0
        ),

        _MetalVertexData(
            R*np.cos((theta+np.pi/4)),
            R*np.sin((theta+np.pi/4)),
            1
        ),
        _MetalVertexData(
            R*np.cos((theta+1*np.pi/3)),
            R*np.sin((theta+1*np.pi/3)),
            -1
        ),

        _MetalVertexData(
            R*np.cos((theta+1*np.pi/3)+(4*np.pi/3)),
            R*np.sin((theta+1*np.pi/3)+(4*np.pi/3)),
            1
        ),
        _MetalVertexData(
            R*np.cos((theta+1*np.pi/3)+(4*np.pi/3)),
            R*np.sin((theta+1*np.pi/3)+(4*np.pi/3)),
            -1
        ),

        _MetalVertexData(
            R*np.cos((theta+1*np.pi/3)+(2*np.pi/3)),
            R*np.sin((theta+1*np.pi/3)+(2*np.pi/3)),
            1
        ),
        _MetalVertexData(
            R*np.cos((theta+1*np.pi/3)+(2*np.pi/3)),
            R*np.sin((theta+1*np.pi/3)+(2*np.pi/3)),
            -1
        ),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[3]),
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[0], vertex_data[5]),
        EdgeData(vertex_data[0], vertex_data[6]),

        EdgeData(vertex_data[1], vertex_data[5]),
        EdgeData(vertex_data[1], vertex_data[6]),
        EdgeData(vertex_data[1], vertex_data[7]),
        EdgeData(vertex_data[1], vertex_data[8]),

        EdgeData(vertex_data[2], vertex_data[3]),
        EdgeData(vertex_data[2], vertex_data[4]),
        EdgeData(vertex_data[2], vertex_data[7]),
        EdgeData(vertex_data[2], vertex_data[8]),
    )


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
        _MetalVertexData(1, 0, 0),
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(-1, 0, 0),
        _MetalVertexData(0, -1, 0),

        _MetalVertexData(1, 1, 1),
        _MetalVertexData(1, 1, -1),

        _MetalVertexData(1, -1, 1),
        _MetalVertexData(1, -1, -1),

        _MetalVertexData(-1, -1, 1),
        _MetalVertexData(-1, -1, -1),

        _MetalVertexData(-1, 1, 1),
        _MetalVertexData(-1, 1, -1),

    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[0], vertex_data[5]),
        EdgeData(vertex_data[0], vertex_data[6]),
        EdgeData(vertex_data[0], vertex_data[7]),

        EdgeData(vertex_data[1], vertex_data[4]),
        EdgeData(vertex_data[1], vertex_data[5]),
        EdgeData(vertex_data[1], vertex_data[10]),
        EdgeData(vertex_data[1], vertex_data[11]),

        EdgeData(vertex_data[2], vertex_data[10]),
        EdgeData(vertex_data[2], vertex_data[11]),
        EdgeData(vertex_data[2], vertex_data[8]),
        EdgeData(vertex_data[2], vertex_data[9]),

        EdgeData(vertex_data[3], vertex_data[8]),
        EdgeData(vertex_data[3], vertex_data[9]),
        EdgeData(vertex_data[3], vertex_data[6]),
        EdgeData(vertex_data[3], vertex_data[7]),
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
        _MetalVertexData(1, 0, 0),
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(-1, 0, 0),
        _MetalVertexData(0, -1, 0),
        _MetalVertexData(0, 0, 1),
        _MetalVertexData(0, 0, -1),

        _MetalVertexData(1, 1, 0),
        _MetalVertexData(1, -1, 0),
        _MetalVertexData(1, 0, 1),
        _MetalVertexData(1, 0, -1),
        _MetalVertexData(-1, 1, 0),
        _MetalVertexData(-1, -1, 0),
        _MetalVertexData(-1, 0, 1),
        _MetalVertexData(-1, 0, -1),
        _MetalVertexData(0, 1, 1),
        _MetalVertexData(0, 1, -1),
        _MetalVertexData(0, -1, 1),
        _MetalVertexData(0, -1, -1),


    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[6]),
        EdgeData(vertex_data[0], vertex_data[7]),
        EdgeData(vertex_data[0], vertex_data[8]),
        EdgeData(vertex_data[0], vertex_data[9]),

        EdgeData(vertex_data[1], vertex_data[6]),
        EdgeData(vertex_data[1], vertex_data[10]),
        EdgeData(vertex_data[1], vertex_data[14]),
        EdgeData(vertex_data[1], vertex_data[15]),

        EdgeData(vertex_data[2], vertex_data[10]),
        EdgeData(vertex_data[2], vertex_data[11]),
        EdgeData(vertex_data[2], vertex_data[12]),
        EdgeData(vertex_data[2], vertex_data[13]),

        EdgeData(vertex_data[3], vertex_data[7]),
        EdgeData(vertex_data[3], vertex_data[11]),
        EdgeData(vertex_data[3], vertex_data[16]),
        EdgeData(vertex_data[3], vertex_data[17]),

        EdgeData(vertex_data[4], vertex_data[8]),
        EdgeData(vertex_data[4], vertex_data[12]),
        EdgeData(vertex_data[4], vertex_data[14]),
        EdgeData(vertex_data[4], vertex_data[16]),

        EdgeData(vertex_data[5], vertex_data[9]),
        EdgeData(vertex_data[5], vertex_data[13]),
        EdgeData(vertex_data[5], vertex_data[15]),
        EdgeData(vertex_data[5], vertex_data[17]),
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

    vertex_data = (
        _MetalVertexData(1, 0, 0),
        _MetalVertexData(-1, 0, 0),
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(0, -1, 0),
        _MetalVertexData(0.5, 0.5, 0.707),
        _MetalVertexData(0.5, -0.5, 0.707),
        _MetalVertexData(-0.5, 0.5, 0.707),
        _MetalVertexData(-0.5, -0.5, 0.707),
        _MetalVertexData(0.5, 0.5, -0.707),
        _MetalVertexData(0.5, -0.5, -0.707),
        _MetalVertexData(-0.5, 0.5, -0.707),
        _MetalVertexData(-0.5, -0.5, -0.707),

        _MetalVertexData(1, 0.35, 0.35),
        _MetalVertexData(1, 0.35, -0.35),
        _MetalVertexData(1, -0.35, 0.35),
        _MetalVertexData(1, -0.35, -0.35),

        _MetalVertexData(-1, 0.35, 0.35),
        _MetalVertexData(-1, 0.35, -0.35),
        _MetalVertexData(-1, -0.35, 0.35),
        _MetalVertexData(-1, -0.35, -0.35),

        _MetalVertexData(0.35, 1, 0.35),
        _MetalVertexData(0.35, 1, -0.35),
        _MetalVertexData(-0.35, 1, 0.35),
        _MetalVertexData(-0.35, 1, -0.35),

        _MetalVertexData(0.35, -1, 0.35),
        _MetalVertexData(0.35, -1, -0.35),
        _MetalVertexData(-0.35, -1, 0.35),
        _MetalVertexData(-0.35, -1, -0.35),

        _MetalVertexData(0.5, 0, 0.707),
        _MetalVertexData(-0.5, 0, 0.707),
        _MetalVertexData(0, 0.5, 0.707),
        _MetalVertexData(0, -0.5, 0.707),
        _MetalVertexData(0.5, 0, -0.707),
        _MetalVertexData(-0.5, 0, -0.707),
        _MetalVertexData(0, 0.5, -0.707),
        _MetalVertexData(0, -0.5, -0.707),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[12]),
        EdgeData(vertex_data[0], vertex_data[13]),
        EdgeData(vertex_data[0], vertex_data[14]),
        EdgeData(vertex_data[0], vertex_data[15]),

        EdgeData(vertex_data[1], vertex_data[16]),
        EdgeData(vertex_data[1], vertex_data[17]),
        EdgeData(vertex_data[1], vertex_data[18]),
        EdgeData(vertex_data[1], vertex_data[19]),

        EdgeData(vertex_data[2], vertex_data[20]),
        EdgeData(vertex_data[2], vertex_data[21]),
        EdgeData(vertex_data[2], vertex_data[22]),
        EdgeData(vertex_data[2], vertex_data[23]),

        EdgeData(vertex_data[3], vertex_data[24]),
        EdgeData(vertex_data[3], vertex_data[25]),
        EdgeData(vertex_data[3], vertex_data[26]),
        EdgeData(vertex_data[3], vertex_data[27]),

        EdgeData(vertex_data[4], vertex_data[28]),
        EdgeData(vertex_data[4], vertex_data[30]),
        EdgeData(vertex_data[4], vertex_data[12]),
        EdgeData(vertex_data[4], vertex_data[20]),

        EdgeData(vertex_data[5], vertex_data[14]),
        EdgeData(vertex_data[5], vertex_data[24]),
        EdgeData(vertex_data[5], vertex_data[28]),
        EdgeData(vertex_data[5], vertex_data[31]),

        EdgeData(vertex_data[6], vertex_data[16]),
        EdgeData(vertex_data[6], vertex_data[29]),
        EdgeData(vertex_data[6], vertex_data[30]),
        EdgeData(vertex_data[6], vertex_data[22]),

        EdgeData(vertex_data[7], vertex_data[18]),
        EdgeData(vertex_data[7], vertex_data[26]),
        EdgeData(vertex_data[7], vertex_data[31]),
        EdgeData(vertex_data[7], vertex_data[29]),

        EdgeData(vertex_data[8], vertex_data[13]),
        EdgeData(vertex_data[8], vertex_data[32]),
        EdgeData(vertex_data[8], vertex_data[34]),
        EdgeData(vertex_data[8], vertex_data[21]),

        EdgeData(vertex_data[9], vertex_data[15]),
        EdgeData(vertex_data[9], vertex_data[32]),
        EdgeData(vertex_data[9], vertex_data[35]),
        EdgeData(vertex_data[9], vertex_data[25]),

        EdgeData(vertex_data[10], vertex_data[17]),
        EdgeData(vertex_data[10], vertex_data[23]),
        EdgeData(vertex_data[10], vertex_data[34]),
        EdgeData(vertex_data[10], vertex_data[33]),

        EdgeData(vertex_data[11], vertex_data[19]),
        EdgeData(vertex_data[11], vertex_data[33]),
        EdgeData(vertex_data[11], vertex_data[27]),
        EdgeData(vertex_data[11], vertex_data[35]),
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

    coord1 = 0.414
    coord2 = -0.414
    coord3 = coord1 - (coord1-coord2)/2
    coord4 = -1
    coord5 = 1
    coord6 = coord5 - (coord5-coord1)/2

    vertex_data = (
        _MetalVertexData(coord1, coord1, coord4),
        _MetalVertexData(coord1, coord1, coord5),
        _MetalVertexData(coord1, coord4, coord2),
        _MetalVertexData(coord1, coord5, coord2),
        _MetalVertexData(coord1, coord2, coord4),
        _MetalVertexData(coord1, coord2, coord5),
        _MetalVertexData(coord1, coord4, coord1),
        _MetalVertexData(coord1, coord5, coord1),
        _MetalVertexData(coord2, coord1, coord4),
        _MetalVertexData(coord2, coord1, coord5),
        _MetalVertexData(coord2, coord2, coord4),
        _MetalVertexData(coord2, coord2, coord5),
        _MetalVertexData(coord2, coord4, coord1),
        _MetalVertexData(coord2, coord5, coord1),
        _MetalVertexData(coord2, coord4, coord2),
        _MetalVertexData(coord2, coord5, coord2),
        _MetalVertexData(coord4, coord1, coord1),
        _MetalVertexData(coord4, coord2, coord1),
        _MetalVertexData(coord4, coord2, coord2),
        _MetalVertexData(coord4, coord1, coord2),
        _MetalVertexData(coord5, coord1, coord1),
        _MetalVertexData(coord5, coord2, coord1),
        _MetalVertexData(coord5, coord2, coord2),
        _MetalVertexData(coord5, coord1, coord2),

        _MetalVertexData(coord1, coord3, coord4),
        _MetalVertexData(coord1, coord4, coord3),
        _MetalVertexData(coord1, coord3, coord5),
        _MetalVertexData(coord1, coord5, coord3),
        _MetalVertexData(coord2, coord3, coord4),
        _MetalVertexData(coord2, coord4, coord3),
        _MetalVertexData(coord2, coord3, coord5),
        _MetalVertexData(coord2, coord5, coord3),

        _MetalVertexData(coord3, coord1, coord4),
        _MetalVertexData(coord4, coord1, coord3),
        _MetalVertexData(coord3, coord1, coord5),
        _MetalVertexData(coord5, coord1, coord3),
        _MetalVertexData(coord3, coord2, coord4),
        _MetalVertexData(coord4, coord2, coord3),
        _MetalVertexData(coord3, coord2, coord5),
        _MetalVertexData(coord5, coord2, coord3),

        _MetalVertexData(coord3, coord4, coord1),
        _MetalVertexData(coord4, coord3, coord1),
        _MetalVertexData(coord3, coord5, coord1),
        _MetalVertexData(coord5, coord3, coord1),
        _MetalVertexData(coord3, coord4, coord2),
        _MetalVertexData(coord4, coord3, coord2),
        _MetalVertexData(coord3, coord5, coord2),
        _MetalVertexData(coord5, coord3, coord2),

        _MetalVertexData(coord1, coord6, coord6),
        _MetalVertexData(coord1, coord6, -coord6),
        _MetalVertexData(coord1, -coord6, coord6),
        _MetalVertexData(coord1, -coord6, -coord6),
        _MetalVertexData(coord2, coord6, coord6),
        _MetalVertexData(coord2, coord6, -coord6),
        _MetalVertexData(coord2, -coord6, coord6),
        _MetalVertexData(coord2, -coord6, -coord6),

        _MetalVertexData(coord6, coord1, coord6),
        _MetalVertexData(coord6, coord1, -coord6),
        _MetalVertexData(-coord6, coord1, coord6),
        _MetalVertexData(-coord6, coord1, -coord6),
        _MetalVertexData(coord6, coord2, coord6),
        _MetalVertexData(coord6, coord2, -coord6),
        _MetalVertexData(-coord6, coord2, coord6),
        _MetalVertexData(-coord6, coord2, -coord6),

        _MetalVertexData(coord6, coord6, coord1),
        _MetalVertexData(coord6, -coord6, coord1),
        _MetalVertexData(-coord6, coord6, coord1),
        _MetalVertexData(-coord6, -coord6, coord1),
        _MetalVertexData(coord6, coord6, coord2),
        _MetalVertexData(coord6, -coord6, coord2),
        _MetalVertexData(-coord6, coord6, coord2),
        _MetalVertexData(-coord6, -coord6, coord2),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[57]),
        EdgeData(vertex_data[0], vertex_data[32]),
        EdgeData(vertex_data[0], vertex_data[49]),
        EdgeData(vertex_data[0], vertex_data[24]),

        EdgeData(vertex_data[1], vertex_data[26]),
        EdgeData(vertex_data[1], vertex_data[56]),
        EdgeData(vertex_data[1], vertex_data[48]),
        EdgeData(vertex_data[1], vertex_data[34]),

        EdgeData(vertex_data[2], vertex_data[44]),
        EdgeData(vertex_data[2], vertex_data[25]),
        EdgeData(vertex_data[2], vertex_data[69]),
        EdgeData(vertex_data[2], vertex_data[51]),

        EdgeData(vertex_data[3], vertex_data[68]),
        EdgeData(vertex_data[3], vertex_data[49]),
        EdgeData(vertex_data[3], vertex_data[27]),
        EdgeData(vertex_data[3], vertex_data[46]),

        EdgeData(vertex_data[4], vertex_data[51]),
        EdgeData(vertex_data[4], vertex_data[36]),
        EdgeData(vertex_data[4], vertex_data[24]),
        EdgeData(vertex_data[4], vertex_data[61]),

        EdgeData(vertex_data[5], vertex_data[50]),
        EdgeData(vertex_data[5], vertex_data[60]),
        EdgeData(vertex_data[5], vertex_data[26]),
        EdgeData(vertex_data[5], vertex_data[38]),

        EdgeData(vertex_data[6], vertex_data[40]),
        EdgeData(vertex_data[6], vertex_data[25]),
        EdgeData(vertex_data[6], vertex_data[65]),
        EdgeData(vertex_data[6], vertex_data[50]),

        EdgeData(vertex_data[7], vertex_data[64]),
        EdgeData(vertex_data[7], vertex_data[48]),
        EdgeData(vertex_data[7], vertex_data[27]),
        EdgeData(vertex_data[7], vertex_data[42]),

        EdgeData(vertex_data[8], vertex_data[53]),
        EdgeData(vertex_data[8], vertex_data[32]),
        EdgeData(vertex_data[8], vertex_data[28]),
        EdgeData(vertex_data[8], vertex_data[59]),

        EdgeData(vertex_data[9], vertex_data[34]),
        EdgeData(vertex_data[9], vertex_data[52]),
        EdgeData(vertex_data[9], vertex_data[30]),
        EdgeData(vertex_data[9], vertex_data[58]),

        EdgeData(vertex_data[10], vertex_data[63]),
        EdgeData(vertex_data[10], vertex_data[28]),
        EdgeData(vertex_data[10], vertex_data[36]),
        EdgeData(vertex_data[10], vertex_data[55]),

        EdgeData(vertex_data[11], vertex_data[38]),
        EdgeData(vertex_data[11], vertex_data[54]),
        EdgeData(vertex_data[11], vertex_data[62]),
        EdgeData(vertex_data[11], vertex_data[30]),

        EdgeData(vertex_data[12], vertex_data[67]),
        EdgeData(vertex_data[12], vertex_data[54]),
        EdgeData(vertex_data[12], vertex_data[29]),
        EdgeData(vertex_data[12], vertex_data[40]),

        EdgeData(vertex_data[13], vertex_data[42]),
        EdgeData(vertex_data[13], vertex_data[31]),
        EdgeData(vertex_data[13], vertex_data[66]),
        EdgeData(vertex_data[13], vertex_data[52]),

        EdgeData(vertex_data[14], vertex_data[71]),
        EdgeData(vertex_data[14], vertex_data[55]),
        EdgeData(vertex_data[14], vertex_data[44]),
        EdgeData(vertex_data[14], vertex_data[29]),

        EdgeData(vertex_data[15], vertex_data[46]),
        EdgeData(vertex_data[15], vertex_data[31]),
        EdgeData(vertex_data[15], vertex_data[70]),
        EdgeData(vertex_data[15], vertex_data[53]),

        EdgeData(vertex_data[16], vertex_data[66]),
        EdgeData(vertex_data[16], vertex_data[58]),
        EdgeData(vertex_data[16], vertex_data[41]),
        EdgeData(vertex_data[16], vertex_data[33]),

        EdgeData(vertex_data[17], vertex_data[41]),
        EdgeData(vertex_data[17], vertex_data[37]),
        EdgeData(vertex_data[17], vertex_data[67]),
        EdgeData(vertex_data[17], vertex_data[62]),

        EdgeData(vertex_data[18], vertex_data[45]),
        EdgeData(vertex_data[18], vertex_data[37]),
        EdgeData(vertex_data[18], vertex_data[71]),
        EdgeData(vertex_data[18], vertex_data[63]),

        EdgeData(vertex_data[19], vertex_data[70]),
        EdgeData(vertex_data[19], vertex_data[59]),
        EdgeData(vertex_data[19], vertex_data[45]),
        EdgeData(vertex_data[19], vertex_data[33]),

        EdgeData(vertex_data[20], vertex_data[43]),
        EdgeData(vertex_data[20], vertex_data[35]),
        EdgeData(vertex_data[20], vertex_data[56]),
        EdgeData(vertex_data[20], vertex_data[64]),

        EdgeData(vertex_data[21], vertex_data[43]),
        EdgeData(vertex_data[21], vertex_data[39]),
        EdgeData(vertex_data[21], vertex_data[65]),
        EdgeData(vertex_data[21], vertex_data[60]),

        EdgeData(vertex_data[22], vertex_data[69]),
        EdgeData(vertex_data[22], vertex_data[61]),
        EdgeData(vertex_data[22], vertex_data[47]),
        EdgeData(vertex_data[22], vertex_data[39]),

        EdgeData(vertex_data[23], vertex_data[47]),
        EdgeData(vertex_data[23], vertex_data[57]),
        EdgeData(vertex_data[23], vertex_data[68]),
        EdgeData(vertex_data[23], vertex_data[35]),

    )

    num_windows = 26
    num_window_types = 2


class M4L6_Oct(MetalCage):
    """
    Represents a M4L6 tetrahedral topology graph with octahedral metal.

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

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    vertex_data = (
        _MetalVertexData(0, 0, np.sqrt(6)/2),
        _MetalVertexData(-1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _MetalVertexData(1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _MetalVertexData(0, 2*np.sqrt(3)/3, -np.sqrt(6)/6)
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[1]),
        EdgeData(vertex_data[0], vertex_data[2]),
        EdgeData(vertex_data[0], vertex_data[3]),
        EdgeData(vertex_data[1], vertex_data[2]),
        EdgeData(vertex_data[1], vertex_data[3]),
        EdgeData(vertex_data[2], vertex_data[3]),
    )

    num_windows = 4
    num_window_types = 1


class M4L6_Oct_Spacer(MetalCage):
    """
    Represents a M4L6 tetrahedral topology graph with octahedral metal.

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

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    _v1, _v2, _v3, _v4 = _vertex_data = (
        _MetalVertexData(0, 0, np.sqrt(6)/2),
        _MetalVertexData(-1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _MetalVertexData(1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _MetalVertexData(0, 2*np.sqrt(3)/3, -np.sqrt(6)/6)
    )

    vertex_data = (
        *_vertex_data,
        _MetalVertexData.init_at_center(_v1, _v2),
        _MetalVertexData.init_at_center(_v1, _v3),
        _MetalVertexData.init_at_center(_v1, _v4),
        _MetalVertexData.init_at_center(_v2, _v3),
        _MetalVertexData.init_at_center(_v2, _v4),
        _MetalVertexData.init_at_center(_v3, _v4),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[0], vertex_data[5]),
        EdgeData(vertex_data[0], vertex_data[6]),
        EdgeData(vertex_data[1], vertex_data[4]),
        EdgeData(vertex_data[1], vertex_data[7]),
        EdgeData(vertex_data[1], vertex_data[8]),
        EdgeData(vertex_data[2], vertex_data[5]),
        EdgeData(vertex_data[2], vertex_data[7]),
        EdgeData(vertex_data[2], vertex_data[9]),
        EdgeData(vertex_data[3], vertex_data[6]),
        EdgeData(vertex_data[3], vertex_data[8]),
        EdgeData(vertex_data[3], vertex_data[9]),
    )

    num_windows = 4
    num_window_types = 1


class M4L4_Oct_Spacer(MetalCage):
    """
    Represents a M4L4 tetrahedral topology graph with octahedral metal.

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

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    _v1, _v2, _v3, _v4 = _vertex_data = (
        _MetalVertexData(0, 0, np.sqrt(6)/2),
        _MetalVertexData(-1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _MetalVertexData(1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _MetalVertexData(0, 2*np.sqrt(3)/3, -np.sqrt(6)/6)
    )

    vertex_data = (
        *_vertex_data,
        _MetalVertexData.init_at_center(_v1, _v2, _v3),
        _MetalVertexData.init_at_center(_v1, _v2, _v4),
        _MetalVertexData.init_at_center(_v1, _v3, _v4),
        _MetalVertexData.init_at_center(_v2, _v3, _v4)
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[0], vertex_data[5]),
        EdgeData(vertex_data[0], vertex_data[6]),
        EdgeData(vertex_data[1], vertex_data[4]),
        EdgeData(vertex_data[1], vertex_data[5]),
        EdgeData(vertex_data[1], vertex_data[7]),
        EdgeData(vertex_data[2], vertex_data[4]),
        EdgeData(vertex_data[2], vertex_data[6]),
        EdgeData(vertex_data[2], vertex_data[7]),
        EdgeData(vertex_data[3], vertex_data[5]),
        EdgeData(vertex_data[3], vertex_data[6]),
        EdgeData(vertex_data[3], vertex_data[7]),
    )

    num_windows = 4
    num_window_types = 1


class M8L6_Oct_Face(MetalCage):
    """
    Represents a M8L6 cube topology graph with octahedral metal.

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

    _v1, _v2, _v3, _v4, _v5, _v6, _v7, _v8 = _vertex_data = (
        _MetalVertexData(1, 1, 1),
        _MetalVertexData(1, -1, 1),
        _MetalVertexData(-1, -1, 1),
        _MetalVertexData(-1, 1, 1),
        _MetalVertexData(1, 1, -1),
        _MetalVertexData(1, -1, -1),
        _MetalVertexData(-1, -1, -1),
        _MetalVertexData(-1, 1, -1),
    )

    vertex_data = (
        *_vertex_data,
        _MetalVertexData.init_at_center(_v1, _v2, _v3, _v4),
        _MetalVertexData.init_at_center(_v1, _v2, _v5, _v6),
        _MetalVertexData.init_at_center(_v1, _v4, _v5, _v8),
        _MetalVertexData.init_at_center(_v3, _v4, _v7, _v8),
        _MetalVertexData.init_at_center(_v5, _v6, _v7, _v8),
        _MetalVertexData.init_at_center(_v2, _v3, _v6, _v7)
    )

    # Ordered with bonds of faces going clockwise.
    edge_data = (
        EdgeData(vertex_data[0], vertex_data[8]),
        EdgeData(vertex_data[1], vertex_data[8]),
        EdgeData(vertex_data[2], vertex_data[8]),
        EdgeData(vertex_data[3], vertex_data[8]),

        EdgeData(vertex_data[4], vertex_data[9]),
        EdgeData(vertex_data[5], vertex_data[9]),
        EdgeData(vertex_data[1], vertex_data[9]),
        EdgeData(vertex_data[0], vertex_data[9]),

        EdgeData(vertex_data[4], vertex_data[10]),
        EdgeData(vertex_data[0], vertex_data[10]),
        EdgeData(vertex_data[3], vertex_data[10]),
        EdgeData(vertex_data[7], vertex_data[10]),

        EdgeData(vertex_data[3], vertex_data[11]),
        EdgeData(vertex_data[2], vertex_data[11]),
        EdgeData(vertex_data[6], vertex_data[11]),
        EdgeData(vertex_data[7], vertex_data[11]),

        EdgeData(vertex_data[5], vertex_data[12]),
        EdgeData(vertex_data[4], vertex_data[12]),
        EdgeData(vertex_data[7], vertex_data[12]),
        EdgeData(vertex_data[6], vertex_data[12]),

        EdgeData(vertex_data[1], vertex_data[13]),
        EdgeData(vertex_data[5], vertex_data[13]),
        EdgeData(vertex_data[6], vertex_data[13]),
        EdgeData(vertex_data[2], vertex_data[13]),
    )


class M6L2L3_Oct(MetalCage):
    """
    Represents a M6L'2L''3 prism topology graph with octahedral metal.

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

    # Vertices of a triangle prisom so that origin is at the origin.
    # Source: http://eusebeia.dyndns.org/4d/prism3
    _v1, _v2, _v3, _v4, _v5, _v6 = _vertex_data = (
        _MetalVertexData(-1, -1/np.sqrt(3), 1),
        _MetalVertexData(1, -1/np.sqrt(3), 1),
        _MetalVertexData(0, 2/np.sqrt(3), 1),

        _MetalVertexData(-1, -1/np.sqrt(3), -1),
        _MetalVertexData(1, -1/np.sqrt(3), -1),
        _MetalVertexData(0, 2/np.sqrt(3), -1),
    )

    vertex_data = (
        *_vertex_data,
        _MetalVertexData.init_at_center(_v1, _v2, _v3),
        _MetalVertexData.init_at_center(_v4, _v5, _v6),

        _MetalVertexData.init_at_center(_v1, _v2, _v4, _v5),
        _MetalVertexData.init_at_center(_v2, _v3, _v5, _v6),
        _MetalVertexData.init_at_center(_v3, _v1, _v6, _v4),
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[6]),
        EdgeData(vertex_data[0], vertex_data[8]),
        EdgeData(vertex_data[0], vertex_data[10]),
        EdgeData(vertex_data[1], vertex_data[6]),
        EdgeData(vertex_data[1], vertex_data[8]),
        EdgeData(vertex_data[1], vertex_data[9]),
        EdgeData(vertex_data[2], vertex_data[6]),
        EdgeData(vertex_data[2], vertex_data[9]),
        EdgeData(vertex_data[2], vertex_data[10]),
        EdgeData(vertex_data[3], vertex_data[7]),
        EdgeData(vertex_data[3], vertex_data[8]),
        EdgeData(vertex_data[3], vertex_data[10]),
        EdgeData(vertex_data[4], vertex_data[7]),
        EdgeData(vertex_data[4], vertex_data[8]),
        EdgeData(vertex_data[4], vertex_data[9]),
        EdgeData(vertex_data[5], vertex_data[7]),
        EdgeData(vertex_data[5], vertex_data[9]),
        EdgeData(vertex_data[5], vertex_data[10]),
    )
