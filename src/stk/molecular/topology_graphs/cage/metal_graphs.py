"""
Defines metal-cage topology graphs.

"""

import logging

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
        _MetalVertexData(0, 1, 0),
        _MetalVertexData(0, -1, 0),

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
