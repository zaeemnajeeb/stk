import os
from os.path import join
import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit

test_dir = 'metal_cage_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _alignment(vertex, building_block, edges):
    fg_position = building_block.get_centroid(
        atom_ids=building_block.func_groups[0].get_bonder_ids()
    )
    fg_vector = stk.normalize_vector(
        fg_position - vertex.get_position()
    )

    def inner(edge_id):
        edge_vector = (
            edges[edge_id].get_position() - vertex.get_position()
        )
        return fg_vector @ stk.normalize_vector(edge_vector)

    return inner


def _test_metal_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
    assert np.allclose(
        a=bb.get_centroid(),
        b=vertex.get_position(),
        atol=1e-8,
    )


def _test_cap_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)

    assert np.allclose(
        a=bb.get_centroid(),
        b=vertex.get_position(),
        atol=1e-8,
    )
    aligned = max(
        vertex.get_edge_ids(),
        key=_alignment(vertex, bb, edges)
    )
    assert aligned is vertex._edge_ids[vertex.get_aligner_edge()]


def _test_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)

    assert np.allclose(
        a=bb.get_centroid(bb.get_bonder_ids()),
        b=vertex.get_position(),
        atol=1e-8,
    )
    aligned = max(
        vertex.get_edge_ids(),
        key=_alignment(vertex, bb, edges)
    )
    assert aligned is vertex._edge_ids[vertex.get_aligner_edge()]


def _test_assignment(vertex, bb, vertices, edges):
    assignments = vertex.assign_func_groups_to_edges(
        building_block=bb,
        vertices=vertices,
        edges=edges
    )
    assert (
        assignments[0] == vertex._edge_ids[vertex.get_aligner_edge()]
    )


def test_vertex(
    tmp_bident,
    tmp_metal
):
    topology_graphs = (
        stk.metal_cage.M2L4_Lantern(),
        stk.metal_cage.M4L8_sqpl()
    )
    building_blocks = {
        2: tmp_bident,
        4: tmp_metal
    }
    for topology_graph in topology_graphs:
        vertices = topology_graph.vertices
        edges = topology_graph.edges
        for vertex in topology_graph.vertices:
            num_edges = vertex.get_num_edges()
            bb = building_blocks[num_edges]
            if num_edges == 4:
                _test_metal_placement(vertex, bb, vertices, edges)
            else:
                _test_placement(vertex, bb, vertices, edges)
            _test_assignment(vertex, bb, vertices, edges)


def _build_metal():
    m = rdkit.MolFromSmiles('[Pd+2]')
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
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
    return metal


def test_metal_definition():
    metal = _build_metal()
    # Test metal position.
    assert np.all(np.isclose(
        metal.get_position_matrix()[0],
        np.array([[0.0, 0.0, 0.0]]),
        rtol=0
    ))
    metal.set_position_matrix(np.array([[1.0, 1.0, 1.0]]))
    assert np.all(np.isclose(
        metal.get_position_matrix()[0],
        np.array([[1.0, 1.0, 1.0]]),
        rtol=0
    ))

    # Test number of FGs and assignment.
    target_fg = stk.FunctionalGroup(
        atoms=tuple([metal.atoms[0]]),
        bonders=tuple([metal.atoms[0]]),
        deleters=(),
        fg_type=stk.FGType(
            name='metal',
            func_group_smarts='',
            bonder_smarts=[],
            deleter_smarts=([])
        ),
    )

    for fg in metal.func_groups:
        assert fg.atoms == target_fg.atoms
        assert fg.bonders == target_fg.bonders
        assert fg.deleters == target_fg.deleters

    # Test metal element.
    for atom in metal.atoms:
        assert atom.atomic_number == 46


def _test_construction(cycle, num_expected_bbs):
    name = cycle.topology_graph.__class__.__name__
    cycle.write(join(test_dir, f'{name}_{len(num_expected_bbs)}.mol'))

    for bb in cycle.get_building_blocks():
        assert cycle.building_block_counter[bb] == num_expected_bbs[bb]

    num_deleters = sum(
        len(fg.deleters)*cycle.building_block_counter[bb]
        for bb in cycle.get_building_blocks() for fg in bb.func_groups
    )
    num_bb_atoms = sum(
        len(bb.atoms)*cycle.building_block_counter[bb]
        for bb in cycle.get_building_blocks()
    )
    num_bb_bonds = sum(
        len(bb.bonds)*cycle.building_block_counter[bb]
        for bb in cycle.get_building_blocks()
    )
    # Check that the correct number of bonds got made.
    assert (
        len(cycle.construction_bonds) ==
        len(cycle.topology_graph.edges)
    )
    # Check correct total number of atoms.
    assert len(cycle.atoms) == num_bb_atoms - num_deleters
    # Check correct total number of bonds.
    assert (
        len(cycle.bonds) ==
        num_bb_bonds + len(cycle.construction_bonds) - num_deleters
    )


def test_m2l4_lantern():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    m4l8 = stk.metal_cage.M4L8_sqpl()
    lantern = stk.ConstructedMolecule(
        building_blocks=[metal, ligand1],
        topology_graph=m4l8,
        building_block_vertices={
            metal: m4l8.vertices[0:4],
            ligand1: m4l8.vertices[4:]
        }
    )
    num_expected_bbs = {
        metal: 4,
        ligand1: 8,
    }

    _test_construction(lantern, num_expected_bbs)


def test_m4l8_sqpl():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    m2l4_lantern = stk.metal_cage.M2L4_Lantern()
    lantern = stk.ConstructedMolecule(
        building_blocks=[metal, ligand1],
        topology_graph=m2l4_lantern,
        building_block_vertices={
            metal: m2l4_lantern.vertices[0:2],
            ligand1: m2l4_lantern.vertices[2:]
        }
    )
    num_expected_bbs = {
        metal: 2,
        ligand1: 4,
    }

    _test_construction(lantern, num_expected_bbs)


def test_heter_m2l4_lantern():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )
    ligand2 = stk.BuildingBlock(
        'COc1c(OC)c2ccc(-c3ccncc3)cc2c2cc(-c3ccncc3)ccc12',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
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
    num_expected_bbs = {
        metal: 2,
        ligand1: 2,
        ligand2: 2,
    }

    _test_construction(hetero_lantern, num_expected_bbs)
