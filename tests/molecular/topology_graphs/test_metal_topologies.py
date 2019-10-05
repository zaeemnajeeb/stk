import os
from os.path import join
import numpy as np
import stk
from rdkit.Chem import AllChem as rdkit

test_dir = 'metal_topology_tests_output'
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
    print('>----')
    print(assignments, vertex._edge_ids[vertex.get_aligner_edge()])
    print(
        assignments[0] == vertex._edge_ids[vertex.get_aligner_edge()]
    )
    print(vertex)
    print(bb)
    print('^----')
    assert (
        assignments[0] == vertex._edge_ids[vertex.get_aligner_edge()]
    )


def test_vertex(
    tmp_monodent,
    tmp_bident,
    tmp_metal
):
    topology_graphs = (
        stk.cage.SquarePlanarMonodentate(),
        stk.cage.SquarePlanarBidentate(),
        stk.cage.M2L4_Lantern(),
        stk.cage.M4L8_sqpl(),
        stk.cage.M6L12_cube(),
        stk.cage.M12L24_sqpl(),
        stk.cage.M24L48_sqpl(),
    )
    building_blocks = {
        1: tmp_monodent,
        2: tmp_bident,
        4: tmp_metal
    }
    for topology_graph in topology_graphs:
        vertices = topology_graph.vertices
        edges = topology_graph.edges
        for vertex in topology_graph.vertices:
            num_edges = vertex.get_num_edges()
            bb = building_blocks[num_edges]
            if num_edges == 1:
                _test_cap_placement(vertex, bb, vertices, edges)
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


def _build_N_atom():
    m = rdkit.MolFromSmiles('N')
    m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
    n_atom = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=['metal_bound_N'],
    )
    return n_atom


def _build_sqpl_metal_centre():
    metal = _build_metal()
    n_atom = _build_N_atom()
    sqpl = stk.metal_centre.SquarePlanar()
    sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, n_atom],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            n_atom: sqpl.vertices[1:]
        }
    )
    return stk.BuildingBlock.init_from_molecule(
        sqpl_complex,
        functional_groups=['metal_bound_N']
    )


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
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    top = stk.cage.M2L4_Lantern()
    cage = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=top,
        building_block_vertices={
            metal_centre: top.vertices[0:2],
            ligand1: top.vertices[2:]
        }
    )
    num_expected_bbs = {
        metal_centre: 2,
        ligand1: 4,
    }

    _test_construction(cage, num_expected_bbs)


def test_m4l8_sqpl():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    top = stk.cage.M4L8_sqpl()
    cage = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=top,
        building_block_vertices={
            metal_centre: top.vertices[0:4],
            ligand1: top.vertices[4:]
        }
    )
    num_expected_bbs = {
        metal_centre: 4,
        ligand1: 8,
    }

    _test_construction(cage, num_expected_bbs)


def test_m6l12_sqpl():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    top = stk.cage.M6L12_cube()
    cage = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=top,
        building_block_vertices={
            metal_centre: top.vertices[0:6],
            ligand1: top.vertices[6:]
        }
    )
    num_expected_bbs = {
        metal_centre: 6,
        ligand1: 12,
    }

    _test_construction(cage, num_expected_bbs)


def test_m12l24_sqpl():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    top = stk.cage.M12L24_sqpl()
    cage = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=top,
        building_block_vertices={
            metal_centre: top.vertices[0:12],
            ligand1: top.vertices[12:]
        }
    )
    num_expected_bbs = {
        metal_centre: 12,
        ligand1: 24,
    }

    _test_construction(cage, num_expected_bbs)


def test_m24l48_sqpl():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    top = stk.cage.M24L48_sqpl()
    cage = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=top,
        building_block_vertices={
            metal_centre: top.vertices[0:24],
            ligand1: top.vertices[24:]
        }
    )
    num_expected_bbs = {
        metal_centre: 24,
        ligand1: 48,
    }

    _test_construction(cage, num_expected_bbs)


def test_heter_m2l4_lantern():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )
    ligand2 = stk.BuildingBlock(
        'COc1c(OC)c2ccc(-c3ccncc3)cc2c2cc(-c3ccncc3)ccc12',
        functional_groups=['pyridine_N_metal']
    )

    # Do construction.
    top = stk.cage.M2L4_Lantern()
    cage = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1, ligand2],
        topology_graph=top,
        building_block_vertices={
            metal_centre: top.vertices[0:2],
            ligand1: top.vertices[2:4],
            ligand2: top.vertices[4:]
        }
    )
    num_expected_bbs = {
        metal_centre: 2,
        ligand1: 2,
        ligand2: 2,
    }

    _test_construction(cage, num_expected_bbs)


def test_sqpl_bidentate_construction():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )

    # Do construction.
    sqpl = stk.cage.SquarePlanarBidentate()
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal_centre: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
    num_expected_bbs = {
        metal_centre: 1,
        ligand1: 2,
    }

    _test_construction(pdl2_sqpl_complex, num_expected_bbs)


def test_sqpl_monodentate_construction():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand1.func_groups = tuple(i for i in [ligand1.func_groups[0]])
    assert len(ligand1.func_groups) == 1

    # Do construction.
    sqpl = stk.cage.SquarePlanarMonodentate()
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal_centre: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
    num_expected_bbs = {
        metal_centre: 1,
        ligand1: 4,
    }

    _test_construction(pdl2_sqpl_complex, num_expected_bbs)


def test_unsat_sqpl_monodentate_construction():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand1.func_groups = tuple(i for i in [ligand1.func_groups[0]])
    assert len(ligand1.func_groups) == 1

    # Do construction.
    sqpl = stk.cage.SquarePlanarMonodentate(
        unsaturated_vertices=[3, 4]
    )
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal_centre: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
    num_expected_bbs = {
        metal_centre: 1,
        ligand1: 2,
    }

    _test_construction(pdl2_sqpl_complex, num_expected_bbs)
