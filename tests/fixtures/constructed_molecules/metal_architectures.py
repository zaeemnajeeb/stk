import pytest
import stk
from os.path import join
from rdkit.Chem import AllChem as rdkit

data_dir = join('..', 'data')


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


@pytest.fixture(scope='function')
def tmp_bidentage_sqpl():
    metal_centre = _build_sqpl_metal_centre()
    ligand1 = stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )

    # Do construction.
    sqpl = stk.cage.SquarePlanarBidentate()
    return stk.ConstructedMolecule(
        building_blocks=[metal_centre, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal_centre: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
