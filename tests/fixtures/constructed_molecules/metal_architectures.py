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


@pytest.fixture(scope='function')
def tmp_bidentage_sqpl():
    metal = _build_metal()
    ligand1 = stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )

    # Do construction.
    sqpl = stk.metal_complex.SquarePlanarBidentate()
    return stk.ConstructedMolecule(
        building_blocks=[metal, ligand1],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            ligand1: sqpl.vertices[1:]
        }
    )
