#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.
import stk
from rdkit.Chem import AllChem as Chem


def main():
    optimizer = stk.MetalOptimizer(scale=2)
    print('testing script')
    m = Chem.MolFromSmiles('[Pb+2]')
    m.AddConformer(Chem.Conformer(m.GetNumAtoms()))
    metal = stk.BuildingBlock.init_from_rdkit_mol(
        m,
        functional_groups=None,
    )
    # metal.set_position_matrix(np.array([0, 0, 0]))
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
    ligand = stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )
    ligand.write('lig.pdb')

    sqpl = stk.metal_complex.SquarePlanarBidentate()
    print(sqpl)
    print('--------------------------------------------------------')
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, ligand],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            ligand: sqpl.vertices[1:]
        }
    )
    print('--------------------------------------------------------')
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    pdl2_sqpl_complex.write('metal_complex.mol')
    pdl2_sqpl_complex.write('metal_complex.pdb')
    print('--------------------------------------------------------')
    optimizer.optimize(
        mol=pdl2_sqpl_complex,
        # xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        # output_dir='comp1',
        # num_cores=1,
        # charge=2,
        # num_unpaired_electrons=0
    )
    print('--------------------------------------------------------')
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    input()
    pdl2_sqpl_complex.write('metal_complex_opt.mol')
    pdl2_sqpl_complex.write('metal_complex_opt.pdb')

    ligand = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand.func_groups = tuple(i for i in [ligand.func_groups[0]])
    ligand.write('lig1.pdb')
    input()

    sqpl = stk.metal_complex.SquarePlanarMonodentate(
        # unsatured_vertices=[3, 4]
    )
    print(sqpl)
    print('--------------------------------------------------------')
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, ligand],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]]),
            ligand: sqpl.vertices[1:]
        }
    )
    print('--------------------------------------------------------')
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    pdl2_sqpl_complex.write('metal_complex_1.mol')
    pdl2_sqpl_complex.write('metal_complex_1.pdb')
    optimizer.optimize(
        mol=pdl2_sqpl_complex,
        # xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        # output_dir='comp2',
        # num_cores=1,
        # charge=2,
        # num_unpaired_electrons=0
    )
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    pdl2_sqpl_complex.write('metal_complex_1_opt.mol')
    pdl2_sqpl_complex.write('metal_complex_1_opt.pdb')
    import sys
    sys.exit()

    sqpl = stk.metal_complex.SquarePlanarMonodentate(
        unsatured_vertices=[1, 2, 3, 4]
    )
    print(sqpl)
    print('--------------------------------------------------------')
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal],
        topology_graph=sqpl,
        building_block_vertices={
            metal: tuple([sqpl.vertices[0]])
        }
    )
    print('--------------------------------------------------------')
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    pdl2_sqpl_complex.write('metal_complex_0.mol')
    pdl2_sqpl_complex.write('metal_complex_0.pdb')
    optimizer.optimize(
        mol=pdl2_sqpl_complex,
        # xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        # output_dir='comp3',
        # num_cores=1,
        # charge=2,
        # num_unpaired_electrons=0
    )
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    pdl2_sqpl_complex.write('metal_complex_0_opt.mol')
    pdl2_sqpl_complex.write('metal_complex_0_opt.pdb')

    print('testing script done')


if __name__ == "__main__":
    main()
