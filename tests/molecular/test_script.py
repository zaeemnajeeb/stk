#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.
import stk
from rdkit.Chem import AllChem as Chem


def main():
    print('testing script')
    m = Chem.MolFromSmiles('[Pd+2]')
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
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )
    print(ligand.func_groups)
    ligand.write('lig_cage.pdb')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='ligopt',
        num_cores=1,
        charge=0,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(ligand)
    m2l4_lantern = stk.metal_organic_cage.M2L4_Lantern()
    print(m2l4_lantern)
    print('--------------------------------------------------------')
    lantern = stk.ConstructedMolecule(
        building_blocks=[metal, ligand],
        topology_graph=m2l4_lantern,
        building_block_vertices={
            metal: m2l4_lantern.vertices[0:2],
            ligand: m2l4_lantern.vertices[2:]
        }
    )
    print('--------------------------------------------------------')
    print(lantern)
    lantern.write('lantern.mol')
    lantern.write('lantern.pdb')
    print('--------------------------------------------------------')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='compL',
        num_cores=4,
        charge=4,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    # Ramp up force constant.
    lantern.write('lantern_opt_1.mol')
    optimizer = stk.MetalOptimizer(
        scale=2,
        force_constant=1.0e2,
        prearrange=True,
        restrict_all_bonds=True,
        restrict_orientation=True,
        res_steps=9
    )
    optimizer.optimize(mol=lantern)
    lantern.write('lantern_opt_2.mol')
    xtb_opt.optimize(mol=lantern)
    print('--------------------------------------------------------')
    print(lantern)
    lantern.write('lantern_opt.mol')
    lantern.write('lantern_opt.pdb')

    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )
    ligand2 = stk.BuildingBlock(
        'COc1c(OC)c2ccc(-c3ccncc3)cc2c2cc(-c3ccncc3)ccc12',
        functional_groups=['pyridine_N_metal']
    )
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='lig1opt',
        num_cores=1,
        charge=0,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(ligand1)
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='lig2opt',
        num_cores=1,
        charge=0,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(ligand2)
    print(ligand1.func_groups)
    ligand1.write('lig1_cage.pdb')
    print(ligand2.func_groups)
    ligand2.write('lig2_cage.pdb')

    m2l4_lantern = stk.metal_organic_cage.M2L4_Lantern()
    print(m2l4_lantern)
    print('--------------------------------------------------------')
    hetero_lantern = stk.ConstructedMolecule(
        building_blocks=[metal, ligand1, ligand2],
        topology_graph=m2l4_lantern,
        building_block_vertices={
            metal: m2l4_lantern.vertices[0:2],
            ligand1: m2l4_lantern.vertices[2:4],
            ligand2: m2l4_lantern.vertices[4:]
        }
    )
    print('--------------------------------------------------------')
    print(hetero_lantern)
    hetero_lantern.write('hetero_lantern.mol')
    hetero_lantern.write('hetero_lantern.pdb')
    print('--------------------------------------------------------')
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='compHL',
        num_cores=4,
        charge=4,
        num_unpaired_electrons=0,
        max_runs=1,
        electronic_temperature=600,
        calculate_hessian=False,
        unlimited_memory=True
    )
    # Ramp up force constant.
    hetero_lantern.write('hetero_lantern_opt_1.mol')
    optimizer = stk.MetalOptimizer(
        scale=2,
        force_constant=1.0e2,
        prearrange=True,
        restrict_all_bonds=True,
        restrict_orientation=True,
        res_steps=9
    )
    optimizer.optimize(mol=hetero_lantern)
    hetero_lantern.write('hetero_lantern_opt_2.mol')
    xtb_opt.optimize(mol=hetero_lantern)
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='compHL',
        num_cores=4,
        charge=4,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    xtb_opt.optimize(mol=hetero_lantern)
    print('--------------------------------------------------------')
    print(hetero_lantern)
    hetero_lantern.write('hetero_lantern_opt.mol')
    hetero_lantern.write('hetero_lantern_opt.pdb')

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
    optimizer = stk.MetalOptimizer(
        scale=2,
        force_constant=1.0e2,
        prearrange=True,
        restrict_all_bonds=True,
        restrict_orientation=True,
        res_steps=9
    )
    xtb_opt = stk.XTB(
        xtb_path='/home/atarzia/software/xtb_190418/bin/xtb',
        output_dir='comp1',
        num_cores=1,
        charge=2,
        num_unpaired_electrons=0,
        max_runs=1,
        calculate_hessian=False,
        unlimited_memory=True
    )
    optimizer.optimize(mol=pdl2_sqpl_complex)
    pdl2_sqpl_complex.write('metal_complex_opt_2.mol')
    xtb_opt.optimize(mol=pdl2_sqpl_complex)
    print('--------------------------------------------------------')
    print(pdl2_sqpl_complex)
    print(pdl2_sqpl_complex.get_position_matrix()[:5])
    pdl2_sqpl_complex.write('metal_complex_opt.mol')
    pdl2_sqpl_complex.write('metal_complex_opt.pdb')

    ligand = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    # Handle multiple functional groups.
    ligand.func_groups = tuple(i for i in [ligand.func_groups[0]])
    ligand.write('lig1.pdb')

    sqpl = stk.metal_complex.SquarePlanarMonodentate(
        unsatured_vertices=[3, 4]
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
    print('--------------------------------------------------------')
    optimizer = stk.MetalOptimizer(
        scale=2,
        force_constant=1.0e2,
        prearrange=True,
        restrict_all_bonds=True,
        restrict_orientation=True,
        res_steps=6
    )
    optimizer.optimize(mol=pdl2_sqpl_complex)
    print('--------------------------------------------------------')
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

    print('testing script done')


if __name__ == "__main__":
    main()
