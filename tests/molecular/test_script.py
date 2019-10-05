#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.
import stk
from rdkit.Chem import AllChem as Chem
from collections import Counter


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

    m2 = Chem.MolFromSmiles('[Zn+2]')
    m2.AddConformer(Chem.Conformer(m2.GetNumAtoms()))
    metal2 = stk.BuildingBlock.init_from_rdkit_mol(
        m2,
        functional_groups=None,
    )
    metal2 = stk.assign_metal_fgs(
        building_block=metal2,
        coordination_info=metal_coord_info
    )

    ligand1 = stk.BuildingBlock(
        'c1cc(-c2ccncc2)cc(-c2ccncc2)c1',
        functional_groups=['pyridine_N_metal']
    )
    ligand2 = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)s2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    ligand3 = stk.BuildingBlock(
        'C(#Cc1ccc2oc3ccc(C#Cc4ccncc4)cc3c2c1)c1ccncc1',
        functional_groups=['pyridine_N_metal']
    )
    ligands = [ligand1, ligand2, ligand3]

    n_metals = [2, 4, 6, 12, 24]
    topologies = [
        stk.metal_cage.M2L4_Lantern(),
        stk.metal_cage.M4L8_sqpl(),
        stk.metal_cage.M6L12_cube(),
        stk.metal_cage.M12L24_sqpl(),
        stk.metal_cage.M24L48_sqpl()
    ]
    print('--------------------------------------------------------')
    for j, l in enumerate(ligands):
        for i, top in enumerate(topologies):
            print(top)
            cage = stk.ConstructedMolecule(
                building_blocks=[metal, l],
                topology_graph=top,
                building_block_vertices={
                    metal: top.vertices[0: n_metals[i]],
                    l: top.vertices[n_metals[i]:],
                    # metal2: top.vertices[n_metals[i]:],
                }
            )
            print('-------------------------------------------------')
            cage.write(f'cage_{j}_{i}.mol')
            cage.write(f'cage_{j}_{i}.pdb')
            cage.dump(f'cage_{j}_{i}.json')
            print('-------------------------------------------------')


if __name__ == "__main__":
    main()
