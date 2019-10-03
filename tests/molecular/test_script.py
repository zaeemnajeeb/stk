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

    ligand1 = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )

    n_metals = [2, 4, 6, 12, 24]
    topologies = [
        stk.metal_cage.M2L4_Lantern(),
        stk.metal_cage.M4L8_sqpl(),
        stk.metal_cage.M6L12_cube(),
        stk.metal_cage.M12L24_sqpl(),
        stk.metal_cage.M24L48_sqpl()
    ]
    print('--------------------------------------------------------')
    for i, top in enumerate(topologies):
        print(top)
        cage = stk.ConstructedMolecule(
            building_blocks=[metal, ligand1],
            topology_graph=top,
            building_block_vertices={
                metal: top.vertices[0: n_metals[i]],
                ligand1: top.vertices[n_metals[i]:],
            }
        )
        print('-----------------------------------------------------')
        cage.write(f'cage_{i}.mol')
        cage.write(f'cage_{i}.pdb')
        cage.dump(f'cage_{i}.json')
        print('-----------------------------------------------------')


if __name__ == "__main__":
    main()
