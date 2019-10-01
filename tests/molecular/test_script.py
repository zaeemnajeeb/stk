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
    ligand2 = stk.BuildingBlock(
        'COc1c(OC)c2ccc(-c3ccncc3)cc2c2cc(-c3ccncc3)ccc12',
        functional_groups=['pyridine_N_metal']
    )
    m2l4_lantern = stk.metal_cage.M2L4_Lantern()
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

    ligand = stk.BuildingBlock(
        'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
        functional_groups=['pyridine_N_metal']
    )
    m2l4_lantern = stk.metal_cage.M2L4_Lantern()
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

    ligand = stk.BuildingBlock(
        'NCCN',
        functional_groups=['amine_metal']
    )
    ligand.write('lig.pdb')
    print('testing script done')


if __name__ == "__main__":
    main()
