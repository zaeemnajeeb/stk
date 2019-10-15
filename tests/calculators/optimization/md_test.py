#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

import sys
import stk
import glob
from rdkit.Chem import AllChem as Chem


def add_atom_charge_flags(atom, atomkey):
    total_valence = Chem.Atom.GetTotalValence(atom)
    print('tv', total_valence)
    formal_charge = Chem.Atom.GetFormalCharge(atom)
    print('fc', formal_charge)

    atnum = int(atom.GetAtomicNum())
    print('atn', atnum)
    # Mg.
    if atnum == 12:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')
    # Al.
    elif atnum == 13:
        if total_valence != 3:
            sys.exit('charge warning')

    # Si.
    elif atnum == 14:
        if total_valence != 4:
            sys.exit('charge warning')

    # P.
    elif atnum == 15:
        if total_valence == 3:
            atomkey += '+3'
        elif total_valence == 5:
            atomkey += '+5'
        else:
            sys.exit('charge warning')

    # S.
    elif atnum == 16:
        if total_valence == 2:
            atomkey += '+2'
        elif total_valence == 4:
            atomkey += '+4'
        elif total_valence == 6:
            atomkey += '+6'
        else:
            sys.exit('charge warning')

    # Zn.
    elif atnum == 30:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')

    # Ga.
    elif atnum == 31:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # As.
    elif atnum == 33:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # Se.
    elif atnum == 34:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')

    # Cd.
    elif atnum == 48:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')

    # In.
    elif atnum == 49:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # Sb.
    elif atnum == 51:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # Te.
    elif atnum == 52:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')

    # Hg.
    elif atnum == 80:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')

    # Tl.
    elif atnum == 81:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # Pb.
    elif atnum == 82:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # Bi.
    elif atnum == 83:
        if total_valence == 3:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    # Po.
    elif atnum == 84:
        if total_valence == 2:
            atomkey += '+2'
        else:
            sys.exit('charge warning')

    # Lanthanides.
    elif atnum >= 57 and atnum <= 71:
        if total_valence == 6:
            atomkey += '+3'
        else:
            sys.exit('charge warning')

    return atomkey


def get_atom_label(atom):
    atnum = int(atom.GetAtomicNum())
    print(atnum)
    atomkey = atom.GetSymbol()
    print(atomkey)
    if len(atomkey) == 1:
        atomkey += '_'
    print(atomkey)

    table = Chem.GetPeriodicTable()

    chk1 = (Chem.PeriodicTable.GetDefaultValence(table, atnum) == -1)
    chk2 = (Chem.PeriodicTable.GetNOuterElecs(table, atnum) != 1)
    chk3 = (Chem.PeriodicTable.GetNOuterElecs(table, atnum) != 7)
    chk4 = chk2 and chk3
    print(chk1, chk2, chk3, chk4)
    if chk1 or chk4:
        print(Chem.PeriodicTable.GetDefaultValence(table, atnum))
        print(Chem.PeriodicTable.GetNOuterElecs(table, atnum))
        print('-----')
        print(Chem.Atom.GetHybridization(atom))
        hybrid = Chem.Atom.GetHybridization(atom)
        print(hybrid)
        if atnum == 84:
            print('84')
            atomkey += '3'
            print(atomkey)
            if hybrid != Chem.HybridizationType.SP3:
                print('warning 84')
        elif atnum == 80:
            print('80')
            atomkey += '1'
            print(atomkey)
            if hybrid != Chem.HybridizationType.SP:
                print('warning 80')
        else:
            if hybrid == Chem.HybridizationType.SP:
                atomkey += '1'
                print(atomkey)
            elif hybrid == Chem.HybridizationType.SP2:
                chk1a = Chem.Atom.GetIsAromatic(atom)
                bonds = Chem.Atom.GetBonds(atom)
                conjugated = False
                for bond in bonds:
                    print(bond)
                    print(Chem.Bond.GetIsConjugated(bond))
                    if Chem.Bond.GetIsConjugated(bond):
                        conjugated = True
                        break
                chk2a = conjugated
                chk3a = atnum in [6, 7, 8, 16]
                chk4a = (chk1a or chk2a)
                if chk4a and chk3a:
                    atomkey += 'R'
                    print(atomkey)
                else:
                    atomkey += '2'
                    print(atomkey)
            elif hybrid == Chem.HybridizationType.SP3:
                atomkey += '3'
                print(atomkey)
            elif hybrid == Chem.HybridizationType.SP3D:
                atomkey += '5'
                print(atomkey)
            elif hybrid == Chem.HybridizationType.SP3D2:
                atomkey += '6'
                print(atomkey)
            else:
                sys.exit('final warning. not recog')
    atomkey = add_atom_charge_flags(atom, atomkey)
    return atomkey


def has_h(bond):
    """
    Check if a bond has a H atom.

    Parameters
    ----------
    bond : :class:`stk.Bond`
        Bond to test if it has a H atom.

    Returns
    -------
    :class:`bool`
        Returns `True` if bond has H atom.
    """
    if bond.atom1.atomic_number == 1:
        return True
    if bond.atom2.atomic_number == 1:
        return True
    return False


def has_M(bond, metal_atoms):
    """
    Check if a bond has a metal atom.

    Parameters
    ----------
    bond : :class:`stk.Bond`
        Bond to test if it has a metal atom.

    metal_atoms : :class:`.list` of :class:`stk.Atom`
        List of metal atoms.

    Returns
    -------
    :class:`bool`
        Returns `True` if bond has metal atom.

    """
    if bond.atom1 in metal_atoms:
        return True
    if bond.atom2 in metal_atoms:
        return True
    return False


def write_gulp_file(mol, atom_labels, filename, metal_atoms):

    type_translator = {}
    types = set([atom_labels[i] for i in atom_labels])
    for t in types:
        not_in = True
        count = 0
        while not_in:
            count += 1
            first = t[0]
            name = f'{first}{count}'
            print(name, first, count, not_in)
            if name not in type_translator:
                type_translator[t] = name
                not_in = False

    print(types)
    print(type_translator)

    top_line = 'opti conv noautobond fix molmec cartesian\n'

    position_section = '\ncartesian\n'
    for atom in mol.atoms:
        atom_type = type_translator[atom_labels[atom.id]]
        position = mol.get_center_of_mass(atom_ids=[atom.id])
        posi_string = (
            f'{atom_type} core {round(position[0], 5)} '
            f'{round(position[1], 5)} {round(position[2], 5)}\n'
        )
        position_section += posi_string

    bond_section = '\n'
    for bond in mol.bonds:
        atom_types = [
            atom_labels[i.id]
            for i in [bond.atom1, bond.atom2]
        ]

        # Set bond orders.
        if has_h(bond):
            # H has bond order of 1.
            bond_type = ''
        elif has_M(bond, metal_atoms):
            bond_type = 'half'
        elif '_R' in atom_types[0] and '_R' in atom_types[1]:
            bond_type = 'resonant'
        elif bond.order == 1:
            bond_type = ''
        elif bond.order == 2:
            bond_type = 'double'
        elif bond.order == 3:
            bond_type = 'triple'
        print(bond, bond_type, atom_types)

        string = (
            f'connect {bond.atom1.id+1} {bond.atom2.id+1} '
            f'{bond_type}'
        )
        bond_section += string+'\n'

    species_section = '\nspecies\n'
    for spec in type_translator:
        name = type_translator[spec]
        species_section += f'{name} {spec}\n'

    library = '\nlibrary uff4mof.lib\n'

    output_section = (
        '\n'
        # 'output xyz final_\n'
        'output movie xyz steps_.xyz\n'
    )

    with open(filename, 'w') as f:
        f.write(top_line)
        f.write(position_section)
        f.write(bond_section)
        f.write(species_section)
        f.write(library)
        f.write(output_section)


def main():
    for file in glob.glob('*.json'):
        print(f'doing {file}')
        cage = stk.ConstructedMolecule.load(file)
        cage.write('structure.xyz')
        print(cage)
        gulp_opt = stk.GulpMetalOptimizer()
        metal_atoms = gulp_opt.get_metal_atoms(cage)
        metal_ids = [i.id for i in metal_atoms]
        metal_bonds, _ = gulp_opt.get_metal_bonds(cage, metal_atoms)
        edit_mol = gulp_opt.to_rdkit_mol_no_metals(
            mol=cage,
            metal_atoms=metal_atoms,
            metal_bonds=metal_bonds
        )
        print(edit_mol)

        # Get cage forcefield parameters.
        # Requires removal of H atoms.
        Chem.SanitizeMol(edit_mol)
        ff = Chem.UFFGetMoleculeForceField(edit_mol)
        atom_labels = {}
        for i in range(edit_mol.GetNumAtoms()):
            if i in metal_ids:
                print('is metal')
                print(i, cage.atoms[i])
                atom_labels[i] = 'metal'
            else:
                atom = edit_mol.GetAtomWithIdx(i)
                atom_label = get_atom_label(atom)
                print(cage.atoms[i], atom_label)
                atom_labels[i] = atom_label
            # input()

        # Write UFF4MOF specific forcefield parameters.
        # Metals.
        for atomid in atom_labels:
            if atom_labels[atomid] == 'metal':
                print(atom_labels[atomid])
                atom_labels[atomid] = 'Pd4+2'
                print(atom_labels[atomid])

        # Metal binder atoms of specific forcefields.
        # Check functional groups.

        # Write GULP file.
        write_gulp_file(
            mol=cage,
            atom_labels=atom_labels,
            filename=file.replace('.json', '.gin'),
            metal_atoms=metal_atoms
        )

        # Run.

        # Update from output.


if __name__ == "__main__":
    main()
