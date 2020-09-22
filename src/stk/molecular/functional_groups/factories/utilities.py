import rdkit.Chem.AllChem as rdkit


def _get_atom_ids(query, molecule):
    """
    Yield the ids of atoms in `molecule` which match `query`.

    Multiple substructures in `molecule` can match `query` and
    therefore each set is yielded as a group.

    Parameters
    ----------
    query : :class:`str`
        A SMARTS string used to query atoms.

    molecule : :class:`.Molecule`
        A molecule whose atoms should be queried.

    Yields
    ------
    :class:`tuple` of :class:`int`
        The ids of atoms in `molecule` which match `query`.

    """

    rdkit_mol = molecule.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    yield from rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(query),
    )

def _find_R_ids(molecule, deleters, bonder, previous_deleters = set()):
    if len(previous_deleters) is len(deleters):
        return deleters
    else:
        new=set()
        for atoms in deleters.difference(previous_deleters):
            for bonds in molecule.get_bonds():
                if atoms is bonds.get_atom1().get_id():
                    new.add(bonds.get_atom2().get_id())
                elif atoms is bonds.get_atom2().get_id():
                    new.add(bonds.get_atom1().get_id())
        new.discard(bonder)  
        previous_deleters=deleters.copy()
        deleters=deleters.union(new)
        return _find_R_ids(
            molecule=molecule,
            deleters=deleters,
            previous_deleters=previous_deleters,
            bonder=bonder
            )