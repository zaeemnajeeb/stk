"""
Defines optimizers for metallo-architectures.

"""

from itertools import combinations
import rdkit.Chem.AllChem as rdkit
import logging
import numpy as np

from ...utilities import vector_angle
from .macromodel import MacroModelMD
from ...molecular.molecules.building_block import BuildingBlock
from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class MetalOptimizer(Optimizer):
    """
    Applies optimizers that maintain metal centre coordination.

    Examples
    --------

    FULL optimize

    Just rest opt


    """

    def __init__(self, scale, force_constant, prearrange=True,
                 restrict_all_bonds=False, restrict_orientation=False,
                 res_steps=False,
                 use_cache=False):
        """
        Initialize a :class:`MetalOptimizer` instance.

        Parameters
        ----------
        scale : :class:`float`
            Distance to place ligand binder atoms from metal.

        force_constant : :class:`float`
            Force_constant to use for restricted bonds.

        prearrange : :class:`bool`, optional
            `True` to prearrange functional groups around metal centre.

        restrict_all_bonds : :class:`bool`, optional
            `True` to restrict all bonds except for ligand-FG bonds.

        restrict_orientation : :class:`bool`, optional
            `True` to restrict metal complex FG angles relative to
            topology centre of mass.

        res_steps : :class:`bool`, optional


        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """
        self._scale = scale
        self._force_constant = force_constant
        self._prearrange = prearrange
        self._restrict_all_bonds = restrict_all_bonds
        self._restrict_orientation = restrict_orientation
        self._res_steps = res_steps
        super().__init__(use_cache=use_cache)

    def prearrange_fgs(self, mol):
        """
        Prearrange the metal interacting functional groups.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """
        bbs = [i for i in mol.get_building_blocks()]
        metal_bbs = [
            i for i in bbs if i.func_groups
            if 'metal' in [j.fg_type.name for j in i.func_groups]
        ]
        bb_vertices = mol.building_block_vertices
        metal_vertices = [i for j in metal_bbs for i in bb_vertices[j]]
        metal_vertex_ids = [i.id for i in metal_vertices]
        topo_graph = mol.topology_graph

        for vertex_id in metal_vertex_ids:
            vertex = topo_graph.vertex_edge_assignments[vertex_id]
            for edge in vertex.edges:
                fg1, fg2 = edge.get_func_groups()
                fg1_name = fg1.fg_type.name
                if fg1_name == 'metal':
                    # Get new position as defined by edge position.
                    edge_vect = edge._position - vertex._position
                    edge_vect = edge_vect / np.linalg.norm(
                        edge._position - vertex._position
                    )
                    edge_vect = edge_vect * self._scale
                    new_position = edge_vect + vertex._position

                    # Get atoms to move.
                    if len(fg2.bonders) > 1:
                        raise ValueError(
                            'Metal interacting FGs should only have'
                            f'1 bonder. {fg2} does not.'
                        )
                    atom_to_move = fg2.bonders[0]
                    # Move atoms to new position.
                    pos_matrix = mol.get_position_matrix()
                    pos_matrix[atom_to_move.id] = new_position
                    mol.set_position_matrix(pos_matrix)

    def to_rdkit_mol_no_metals(self, mol, metal_atoms, metal_bonds):
        """
        Write rdkit.Mol with metals replaced by H atoms.

        """
        edit_mol = rdkit.EditableMol(rdkit.Mol())
        for atom in mol.atoms:
            if atom in metal_atoms:
                # In place of metals, add H's that will be constrained.
                # This allows the atom ids to not be changed.
                rdkit_atom = rdkit.Atom(1)
                rdkit_atom.SetFormalCharge(0)
            else:
                rdkit_atom = rdkit.Atom(atom.atomic_number)
                rdkit_atom.SetFormalCharge(atom.charge)
            edit_mol.AddAtom(rdkit_atom)

        for bond in mol.bonds:
            if bond in metal_bonds:
                # Do not add bonds to metal atoms (replaced with H's).
                continue
            edit_mol.AddBond(
                beginAtomIdx=bond.atom1.id,
                endAtomIdx=bond.atom2.id,
                order=rdkit.BondType(bond.order)
            )

        edit_mol = edit_mol.GetMol()
        rdkit_conf = rdkit.Conformer(len(mol.atoms))
        for atom_id, atom_coord in enumerate(mol._position_matrix.T):
            rdkit_conf.SetAtomPosition(atom_id, atom_coord)
            edit_mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
        edit_mol.AddConformer(rdkit_conf)

        return edit_mol

    def restricted_optimization(self, mol, metal_atoms, metal_bonds,
                                ids_to_metals,
                                rel_distance=None,
                                force_constant=None,
                                input_constraints=None):
        """
        Optimize `mol` with restrictions on metal-ligand bonds.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        rel_distance : :class:`.Molecule`
            The molecule to be optimized.

        force_constant : :class:`.Molecule`
            The molecule to be optimized.

        input_constraints : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Write rdkit molecule with metal atoms and bonds deleted.
        edit_mol = self.to_rdkit_mol_no_metals(
            mol,
            metal_atoms=metal_atoms,
            metal_bonds=metal_bonds
        )

        # Constrain selected bonds, angles and dihedrals and metal
        # atoms.
        rdkit.SanitizeMol(edit_mol)
        ff = rdkit.UFFGetMoleculeForceField(edit_mol)

        # Constrain the metal centre.
        self.apply_metal_centre_constraints(
            mol,
            ff,
            metal_bonds,
            metal_atoms
        )

        # Add angular constraints that enforce relative orientation
        # between metal complexes + the topology centre of mass.
        if self._restrict_orientation:
            self.apply_orientation_restriction(
                ff,
                mol,
                metal_bonds,
                metal_atoms
            )

        # Constrain all bonds and angles based on input structure
        # except for:
        # (1) bonds including metals
        # (2) bonds including atoms bonded to metals
        if self._restrict_all_bonds:
            for const in input_constraints:
                constraint = input_constraints[const]
                if constraint['type'] == 'bond':
                    # Add distance constraints in place of metal bonds.
                    ff.UFFAddDistanceConstraint(
                        idx1=constraint['idx1'],
                        idx2=constraint['idx2'],
                        relative=True,
                        minLen=1.0,
                        maxLen=1.0,
                        forceConstant=constraint['fc']
                    )
                elif constraint['type'] == 'angle':
                    ff.UFFAddAngleConstraint(
                        idx1=constraint['idx1'],
                        idx2=constraint['idx2'],
                        idx3=constraint['idx3'],
                        relative=False,
                        minAngleDeg=constraint['angle'],
                        maxAngleDeg=constraint['angle'],
                        forceConstant=constraint['fc']
                    )
                elif constraint['type'] == 'torsion':
                    ff.UFFAddTorsionConstraint(
                        idx1=constraint['idx1'],
                        idx2=constraint['idx2'],
                        idx3=constraint['idx3'],
                        idx4=constraint['idx4'],
                        relative=False,
                        minDihedralDeg=constraint['torsion'],
                        maxDihedralDeg=constraint['torsion'],
                        forceConstant=constraint['fc']
                    )

        # For bonds (2), a weak force constant is applied to minimize
        # to rel_distance.
        if rel_distance is not None and force_constant is not None:
            for bond in mol.bonds:
                idx1 = bond.atom1.id
                idx2 = bond.atom2.id
                if idx1 in ids_to_metals or idx2 in ids_to_metals:
                    # Do not restrict H atom distances.
                    if self.has_h(bond):
                        continue
                    if self.has_M(bond, metal_atoms):
                        continue
                    # Add distance constraints in place of metal bonds.
                    ff.UFFAddDistanceConstraint(
                        idx1=idx1,
                        idx2=idx2,
                        relative=True,
                        minLen=rel_distance,
                        maxLen=rel_distance,
                        forceConstant=force_constant
                    )

        # Perform UFF optimization with rdkit.
        ff.Minimize(maxIts=50)

        # Update stk molecule from optimized molecule. This should
        # only modify atom positions, which means metal atoms will be
        # reinstated.
        mol.update_from_rdkit_mol(edit_mol)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Find all metal atoms and atoms they are bonded to.
        metal_a_no = list(range(21, 31))
        metal_a_no += list(range(39, 49))+list(range(72, 81))
        metal_atoms = []
        for atom in mol.atoms:
            if atom.atomic_number in metal_a_no:
                metal_atoms.append(atom)

        metal_bonds = []
        ids_to_metals = []
        for bond in mol.bonds:
            if bond.atom1 in metal_atoms:
                metal_bonds.append(bond)
                ids_to_metals.append(bond.atom2.id)
            elif bond.atom2 in metal_atoms:
                metal_bonds.append(bond)
                ids_to_metals.append(bond.atom1.id)

        # Second step is to perform a forcefield optimisation that
        # only optimises non metal atoms that are not bonded to the
        # metal.
        # This is done in loops, where the metal-ligand bonds are
        # slowly relaxed to normal values.
        # The rest of the ligand is constrained to the input value.
        input_constraints = self.get_input_constraints(
            mol,
            ids_to_metals,
            metal_atoms,
            include_bonders=False
        )

        # First step is to pre-arrange the metal centre based on the
        # MetalComplex topology.
        if self._prearrange:
            self.prearrange_fgs(mol)

    def get_input_constraints(self, mol, ids_to_metals, metal_atoms,
                              include_bonders=False):
        """
        Get a series of constraint definitions based on mol.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        constraints : :class:`dict`
            Dictionary of all constraints of the form:
                bonds: {KEY: {VALUES}}
                angles: {KEY: {VALUES}}
                torsions: {KEY: {VALUES}}

        """
        constraints = {}
        # Bond constraints.
        for bond in mol.bonds:
            idx1 = bond.atom1.id
            idx2 = bond.atom2.id
            if idx1 in ids_to_metals or idx2 in ids_to_metals:
                continue
            # Do not restrict H atom distances.
            if self.has_h(bond):
                continue
            if self.has_M(bond, metal_atoms):
                continue
            constraints[bond] = {
                'idx1': idx1,
                'idx2': idx2,
                'type': 'bond',
                'fc': 1.0e2
            }

        # Angle constraints
        # Add angle constraints to angles containing H atoms.
        for bonds in combinations(mol.bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.atom1, bond1.atom2]
            bond2_atoms = [bond2.atom1, bond2.atom2]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) > 3:
                continue
            # Only want angles with at least 1 X-H bond.
            # if not self.has_h(bond1) and not self.has_h(bond2):
            #     continue
            # Do not want any metal containing bonds.
            if any([self.has_M(i, metal_atoms) for i in bonds]):
                continue
            for atom in pres_atoms:
                if atom in bond1_atoms and atom in bond2_atoms:
                    idx2 = atom.id
                elif atom in bond1_atoms:
                    idx1 = atom.id
                elif atom in bond2_atoms:
                    idx3 = atom.id
            if not include_bonders:
                # Do not restrict angles of bonds to atoms bonded to
                # metals.
                bonded_to_metals = [
                    idx1 in ids_to_metals,
                    idx2 in ids_to_metals,
                    idx3 in ids_to_metals
                ]
                if any(bonded_to_metals):
                    continue
            pos1 = [
                i for i in mol.get_atom_coords(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_coords(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_coords(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            if self.has_h(bond1) or self.has_h(bond2):
                # Increased force constant for angles with H atoms.
                constraints[(bond1, bond2)] = {
                    'idx1': idx1,
                    'idx2': idx2,
                    'idx3': idx3,
                    'type': 'angle',
                    'angle': np.degrees(angle),
                    'fc': 1.0e2
                }
            else:
                constraints[(bond1, bond2)] = {
                    'idx1': idx1,
                    'idx2': idx2,
                    'idx3': idx3,
                    'type': 'angle',
                    'angle': np.degrees(angle),
                    'fc': 1.0e1
                }

            # # Dihedral constraints.
            # for bond3 in mol.bonds:
            #     if bond3 == bond2 or bond3 == bond1:
            #         continue
            #     # Want bonds connected to idx3 atom.
            #     bond3_atoms = [bond3.atom1, bond3.atom2]
            #     if bond3.atom1.id == idx3:
            #         idx4 = bond3.atom2.id
            #     elif bond3.atom2.id == idx3:
            #         idx4 = bond3.atom1.id
            #     else:
            #         continue
            #     pres_atoms = list(set(
            #         bond1_atoms + bond2_atoms + bond3_atoms
            #     ))
            #     # If there are more than 3 atoms, implies at least two
            #     # independant bonds.
            #     if len(pres_atoms) > 4:
            #         continue
            #     print('p', bond1, bond2, bond3)
            #     # Do not want any metal containing bonds.
            #     if self.has_M(bond3, metal_atoms):
            #         continue
            #
            #     # Do not restrict angles of bonds to atoms bonded to
            #     # metals.
            #     bonded_to_metals = [
            #         idx1 in ids_to_metals,
            #         idx2 in ids_to_metals,
            #         idx3 in ids_to_metals,
            #         idx4 in ids_to_metals
            #     ]
            #     if any(bonded_to_metals):
            #         continue
            #     print(idx1, idx2, idx3, idx4)
            #     pos1 = [
            #         i for i in mol.get_atom_coords(atom_ids=[idx1])
            #     ][0]
            #     pos2 = [
            #         i for i in mol.get_atom_coords(atom_ids=[idx2])
            #     ][0]
            #     pos3 = [
            #         i for i in mol.get_atom_coords(atom_ids=[idx3])
            #     ][0]
            #     pos4 = [
            #         i for i in mol.get_atom_coords(atom_ids=[idx4])
            #     ][0]
            #     torsion = get_dihedral(pos1, pos2, pos3, pos4)
            #     print(bond1, bond2, bond3, torsion)
            #     constraints[(bond1, bond2)] = {
            #         'idx1': idx1,
            #         'idx2': idx2,
            #         'idx3': idx3,
            #         'idx4': idx4,
            #         'type': 'torsion',
            #         'torsion': torsion,
            #         'fc': 0
            #     }
        return constraints

    def has_h(self, bond):
        """
        Check if a bond has a H atom.

        """
        if bond.atom1.atomic_number == 1:
            return True
        if bond.atom2.atomic_number == 1:
            return True
        return False

    def has_M(self, bond, metal_atoms):
        """
        Check if a bond has a metal atom.

        """
        if bond.atom1 in metal_atoms:
            return True
        if bond.atom2 in metal_atoms:
            return True
        return False

    def apply_metal_centre_constraints(self, mol, ff, metal_bonds,
                                       metal_atoms):
        """
        Applies UFF metal centre constraints.

        """
        metal_dist = max(
            bb.get_maximum_diameter()
            for bb in mol.building_block_vertices
        )
        for atoms in combinations(metal_atoms, r=2):
            ff.UFFAddDistanceConstraint(
                idx1=atoms[0].id,
                idx2=atoms[1].id,
                relative=False,
                minLen=metal_dist,
                maxLen=metal_dist,
                forceConstant=0.25e1
            )

        # Add constraints to UFF to hold metal geometry in place.
        for bond in metal_bonds:
            idx1 = bond.atom1.id
            idx2 = bond.atom2.id
            # Add distance constraints in place of metal bonds.
            ff.UFFAddDistanceConstraint(
                idx1=idx1,
                idx2=idx2,
                relative=False,
                minLen=self._scale,
                maxLen=self._scale,
                forceConstant=self._force_constant
            )

        # Also implement angular constraints.
        for bonds in combinations(metal_bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.atom1, bond1.atom2]
            bond2_atoms = [bond2.atom1, bond2.atom2]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) > 3:
                continue
            for atom in pres_atoms:
                if atom in bond1_atoms and atom in bond2_atoms:
                    idx2 = atom.id
                elif atom in bond1_atoms:
                    idx1 = atom.id
                elif atom in bond2_atoms:
                    idx3 = atom.id
            pos1 = [
                i for i in mol.get_atom_coords(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_coords(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_coords(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=idx3,
                relative=False,
                minAngleDeg=np.degrees(angle)-2,
                maxAngleDeg=np.degrees(angle)+2,
                forceConstant=1.0e4
            )

    def apply_orientation_restriction(self, ff, mol, metal_bonds,
                                      metal_atoms):
        """
        Applies UFF relative orientation restrcitions.

        """
        # Add a fixed point.
        COM = mol.get_center_of_mass()
        nidx = ff.AddExtraPoint(
            x=COM[0],
            y=COM[1],
            z=COM[2],
            fixed=True
        )
        ff.Initialize()
        for bond in metal_bonds:
            if bond.atom1 in metal_atoms:
                idx2 = bond.atom1.id
                idx1 = bond.atom2.id
            elif bond.atom2 in metal_atoms:
                idx1 = bond.atom1.id
                idx2 = bond.atom2.id
            pos1 = [
                i for i in mol.get_atom_coords(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_coords(atom_ids=[idx2])
            ][0]
            pos3 = COM
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=nidx-1,
                relative=False,
                minAngleDeg=np.degrees(angle)-2,
                maxAngleDeg=np.degrees(angle)+2,
                forceConstant=1.0e4
            )
            # Add fixed points at all other vertex positions (this
            # can be achieved using the COM of all building_blocks)
            # and enforce the relative orientation to all of them
            # also.
            # bb_atoms_ids = {}
            # for atom in mol.atoms:
            #     # Skip H atoms.
            #     if atom.atomic_number == 1:
            #         continue
            #     atom_bb_fgs = list(set([
            #         i.fg_type.name
            #         for i in atom.building_block.func_groups
            #     ]))
            #     print(atom, atom_bb_fgs)
            #     # Skip metal building blocks.
            #     if 'metal' in atom_bb_fgs:
            #         continue
            #     bb_key = (
            #         atom.building_block,
            #         atom.building_block_id
            #     )
            #     print(bb_key)
            #     if bb_key not in bb_atoms_ids:
            #         bb_atoms_ids[bb_key] = []
            #     bb_atoms_ids[bb_key].append(atom.id)
            #
            # print(bb_atoms_ids)
            # input()
            # for bb_key in bb_atoms_ids:
            #     # Add a fixed point.
            #     COM = mol.get_center_of_mass(
            #         atom_ids=bb_atoms_ids[bb_key]
            #     )
            #     print(COM)
            #     input(bb_key)
            #     nidx = ff.AddExtraPoint(
            #         x=COM[0],
            #         y=COM[1],
            #         z=COM[2],
            #         fixed=True
            #     )
            #     ff.Initialize()
            #     for bond in metal_bonds:
            #         if bond.atom1 in metal_atoms:
            #             idx2 = bond.atom1.id
            #             idx1 = bond.atom2.id
            #         elif bond.atom2 in metal_atoms:
            #             idx1 = bond.atom1.id
            #             idx2 = bond.atom2.id
            #         pos1 = [
            #             i for i in mol.get_atom_coords(
            #                 atom_ids=[idx1]
            #             )
            #         ][0]
            #         pos2 = [
            #             i for i in mol.get_atom_coords(
            #                 atom_ids=[idx2]
            #             )
            #         ][0]
            #         pos3 = COM
            #         v1 = pos1 - pos2
            #         v2 = pos3 - pos2
            #         angle = vector_angle(v1, v2)
            #         print(bond)
            #         print(v1, v2, angle)
            #         ff.UFFAddAngleConstraint(
            #             idx1=idx1,
            #             idx2=idx2,
            #             idx3=nidx-1,
            #             relative=False,
            #             minAngleDeg=np.degrees(angle)-2,
            #             maxAngleDeg=np.degrees(angle)+2,
            #             forceConstant=1.0e6
            #         )
