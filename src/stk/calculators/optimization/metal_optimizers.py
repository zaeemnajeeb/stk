"""
Defines optimizers for metallo-architectures.

"""

from itertools import combinations
import rdkit.Chem.AllChem as rdkit
from os.path import join
import logging
import numpy as np
import subprocess as sp
import re
import os
import uuid

from ...utilities import vector_angle
from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class MetalOptimizer(Optimizer):
    """
    Applies forcefield optimizers that can handle metal centres.

    Notes
    -----
    By default, :meth:`optimize` will run a restricted optimization
    using constraints and the UFF. To implement this, metal atoms are
    replaced by noninteracting H atoms, and constraints are applied
    to maintain the metal centre geometry.

    Restrictions are applied to the ligand with respect to its input
    structure. So if that is poorly optimised, then the output will be
    also.

    Attributes
    ----------

    _force_constant : :class:`float`
        Force constant to use for restricted metal-ligand bonds.
        All other force constants are hard-coded/standard values.

    _rel_distance : :class:`float`
        Set the relative distance to optimise the metal-ligand
        bonds to in each optimisation step.

    _res_steps : :class:`int`
        Number of optimisation steps to run. Each optimisation step
        is very short to ensure the structure does not over-react
        to the extreme forces on it.

    _restrict_all_bonds : :class:`bool`
        `True` to restrict all bonds except for ligand-FG bonds.

    _restrict_orientation : :class:`bool`
        `True` to restrict metal complex FG angles relative to
        topology centre of mass.

    _metal_a_no : :class:`list` of :class:`int`
        Atomic numbers of metals based on the periodic table.

    Examples
    --------

    :class:`MetalOptimizer` allows for the restricted optimization of
    :class:`ConstructedMolecule` containing metals.

    .. code-block:: python

        import stk
        from rdkit.Chem import AllChem as rdkit
        import numpy as np

        # Define an stk BuildingBlock with no functional groups and a
        # single metal (Pd with 2+ charge) atom.
        m = rdkit.MolFromSmiles('[Pd+2]')
        m.AddConformer(rdkit.Conformer(m.GetNumAtoms()))
        metal = stk.BuildingBlock.init_from_rdkit_mol(
            m,
            functional_groups=None,
        )

        # Manually set functional group information for the metal
        # centre. The key of this dictionary is the fg.id.
        # Pd2+ is 4 coordinate, so we write 4 metal functiongal groups.
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
            'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
            functional_groups=['pyridine_N_metal']
        )
        # Handle multiple functional groups on the ligand to ensure
        # only one functional group is bound to the metal.
        ligand.func_groups = tuple(i for i in [ligand.func_groups[0]])

        # Construct four-coordinated Pd 2+ square planar complex.
        sqpl = stk.metal_complex.SquarePlanarMonodentate()
        pdl2_sqpl_complex = stk.ConstructedMolecule(
            building_blocks=[metal, ligand],
            topology_graph=sqpl,
            # Assign the metal to vertex 0, although this is not a
            # requirement.
            building_block_vertices={
                metal: tuple([sqpl.vertices[0]]),
                ligand: sqpl.vertices[1:]
            }
        )

        # Define an optimizer, with ligand bonder atoms placed 1.9
        # Angstrom from the metal centres. A weak force_constant is
        # used here to not overshadow all other forces. A slow
        # optimisation of the ligand - bonder bonds is done over
        # res_steps, where each step is a very short RDKit UFF
        # optimisation with constraints.
        # All ligand bonds and angles have restricitons applied based
        # on the input structure.
        optimizer = stk.MetalOptimizer(
            scale=1.9,
            force_constant=1.0e2,
            rel_distance=0.9,
            res_steps=50,
            restrict_all_bonds=True,
            restrict_orientation=True
        )

        # Optimize.
        optimizer.optimize(mol=pdl2_sqpl_complex)

    The optimisation also works for metal cages, although for some
    structures you may want to modify the number of steps or scale.

    .. code-block:: python
        m2l4_lantern = stk.metal_cage.M2L4_Lantern()

        ligand = stk.BuildingBlock(
            'C(#Cc1cccc(C#Cc2cccnc2)c1)c1cccnc1',
            functional_groups=['pyridine_N_metal']
        )

        lantern = stk.ConstructedMolecule(
            building_blocks=[metal, ligand],
            topology_graph=m2l4_lantern,
            building_block_vertices={
                metal: m2l4_lantern.vertices[0:2],
                ligand: m2l4_lantern.vertices[2:]
            }
        )

        optimizer.optimize(mol=lantern)

    """

    def __init__(
        self,
        metal_binder_distance,
        metal_binder_fc,
        binder_ligand_fc,
        ignore_vdw,
        rel_distance,
        res_steps,
        max_iterations,
        do_long_opt,
        restrict_bonds=False,
        restrict_angles=False,
        restrict_orientation=False,
        use_cache=False
    ):
        """
        Initialize a :class:`MetalOptimizer` instance.

        Parameters
        ----------

        force_constant : :class:`float`
            Force constant to use for restricted metal-ligand bonds.
            All other force constants are hard-coded/standard values.

        rel_distance : :class:`float`
            Set the relative distance to optimise the metal-ligand
            bonds to.

        res_steps : :class:`int`
            Number of optimisation steps to run. Each optimisation step
            is very short to ensure the structure does not over-react
            to the extreme forces on it.

        restrict_all_bonds : :class:`bool`, optional
            `True` to restrict all bonds except for ligand-FG bonds.

        restrict_orientation : :class:`bool`, optional
            `True` to restrict metal complex FG angles relative to
            topology centre of mass.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """
        self._metal_binder_distance = metal_binder_distance
        self._metal_binder_fc = metal_binder_fc
        self._binder_ligand_fc = binder_ligand_fc
        self._ignore_vdw = ignore_vdw
        self._rel_distance = rel_distance
        self._restrict_bonds = restrict_bonds
        self._restrict_angles = restrict_angles
        self._restrict_orientation = restrict_orientation
        self._res_steps = res_steps
        self._max_iterations = max_iterations
        self._do_long_opt = do_long_opt

        self.metal_a_no = list(range(21, 31))
        self.metal_a_no += list(range(39, 49))+list(range(72, 81))

        super(MetalOptimizer, self).__init__(use_cache=use_cache)

    def to_rdkit_mol_no_metals(self, mol, metal_atoms, metal_bonds):
        """
        Write :class:`rdkit.Mol` with metals replaced by H atoms.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        Returns
        -------
        edit_mol : :class:`rdkit.Mol`
            RDKit molecule with metal atoms replaced with H atoms.

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

    def get_input_constraints(
        self,
        mol,
        ids_to_metals,
        metal_atoms,
        include_bonders=False
    ):
        """
        Get a series of constraint definitions based on mol.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        ids_to_metals : :class:`.list` of :class:`int`
            List of atom ids of atoms bonded to metals.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        include_bonders : :class:`bool`, optional
            Whether to constrain the angles of atoms bonded to the
            metal atoms. Defaults to `False`

        Returns
        -------
        constraints : :class:`dict`
            Dictionary of all constraints of the form:
                bonds: {stk.bond: {
                    'idx1': :class:`int`,
                    'idx2': :class:`int`,
                    'type': 'bond',
                    'fc': :class:`float`
                }}
                angles: {(stk.bond2, stk.bond): {
                    'idx1': :class:`int`,
                    'idx2': :class:`int`,
                    'idx3': :class:`int`,
                    'type': 'angle',
                    'angle': :class:`float`,
                    'fc': :class:`float`
                }}
                torsions: {(stk.bond, stk.bond, stk.bond): {
                    'idx1': :class:`int`,
                    'idx2': :class:`int`,
                    'idx3': :class:`int`,
                    'idx4': :class:`int`,
                    'type': 'torsion',
                    'torsion': :class:`float`,
                    'fc': :class:`float`
                }}

        """
        constraints = {}
        # Bond constraints. Set to have a force constant of 1E2.
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
            length = mol.get_atom_distance(idx1, idx2)
            constraints[bond] = {
                'idx1': idx1,
                'idx2': idx2,
                'type': 'bond',
                'length': length,
                'fc': 1.0e2
            }

        # Angle constraints
        # Add angle constraints to angles containing H atoms.
        # Set to have a force constant of 1E1, unless the angle
        # includes H atoms, then it has a force constant of 1E2.
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
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_positions(atom_ids=[idx3])
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

        return constraints

    def _restricted_optimization(
        self,
        ff,
        mol,
        edit_mol,
        metal_atoms,
        metal_bonds,
        ids_to_metals,
        input_constraints=None
    ):
        """
        Optimize `mol` with restrictions on metal-ligand bonds.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        ids_to_metals : :class:`.list` of :class:`int`
            List of atom ids of atoms bonded to metals.

        input_constraints : :class:`dict`
            Dictionary of constraints to apply.

        Returns
        -------
        None : :class:`NoneType`

        """

        # For bonds between ligand bonders and the rest of the ligand,
        # a weak force constant is applied to minimize to rel_distance.
        # This is the slow relaxation of the high-force bonds.
        if self._rel_distance is not None:
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
                    distance = mol.get_atom_distance(idx1, idx2)
                    ff.UFFAddDistanceConstraint(
                        idx1=idx1,
                        idx2=idx2,
                        relative=False,
                        minLen=self._rel_distance*distance,
                        maxLen=self._rel_distance*distance,
                        forceConstant=self._binder_ligand_fc
                    )

        # Perform UFF optimization with rdkit.
        ff.Minimize(maxIts=self._max_iterations)

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
        metal_atoms = self.get_metal_atoms(mol)
        metal_bonds, ids_to_metals = self.get_metal_bonds(
            mol,
            metal_atoms
        )

        input_constraints = self.get_input_constraints(
            mol,
            ids_to_metals,
            metal_atoms,
            include_bonders=False
        )

        # Perform a forcefield optimisation that
        # only optimises non metal atoms that are not bonded to the
        # metal.

        # Write rdkit molecule with metal atoms and bonds deleted.
        edit_mol = self.to_rdkit_mol_no_metals(
            mol,
            metal_atoms=metal_atoms,
            metal_bonds=metal_bonds
        )

        # Non-bonded interaction interactions need to be turned on
        # because at the point of intialisation, the molecule are
        # technically separate (not bonded) and are treated as
        # fragments.
        rdkit.SanitizeMol(edit_mol)
        ff = rdkit.UFFGetMoleculeForceField(
            edit_mol,
            ignoreInterfragInteractions=self._ignore_vdw
        )

        # Constrain selected bonds, angles and dihedrals and metal
        # atoms.
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
        for constraint in input_constraints:
            const = input_constraints[constraint]
            if const['type'] == 'bond' and self._restrict_bonds:
                # Add distance constraints in place of metal bonds.
                ff.UFFAddDistanceConstraint(
                    idx1=const['idx1'],
                    idx2=const['idx2'],
                    relative=False,
                    minLen=const['length'],
                    maxLen=const['length'],
                    forceConstant=const['fc']
                )
            elif const['type'] == 'angle' and self._restrict_angles:
                ff.UFFAddAngleConstraint(
                    idx1=const['idx1'],
                    idx2=const['idx2'],
                    idx3=const['idx3'],
                    relative=False,
                    minAngleDeg=const['angle'],
                    maxAngleDeg=const['angle'],
                    forceConstant=const['fc']
                )

        # Optimisation with UFF in RDKit. This method uses constraints
        # on the metal centre to attempt to enforce the metal geometry
        # described by the metal topology.
        # For RDKit, the optimisation is done in loops, where the
        # metal-ligand and adjacent bonds are slowly relaxed to normal
        # values.
        # The rest of the ligands are constrained to the input value.
        for i in range(self._res_steps):
            self._restricted_optimization(
                mol=mol,
                edit_mol=edit_mol,
                ff=ff,
                metal_atoms=metal_atoms,
                metal_bonds=metal_bonds,
                ids_to_metals=ids_to_metals,
                input_constraints=input_constraints
            )

        # Finish with one long optimisation.
        if self._do_long_opt:
            self._max_iterations = 500
            self._restricted_optimization(
                mol=mol,
                edit_mol=edit_mol,
                ff=ff,
                metal_atoms=metal_atoms,
                metal_bonds=metal_bonds,
                ids_to_metals=ids_to_metals,
                input_constraints=input_constraints
            )

    def apply_metal_centre_constraints(
        self,
        mol,
        ff,
        metal_bonds,
        metal_atoms
    ):
        """
        Applies UFF metal centre constraints.

        Parameters
        ----------
        ff : :class:`rdkit.ForceField`
            Forcefield to apply constraints to. Generally use UFF.

        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Add a very weak force constraint on all metal-metal
        # distances.
        #
        # for atoms in combinations(metal_atoms, r=2):
        #     ff.UFFAddDistanceConstraint(
        #         idx1=atoms[0].id,
        #         idx2=atoms[1].id,
        #         relative=True,
        #         minLen=0.8,
        #         maxLen=0.8,
        #         forceConstant=0.25e1
        #     )

        # Add constraints to UFF to hold metal geometry in place.
        for bond in metal_bonds:
            idx1 = bond.atom1.id
            idx2 = bond.atom2.id
            # Add distance constraints in place of metal bonds.
            # Target distance set to a given metal_binder_distance.
            ff.UFFAddDistanceConstraint(
                idx1=idx1,
                idx2=idx2,
                relative=False,
                minLen=self._metal_binder_distance,
                maxLen=self._metal_binder_distance,
                forceConstant=self._metal_binder_fc
            )

        # Also implement angular constraints to all atoms in the
        # metal complex.
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
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_positions(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=idx3,
                relative=False,
                minAngleDeg=np.degrees(angle),
                maxAngleDeg=np.degrees(angle),
                forceConstant=1.0e5
            )

    def apply_orientation_restriction(
        self,
        ff,
        mol,
        metal_bonds,
        metal_atoms
    ):
        """
        Applies UFF relative orientation restrcitions.

        Parameters
        ----------
        ff : :class:`rdkit.ForceField`
            Forcefield to apply constraints to. Generally use UFF.

        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        Returns
        -------
        None : :class:`NoneType`

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
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
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

        # Apply an angular constraint on the binder-metal-metal atoms
        # to maintain the metal centres relative orientation.
        for bonds in combinations(metal_bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.atom1, bond1.atom2]
            bond2_atoms = [bond2.atom1, bond2.atom2]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) < 4:
                continue

            if bond1.atom1 in metal_atoms:
                idx1 = bond1.atom1.id
                idx2 = bond1.atom2.id
            else:
                idx1 = bond1.atom2.id
                idx2 = bond1.atom1.id

            if bond2.atom1 in metal_atoms:
                idx3 = bond2.atom1.id
            else:
                idx3 = bond2.atom2.id

            pos1 = [
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_positions(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=idx3,
                relative=False,
                minAngleDeg=np.degrees(angle),
                maxAngleDeg=np.degrees(angle),
                forceConstant=1.0e4
            )

    def get_metal_atoms(self, mol):
        metal_atoms = []
        for atom in mol.atoms:
            if atom.atomic_number in self.metal_a_no:
                metal_atoms.append(atom)

        return metal_atoms

    def get_metal_bonds(self, mol, metal_atoms):
        metal_bonds = []
        ids_to_metals = []
        for bond in mol.bonds:
            if bond.atom1 in metal_atoms:
                metal_bonds.append(bond)
                ids_to_metals.append(bond.atom2.id)
            elif bond.atom2 in metal_atoms:
                metal_bonds.append(bond)
                ids_to_metals.append(bond.atom1.id)

        return metal_bonds, ids_to_metals

    def has_h(self, bond):
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

    def has_M(self, bond, metal_atoms):
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


class GulpMetalOptimizer(MetalOptimizer):
    """
    Applies forcefield optimizers that can handle metal centres.

    Notes
    -----
    By default, :meth:`optimize` will run a restricted optimization
    using constraints and the UFF4MOF. This forcefield requires some
    explicit metal atom definitions, which are determined by the user.

    Attributes
    ----------


    Examples
    --------

    """

    def __init__(
        self,
        gulp_path,
        metal_FF,
        output_dir=None,
        use_cache=False
    ):
        """
        Initialize a :class:`GulpMetalOptimizer` instance.

        Parameters
        ----------

        """
        self._gulp_path = gulp_path
        self._metal_FF = metal_FF
        self._output_dir = output_dir
        self.metal_a_no = list(range(21, 31))
        self.metal_a_no += list(range(39, 49))+list(range(72, 81))

        super(MetalOptimizer, self).__init__(use_cache=use_cache)

    def _add_atom_charge_flags(self, atom, atomkey):
        """


        https://github.com/rdkit/rdkit/blob/master/Code/
        GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp
        """
        total_valence = rdkit.Atom.GetTotalValence(atom)
        atnum = int(atom.GetAtomicNum())

        # Go through element cases.
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
            hybrid = rdkit.Atom.GetHybridization(atom)
            if hybrid != rdkit.HybridizationType.SP2:
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

    def _get_atom_label(self, atom):
        """


        https://github.com/rdkit/rdkit/blob/master/Code/
        GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp
        """
        atnum = int(atom.GetAtomicNum())
        atomkey = atom.GetSymbol()
        if len(atomkey) == 1:
            atomkey += '_'

        table = rdkit.GetPeriodicTable()

        chk1 = (
            rdkit.PeriodicTable.GetDefaultValence(table, atnum) == -1
        )
        chk2 = (rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 1)
        chk3 = (rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 7)
        chk4 = chk2 and chk3
        if chk1 or chk4:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if atnum == 84:
                atomkey += '3'
                if hybrid != rdkit.HybridizationType.SP3:
                    sys.exit('warning 84')
            elif atnum == 80:
                atomkey += '1'
                if hybrid != rdkit.HybridizationType.SP:
                    sys.exit('warning 80')
            else:
                if hybrid == rdkit.HybridizationType.SP:
                    atomkey += '1'
                elif hybrid == rdkit.HybridizationType.SP2:
                    chk1a = rdkit.Atom.GetIsAromatic(atom)
                    bonds = rdkit.Atom.GetBonds(atom)
                    conjugated = False
                    for bond in bonds:
                        if rdkit.Bond.GetIsConjugated(bond):
                            conjugated = True
                            break
                    chk2a = conjugated
                    chk3a = atnum in [6, 7, 8, 16]
                    chk4a = (chk1a or chk2a)
                    if chk4a and chk3a:
                        atomkey += 'R'
                    else:
                        atomkey += '2'
                elif hybrid == rdkit.HybridizationType.SP3:
                    atomkey += '3'
                elif hybrid == rdkit.HybridizationType.SP3D:
                    atomkey += '5'
                elif hybrid == rdkit.HybridizationType.SP3D2:
                    atomkey += '6'
                else:
                    sys.exit('final warning. not recog')
        atomkey = self._add_atom_charge_flags(atom, atomkey)
        return atomkey

    def _type_translator(self):
        type_translator = {}
        types = set([self.atom_labels[i][0] for i in self.atom_labels])
        for t in types:
            if not t[1].isalpha():
                symb = t[0]
            else:
                symb = t[0:2]
            for i in range(1, 100):
                name = f'{symb}{i}'
                if name in type_translator.values():
                    continue
                else:
                    type_translator[t] = name
                    break

        return type_translator

    def _position_section(self, mol, type_translator):
        position_section = '\ncartesian\n'
        for atom in mol.atoms:
            atom_type = type_translator[self.atom_labels[atom.id][0]]
            position = mol.get_center_of_mass(atom_ids=[atom.id])
            posi_string = (
                f'{atom_type} core {round(position[0], 5)} '
                f'{round(position[1], 5)} {round(position[2], 5)}\n'
            )
            position_section += posi_string

        return position_section

    def _bond_section(self, mol, metal_atoms):
        bond_section = '\n'
        for bond in mol.bonds:
            atom_types = [
                self.atom_labels[i.id][0]
                for i in [bond.atom1, bond.atom2]
            ]

            # Set bond orders.
            if self.has_h(bond):
                # H has bond order of 1.
                bond_type = ''
            elif self.has_M(bond, metal_atoms):
                bond_type = 'half'
            elif '_R' in atom_types[0] and '_R' in atom_types[1]:
                bond_type = 'resonant'
            elif bond.order == 1:
                bond_type = ''
            elif bond.order == 2:
                bond_type = 'double'
            elif bond.order == 3:
                bond_type = 'triple'

            string = (
                f'connect {bond.atom1.id+1} {bond.atom2.id+1} '
                f'{bond_type}'
            )
            bond_section += string+'\n'

        return bond_section

    def _species_section(self, type_translator):
        species_section = '\nspecies\n'
        for spec in type_translator:
            name = type_translator[spec]
            species_section += f'{name} {spec}\n'

        return species_section

    def _write_gulp_file(self, mol, metal_atoms, in_file, output_xyz):

        type_translator = self._type_translator()

        top_line = 'opti conv noautobond fix molmec cartesian\n'

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)

        library = '\nlibrary uff4mof.lib\n'

        output_section = (
            '\n'
            f'output xyz {output_xyz}\n'
            # 'output movie xyz steps_.xyz\n'
        )

        with open(in_file, 'w') as f:
            f.write(top_line)
            f.write(position_section)
            f.write(bond_section)
            f.write(species_section)
            f.write(library)
            f.write(output_section)

    def assign_FF(self, mol):
        """
        Assign forcefield types to molecule.

        Parameters
        ----------

        Returns
        -------

        None : :class:`NoneType`

        """
        metal_atoms = self.get_metal_atoms(mol)
        metal_ids = [i.id for i in metal_atoms]
        metal_bonds, _ = self.get_metal_bonds(mol, metal_atoms)
        edit_mol = self.to_rdkit_mol_no_metals(
            mol=mol,
            metal_atoms=metal_atoms,
            metal_bonds=metal_bonds
        )

        # Get forcefield parameters.
        rdkit.SanitizeMol(edit_mol)
        self.atom_labels = {}
        func_group_atoms = {
            j.id: i.fg_type.name
            for i in mol.func_groups
            for j in i.bonders
        }
        for i in range(edit_mol.GetNumAtoms()):
            if i in metal_ids:
                self.atom_labels[i] = [None, 'metal', None]
            else:
                atom = edit_mol.GetAtomWithIdx(i)
                atom_label = self._get_atom_label(atom)
                if i in func_group_atoms:
                    fg = func_group_atoms[i]
                    self.atom_labels[i] = [atom_label, 'bonder', fg]
                else:
                    self.atom_labels[i] = [atom_label, None, None]

        # Write UFF4MOF specific forcefield parameters.
        # Metals.
        for atomid in self.atom_labels:
            if self.atom_labels[atomid][1] == 'metal':
                self.atom_labels[atomid][0] = self._metal_FF

        # Metal binder atoms of specific forcefields.
        # Check functional groups.

    def _run_gulp(self, in_file, out_file):
        cmd = f'{self._gulp_path} < {in_file}'
        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def _move_generated_files(self, files):
        if not os.path.exists(self._output_dir):
            os.mkdir(self._output_dir)

        for file in files:
            os.rename(file, f'{self._output_dir}/{file}')

    def extract_final_energy(self, file):
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        with open(file, 'r') as f:
            for line in f.readlines():
                if 'Final energy =' in line:
                    string = nums.search(line.rstrip()).group(0)
                    return float(string)

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
        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        in_file = 'gulp_opt.gin'
        out_file = 'gulp_opt.ginout'
        output_xyz = 'gulp_opt.xyz'

        metal_atoms = self.get_metal_atoms(mol)

        # Write GULP file.
        self._write_gulp_file(
            mol=mol,
            metal_atoms=metal_atoms,
            in_file=in_file,
            output_xyz=output_xyz
        )

        # Run.
        self._run_gulp(in_file, out_file)

        # Update from output.
        mol.update_from_file(output_xyz)

        # Move files.
        self._move_generated_files(
            files=[in_file, out_file, output_xyz]
        )


class GulpMDMetalOptimizer(GulpMetalOptimizer):
    """
    Applies forcefield MD that can handle metal centres.

    Notes
    -----
    By default, :meth:`optimize` will run a restricted optimization
    using constraints and the UFF4MOF. This forcefield requires some
    explicit metal atom definitions, which are determined by the user.

    Attributes
    ----------


    Examples
    --------

    """

    def __init__(
        self,
        gulp_path,
        metal_FF,
        output_dir=None,
        integrator='stochastic',
        ensemble='nvt',
        temperature='300',
        equilbration='1.0',
        production='10.0',
        timestep='1.0',
        N_conformers=10,
        use_cache=False
    ):
        """
        Initialize a :class:`GulpMetalOptimizer` instance.

        Parameters
        ----------

        """
        self._gulp_path = gulp_path
        self._metal_FF = metal_FF
        self._output_dir = output_dir
        self._integrator = integrator
        self._ensemble = ensemble
        self._temperature = temperature
        self._equilbration = equilbration
        self._production = production
        self._timestep = timestep
        self._N_conformers = N_conformers
        samples = float(self._production) / float(self._N_conformers)
        self._sample = samples
        self._write = samples

        self.metal_a_no = list(range(21, 31))
        self.metal_a_no += list(range(39, 49))+list(range(72, 81))

        super(MetalOptimizer, self).__init__(use_cache=use_cache)

    def _write_gulp_file(
        self,
        mol,
        metal_atoms,
        in_file,
        output_traj
    ):

        type_translator = self._type_translator()

        top_line = 'md conv noautobond fix molmec cartesian\n'

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)

        library = '\nlibrary uff4mof.lib\n'

        md_section = (
            '\n'
            "mdmaxtemp 10000\n"
            f"integrator {self._integrator}\n"
            f"ensemble {self._ensemble}\n"
            f"temperature {self._temperature}\n"
            f"equilbration {self._equilbration} ps\n"
            f"production {self._production} ps\n"
            f"timestep {self._timestep} fs\n"
            f"sample {self._sample} ps\n"
            f"write {self._write} ps\n"
            '\n'
        )

        output_section = (
            '\n'
            f'output trajectory ascii {output_traj}\n'
        )

        with open(in_file, 'w') as f:
            f.write(top_line)
            f.write(position_section)
            f.write(bond_section)
            f.write(species_section)
            f.write(library)
            f.write(md_section)
            f.write(output_section)

    def _convert_traj_to_xyz(self, output_xyz, output_traj, xyz_traj):

        # Get atom types from an existing xyz file.
        atom_types = []
        with open(output_xyz, 'r') as f:
            for line in f.readlines()[2:]:
                atom_types.append(line.rstrip().split(' ')[0])

        # Read in lines from trajectory file.
        with open(output_traj, 'r') as f:
            lines = f.readlines()

        # Split file using strings.
        id = 0
        s_times = {}
        s_coords = {}
        s_vels = {}
        s_derivs = {}
        s_sites = {}
        coords = []
        vels = []
        derivs = []
        sites = []
        switch = None
        for line in lines:
            li = line.rstrip()
            if '#  Time/KE/E/T' in li:
                switch = 'times'
                s_coords[id] = coords
                coords = []
                s_vels[id] = vels
                vels = []
                s_derivs[id] = derivs
                derivs = []
                s_sites[id] = sites
                sites = []
                id += 1
            elif '#  Coordinates' in li:
                switch = 'coords'
            elif '#  Velocities' in li:
                switch = 'vels'
            elif '#  Derivatives' in li:
                switch = 'derivs'
            elif '#  Site energies' in li:
                switch = 'sites'
            elif switch == 'coords':
                coords.append([i for i in li.split(' ') if i])
            elif switch == 'vels':
                vels.append([i for i in li.split(' ') if i])
            elif switch == 'derivs':
                derivs.append([i for i in li.split(' ') if i])
            elif switch == 'sites':
                sites.append([i for i in li.split(' ') if i])
            elif switch == 'times':
                s_times[id] = [i for i in li.split(' ') if i]
            elif switch is None:
                pass
        # Add final timestep.
        s_coords[id] = coords
        coords = []
        s_vels[id] = vels
        vels = []
        s_derivs[id] = derivs
        derivs = []
        s_sites[id] = sites
        sites = []

        ids = []
        tt = []
        pot_energies = []
        new_lines = []
        for id in s_times:
            times = s_times[id]
            ids.append(id)
            tt.append(float(times[0]))
            pot_energies.append(float(times[2]))
            coords = s_coords[id]
            sites = s_sites[id]
            xyz_string = (
                f'{len(coords)}\n'
                f'{times[0]},{times[1]},{times[2]},{times[3]}\n'
            )

            for i, coord in enumerate(coords):
                site = sites[i][0]
                xyz_string += (
                    f'{atom_types[i]} {coord[0]} {coord[1]} '
                    f'{coord[2]} {site}\n'
                )

            new_lines.append(xyz_string)

        with open(xyz_traj, 'w') as f:
            for line in new_lines:
                f.write(line)

        return atom_types, ids, tt, pot_energies, s_times, s_coords

    def _write_conformer_xyz_file(
        self,
        id,
        filename,
        s_times,
        s_coords,
        atom_types
    ):
        times = s_times[id]
        coords = s_coords[id]
        xyz_string = (
            f'{len(coords)}\n'
            f'{times[0]},{times[1]},{times[2]},{times[3]}\n'
        )
        for i, coord in enumerate(coords):
            xyz_string += (
                f'{atom_types[i]} {coord[0]} {coord[1]} {coord[2]}\n'
            )
        with open(filename, 'w') as f:
            f.write(xyz_string)

    def _save_lowest_energy_conf(
        self,
        mol,
        output_xyz,
        output_traj,
        xyz_traj,
        low_conf_xyz
    ):

        # Convert GULP trajectory file to xyz trajectory.
        results = self._convert_traj_to_xyz(
            output_xyz,
            output_traj,
            xyz_traj
        )
        atom_types, ids, tt, pot_energies, s_times, s_coords = results

        # Optimise all conformers.
        min_energy = 1E10
        for id in ids:
            self._write_conformer_xyz_file(
                id=id,
                filename='temp_conf.xyz',
                s_times=s_times,
                s_coords=s_coords,
                atom_types=atom_types
            )
            mol.update_from_file('temp_conf.xyz')
            gulp_opt = GulpMetalOptimizer(
                gulp_path=self._gulp_path,
                metal_FF=self._metal_FF,
                output_dir=self._output_dir
            )
            gulp_opt.assign_FF(mol)
            gulp_opt.optimize(mol=mol)
            energy = gulp_opt.extract_final_energy(
                file=join(self._output_dir, 'gulp_opt.ginout')
            )
            print(energy, min_energy)

            if energy < min_energy:
                # Find lowest energy conformation and output to XYZ.
                print(min_energy)
                min_energy = energy
                self._write_conformer_xyz_file(
                    id=id,
                    filename=low_conf_xyz,
                    s_times=s_times,
                    s_coords=s_coords,
                    atom_types=atom_types
                )

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
        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        in_file = 'gulp_MD.gin'
        out_file = 'gulp_MD.ginout'
        output_xyz = 'gulp_MD_template.xyz'
        mol.write(output_xyz)
        output_traj = 'gulp_MD.trg'
        xyz_traj = 'gulp_MD_traj.xyz'
        low_conf_xyz = 'low_energy_conf.xyz'

        metal_atoms = self.get_metal_atoms(mol)

        # Write GULP file.
        self._write_gulp_file(
            mol=mol,
            metal_atoms=metal_atoms,
            in_file=in_file,
            output_traj=output_traj
        )

        # Run.
        self._run_gulp(in_file, out_file)

        # Get lowest energy conformer from trajectory.
        self._save_lowest_energy_conf(
            mol=mol,
            output_xyz=output_xyz,
            output_traj=output_traj,
            xyz_traj=xyz_traj,
            low_conf_xyz=low_conf_xyz
        )

        # Update from output.
        mol.update_from_file(low_conf_xyz)

        # Move files.
        self._move_generated_files(
            files=[
                in_file, out_file, output_xyz,
                output_traj, xyz_traj, low_conf_xyz
            ]
        )
