"""
Defines optimizers for metallo-architectures.

"""

from itertools import combinations
import rdkit.Chem.AllChem as rdkit
import logging
import numpy as np

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
    _scale : :class:`float`
        Distance to place ligand binder atoms from metal.

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

    _prearrange : :class:`bool`
        `True` to prearrange functional groups around metal centre.

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
            restrict_orientation=True,
            prearrange=True,
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

    def __init__(self, scale, force_constant, rel_distance,
                 res_steps, prearrange=True, restrict_all_bonds=False,
                 restrict_orientation=False,
                 use_cache=False):
        """
        Initialize a :class:`MetalOptimizer` instance.

        Parameters
        ----------
        scale : :class:`float`
            Distance to place ligand binder atoms from metal.

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

        prearrange : :class:`bool`, optional
            `True` to prearrange functional groups around metal centre.

        restrict_all_bonds : :class:`bool`, optional
            `True` to restrict all bonds except for ligand-FG bonds.

        restrict_orientation : :class:`bool`, optional
            `True` to restrict metal complex FG angles relative to
            topology centre of mass.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """
        self._scale = scale
        self._force_constant = force_constant
        self._rel_distance = rel_distance
        self._prearrange = prearrange
        self._restrict_all_bonds = restrict_all_bonds
        self._restrict_orientation = restrict_orientation
        self._res_steps = res_steps

        self.metal_a_no = list(range(21, 31))
        self.metal_a_no += list(range(39, 49))+list(range(72, 81))

        super(MetalOptimizer, self).__init__(use_cache=use_cache)

    def prearrange_fgs(self, mol):
        """
        Prearrange the metal interacting ligand functional groups.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`ValueError`
            If the metal functional group has two bonder atoms, this
            code will fail and raise this error. This assumption
            implies that each metal-interacting functional group makes
            a single metal-ligand bond.

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

    def get_input_constraints(self, mol, ids_to_metals, metal_atoms,
                              include_bonders=False):
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
            constraints[bond] = {
                'idx1': idx1,
                'idx2': idx2,
                'type': 'bond',
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

    def _restricted_optimization(self, mol, metal_atoms,
                                 metal_bonds,
                                 ids_to_metals,
                                 input_constraints=None):
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

        # For bonds between ligand bonders and the rest of the liagnd,
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
                    ff.UFFAddDistanceConstraint(
                        idx1=idx1,
                        idx2=idx2,
                        relative=True,
                        minLen=self._rel_distance,
                        maxLen=self._rel_distance,
                        forceConstant=self._force_constant
                    )

        # Perform UFF optimization with rdkit.
        ff.Minimize(maxIts=25)

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
        metal_atoms = []
        for atom in mol.atoms:
            if atom.atomic_number in self.metal_a_no:
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

        # Second step is to perform a forcefield optimisation that
        # only optimises non metal atoms that are not bonded to the
        # metal.

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
                metal_atoms=metal_atoms,
                metal_bonds=metal_bonds,
                ids_to_metals=ids_to_metals,
                input_constraints=input_constraints
            )

    def apply_metal_centre_constraints(self, mol, ff, metal_bonds,
                                       metal_atoms):
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
            # Target distance set to _scale.
            ff.UFFAddDistanceConstraint(
                idx1=idx1,
                idx2=idx2,
                relative=False,
                minLen=self._scale,
                maxLen=self._scale,
                forceConstant=self._force_constant
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
                minAngleDeg=np.degrees(angle)-2,
                maxAngleDeg=np.degrees(angle)+2,
                forceConstant=1.0e4
            )

    def apply_orientation_restriction(self, ff, mol, metal_bonds,
                                      metal_atoms):
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
