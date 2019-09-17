"""
Optimizers
==========

#. :class:`.NullOptimizer`
#. :class:`.MMFF`
#. :class:`.UFF`
#. :class:`.ETKDG`
#. :class:`.XTB`
#. :class:`.MacroModelForceField`
#. :class:`.MacroModelMD`
#. :class:`.MOPAC`
#. :class:`.OptimizerSequence`
#. :class:`.CageOptimizerSequence`
#. :class:`.TryCatchOptimizer`
#. :class:`.RaisingOptimizer`
#. :class:`.MetalOptimizer`


Optimizers are objects used to optimize molecules. Each optimizer is
initialized with some settings and can optimize a molecule
with :meth:`~.Optimizer.optimize`.

.. code-block:: python

    import stk

    mol = stk.BuildingBlock('NCCCN', ['amine'])
    mmff = stk.MMFF()
    mmff.optimize(mol)

    # Optimizers also work with ConstructedMolecule objects.
    polymer = stk.ConstructedMolecule(
        building_blocks=[mol],
        topology_graph=stk.polymer.Linear('A', [0], n=3)
    )
    etkdg = stk.ETKDG()
    etkdg.optimize(polymer)

Sometimes it is desirable to chain multiple optimizations, one after
another. For example, before running an optimization, it may be
desirable to embed a molecule first, to generate an initial structure.
:class:`.OptimizerSequence` may be used for this.

.. code-block:: python

    # Create a new optimizer which chains the previously defined
    # mmff and etkdg optimizers.
    optimizer_sequence = stk.OptimizerSequence(etkdg, mmff)

    # Run each optimizer in sequence.
    optimizer_sequence.optimize(polymer)

By default, running :meth:`.Optimizer.optimize` twice on the same
molecule will perform an optimization a second time on a molecule. If
we want to skip optimizations on molecules which have already been
optimized we can use the `use_cache` flag.

.. code-block:: python

    caching_etkdg = stk.ETKDG(use_cache=True)
    # First optimize call runs an optimization.
    caching_etkdg.optimize(polymer)
    # Second call does nothing.
    caching_etkdg.optimize(polymer)

Caching is done on a per :class:`.Optimizer` basis. Just because the
molecule has been cached by one :class:`.Optimizer` instance does not
mean that a different :class:`.Optimizer` instance will no longer
optimize the molecule.

.. _`adding optimizers`:

Making New Optimizers
---------------------

New optimizers can be made by simply making a class which inherits the
:class:`.Optimizer` class. In addition to this, the new class must
define a :meth:`~.Optimizer.optimize` method. The method must take 1
mandatory `mol` parameter. :meth:`~.Optimizer.optimize` will take the
`mol` and change its structure in whatever way it likes. Beyond this
there are no requirements. New optimizers can be added into the
:mod:`.optimizers` submodule or into a new submodule.

"""

import logging
import numpy as np
import rdkit.Chem.AllChem as rdkit
import warnings
from itertools import combinations
import os
from functools import wraps
import subprocess as sp
import uuid
import shutil
from ...utilities import (
    is_valid_xtb_solvent,
    XTBInvalidSolventError,
    XTBExtractor,
    vector_angle,
    get_dihedral
)
import pywindow

logger = logging.getLogger(__name__)


def _add_cache_use(optimize):
    """
    Make :meth:`~Optimizer.optimize` use the :attr:`~Optimizer._cache`.

    Decorates `optimize` so that before running it checks if the
    :class:`.Molecule` has already been optimized by the
    optimizer. If so, and :attr:`~Optimizer.use_cache` is ``True``,
    then the molecule is skipped and no optimization is performed.

    Parameters
    ----------
    optimize : :class:`function`
        A function which is to have skipping added to it.

    Returns
    -------
    :class:`function`
        The decorated function.

    """

    @wraps(optimize)
    def inner(self, mol):
        if self._use_cache and mol in self._cache:
            logger.info(f'Skipping optimization on {mol}.')
        else:
            optimize(self, mol)
            if self._use_cache:
                self._cache.add(mol)

    return inner


class Optimizer:
    """
    A base class for optimizers.

    """

    def __init__(self, use_cache=False):
        """
        Initialize an :class:`Optimizer`.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        # Holds every previously optimized molecule if use_cache is
        # true.
        self._cache = set()
        self._use_cache = use_cache

    def __init_subclass__(cls, **kwargs):
        cls.optimize = _add_cache_use(cls.optimize)
        return super().__init_subclass__(**kwargs)

    def is_caching(self):
        """
        ``True`` if the optimizer has caching turned on.

        Returns
        -------
        :class:`bool`
            ``True`` if the optimizer has caching turned on.

        """

        return self._use_cache

    def add_to_cache(self, mol):
        """
        Add a molecule to the cache.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be added to the cache.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._cache.add(mol)

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

        raise NotImplementedError()


class OptimizerSequence(Optimizer):
    """
    Applies optimizers in sequence.

    Examples
    --------
    Let's say we want to embed a molecule with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        optimizer = stk.OptimizerSequence(stk.ETKDG(), stk.MMFF())
        optimizer.optimize(mol)

    """

    def __init__(self, *optimizers, use_cache=False):
        """
        Initialize a :class:`OptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._optimizers = optimizers
        super().__init__(use_cache=use_cache)

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

        for optimizer in self._optimizers:
            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol}".')
            optimizer.optimize(mol)


class CageOptimizerSequence(Optimizer):
    """
    Applies :class:`Optimizer` objects to a cage.

    Before each :class:`Optimizer` in the sequence is applied to the
    cage, it is checked to see if it is collapsed. If it is
    collapsed, the optimization sequence ends immediately.

    Examples
    --------
    Let's say we want to embed a cage with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCNCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
        cage = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.cage.FourPlusSix()
        )
        optimizer = stk.CageOptimizerSequence(stk.ETKDG(), stk.MMFF())
        optimizer.optimize(cage)

    """

    def __init__(self, *optimizers, use_cache=False):
        """
        Initialize a :class:`CageOptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            The :class:`Optimizers` used in sequence to optimize
            cage molecules.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._optimizers = optimizers
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The cage to be optimized. It must have an attribute
            :attr:`num_windows`, which contains the expected number
            of windows if the cage is not collapsed.

        Returns
        -------
        None : :class:`NoneType`

        """

        for optimizer in self._optimizers:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                loader = pywindow.molecular.Molecule.load_rdkit_mol
                pw_molecule = loader(mol.to_rdkit_mol())
                windows = pw_molecule.calculate_windows()
            logger.debug(f'Windows found: {windows}.')

            if windows is None or len(windows) != mol.num_windows:
                logger.info(f'"{mol}" is collapsed, exiting early.')
                return

            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol}".')
            optimizer.optimize(mol)


class NullOptimizer(Optimizer):
    """
    Does not perform optimizations.

    """

    def optimize(self, mol):
        """
        Do not optimize `mol`.

        This function just returns immediately without changing the
        molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        return


class TryCatchOptimizer(Optimizer):
    """
    Try to optimize with a :class:`Optimizer`, use another on failure.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create some molecules to optimize.
        mol1 = stk.BuildingBlock('NCCN', ['amine'])
        mol2 = stk.BuildingBlock('CCCCC')
        mol3 = stk.BuildingBlock('O=CCCN')

        # Create an optimizer which may fail.
        uff = stk.UFF()

        # Create a backup optimizer.
        mmff = stk.MMFF()

        # Make an optimizer which tries to run raiser and if that
        # raises an error, will run mmff on the molecule instead.
        try_catch = stk.TryCatchOptimizer(
            try_optimizer=uff,
            catch_optimizer=mmff
        )

        # Optimize the molecules. In each case if the optimization with
        # UFF fails, MMFF is used to optimize the molecule instead.
        try_catch.optimize(mol1)
        try_catch.optimize(mol2)
        try_catch.optimzie(mol3)

    """

    def __init__(
        self,
        try_optimizer,
        catch_optimizer,
        use_cache=False
    ):
        """
        Initialize a :class:`TryCatchOptimizer` instance.

        Parameters
        ----------
        try_optimizer : :class:`Optimizer`
            The optimizer which is used initially to try and optimize a
            :class:`.Molecule`.

        catch_optimizer : :class:`Optimizer`
            If :attr:`try_optimizer` raises an error, this optimizer is
            run on the :class:`.Molecule` instead.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._try_optimizer = try_optimizer
        self._catch_optimizer = catch_optimizer
        super().__init__(use_cache=use_cache)

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

        try:
            return self._try_optimizer.optimize(mol)
        except Exception:
            try_name = self._try_optimizer.__class__.__name__
            catch_name = self._catch_optimizer.__class__.__name__
            logger.error(
                f'{try_name} failed, trying {catch_name}.',
                exc_info=True
            )
            return self._catch_optimizer.optimize(mol)


class RaisingOptimizerError(Exception):
    ...


class RaisingOptimizer(Optimizer):
    """
    Raises and optimizes at random.

    This optimizer is used for debugging to simulate optimizations
    which sometimes complete successfully and sometimes
    randomly fail.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        partial_raiser = stk.RaisingOptimizer(
            optimizer=stk.ETKDG(),
            fail_chance=0.75
        )
        # 75 % chance an error will be raised by calling optimize.
        partial_raiser.optimize(mol)

    """

    def __init__(self, optimizer, fail_chance=0.5, use_cache=False):
        """
        Initialize :class:`PartialRaiser`.

        Parameters
        ----------
        optimizer : :class:`Optimizer`
            When the optimizer does not fail, it uses this
            :class:`Optimizer` to optimize molecules.

        fail_chance : :class:`float`, optional
            The probability that the optimizer will raise an error each
            time :meth:`optimize` is used.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._optimizer = optimizer
        self._fail_chance = fail_chance
        super().__init__(use_cache=use_cache)

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

        Raises
        ------
        :class:`RaisingOptimizerError`
            This error is raised at random.

        """

        if np.random.rand() < self._fail_chance:
            raise RaisingOptimizerError('Used RaisingOptimizer.')
        return self._optimizer.optimize(mol)


class MMFF(Optimizer):
    """
    Use the MMFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        mmff = stk.MMFF()
        mmff.optimize(mol)

    """

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

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.MMFFOptimizeMolecule(rdkit_mol)
        mol.update_from_rdkit_mol(rdkit_mol)


class UFF(Optimizer):
    """
    Use the UFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        uff = stk.UFF()
        uff.optimize(mol)

    """

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

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.UFFOptimizeMolecule(rdkit_mol)
        mol.update_from_rdkit_mol(rdkit_mol)


class ETKDG(Optimizer):
    """
    Uses the ETKDG [#]_ v2 algorithm to find an optimized structure.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        etkdg = stk.ETKDG()
        etkdg.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, random_seed=12, use_cache=False):
        """
        Initialize a :class:`ETKDG` instance.

        Parameters
        ----------
        random_seed : :class:`int`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._random_seed = random_seed
        super().__init__(use_cache=use_cache)

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

        params = rdkit.ETKDGv2()
        params.clearConfs = True
        params.random_seed = self._random_seed

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.EmbedMolecule(rdkit_mol, params)
        mol.set_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )


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
            metal_atoms
        )

        # First step is to pre-arrange the metal centre based on the
        # MetalComplex topology.
        if self._prearrange:
            self.prearrange_fgs(mol)

        for i in range(self._res_steps):
            print(i)
            self.restricted_optimization(
                mol=mol,
                metal_atoms=metal_atoms,
                metal_bonds=metal_bonds,
                ids_to_metals=ids_to_metals,
                rel_distance=0.9,
                force_constant=1e2,
                input_constraints=input_constraints
            )
            mol.write(f'opt{i}.mol')

    def get_input_constraints(self, mol, ids_to_metals, metal_atoms):
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
            # # Do not restrict angles of bonds to atoms bonded to
            # # metals.
            # bonded_to_metals = [
            #     idx1 in ids_to_metals,
            #     idx2 in ids_to_metals,
            #     idx3 in ids_to_metals
            # ]
            # if any(bonded_to_metals):
            #     continue
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


class XTBOptimizerError(Exception):
    ...


class XTBConvergenceError(XTBOptimizerError):
    ...


class XTB(Optimizer):
    """
    Uses GFN-xTB [1]_ to optimize molecules.

    Notes
    -----
    When running :meth:`optimize`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised, so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.XTB` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Furthermore, :meth:`optimize` will check that the
    structure is adequately optimized by checking for negative
    frequencies after a Hessian calculation. `max_runs` can be
    provided to the initializer to set the maximum number of
    optimizations which will be attempted at the given
    `opt_level` to obtain an optimized structure. However, we outline
    in the examples how to iterate over `opt_levels` to increase
    convergence criteria and hopefully obtain an optimized structure.
    The presence of negative frequencies can occur even when the
    optimization has converged based on the given `opt_level`.

    Attributes
    ----------
    incomplete : :class:`set` of :class:`.Molecule`
        A :class:`set` of molecules passed to :meth:`optimize` whose
        optimzation was incomplete.

    Examples
    --------
    Note that for :class:`.ConstructedMolecule` objects constructed by
    ``stk``, :class:`XTB` should usually be used in a
    :class:`.OptimizerSequence`. This is because xTB only uses
    xyz coordinates as input and so will not recognize the long bonds
    created during construction. An optimizer which can minimize
    these bonds should be used before :class:`XTB`.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCNCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCCC=O', ['aldehyde'])
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear("AB", [0, 0], 3)
        )

        xtb = stk.OptimizerSequence(
            stk.UFF(),
            stk.XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
        )
        xtb.optimize(polymer)

    By default, all optimizations with xTB are performed using the
    ``--ohess`` flag, which forces the calculation of a numerical
    Hessian, thermodynamic properties and vibrational frequencies.
    :meth:`optimize` will check that the structure is appropriately
    optimized (i.e. convergence is obtained and no negative vibrational
    frequencies are present) and continue optimizing a structure (up to
    `max_runs` times) until this is achieved. This loop, by
    default, will be performed at the same `opt_level`. The
    following example shows how a user may optimize structures with
    tigher convergence criteria (i.e. different `opt_level`)
    until the structure is sufficiently optimized. Furthermore, the
    calculation of the Hessian can be turned off using
    `max_runs` to ``1`` and `calculate_hessian` to ``False``.

    .. code-block:: python

        # Use crude optimization with max_runs=1 because this will
        # not achieve optimization and rerunning it is unproductive.
        xtb_crude = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_crude',
            unlimited_memory=True,
            opt_level='crude',
            max_runs=1,
            calculate_hessian=True
        )
        # Use normal optimization with max_runs == 2.
        xtb_normal = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_normal',
            unlimited_memory=True,
            opt_level='normal',
            max_runs=2
        )
        # Use vtight optimization with max_runs == 2, which should
        # achieve sufficient optimization.
        xtb_vtight = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_vtight',
            unlimited_memory=True,
            opt_level='vtight',
            max_runs=2
        )

        optimizers = [xtb_crude, xtb_normal, xtb_vtight]
        for optimizer in optimizers:
            optimizer.optimize(polymer)
            if polymer not in optimizer.incomplete:
                break

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/setup.html

    """

    def __init__(
        self,
        xtb_path,
        gfn_version=2,
        output_dir=None,
        opt_level='normal',
        max_runs=2,
        calculate_hessian=True,
        num_cores=1,
        electronic_temperature=300,
        solvent=None,
        solvent_grid='normal',
        charge=0,
        num_unpaired_electrons=0,
        unlimited_memory=False,
        use_cache=False
    ):
        """
        Initialize a :class:`XTB` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`int`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        max_runs : :class:`int`, optional
            Maximum number of optimizations to attempt in a row.

        calculate_hessian : :class:`bool`, optional
            Toggle calculation of the hessian and vibrational
            frequencies after optimization. ``True`` is required to
            check that the structure is completely optimized.
            ``False`` will drastically speed up the calculation but
            potentially provide incomplete optimizations and forces
            :attr:`max_runs` to be ``1``.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        electronic_temperature : :class:`int`, optional
            Electronic temperature in Kelvin.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        charge : :class:`int`, optional
            Formal molecular charge.

        num_unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                raise XTBInvalidSolventError(
                    f'No solvent valid for version',
                    f' {gfn_version!r}.'
                )
            if not is_valid_xtb_solvent(gfn_version, solvent):
                raise XTBInvalidSolventError(
                    f'Solvent {solvent!r} is invalid for ',
                    f'version {gfn_version!r}.'
                )

        if not calculate_hessian and max_runs != 1:
            max_runs = 1
            logger.warning(
                'Requested that hessian calculation was skipped '
                'but the number of optimizations requested was '
                'greater than 1. The number of optimizations has been '
                'set to 1.'
            )

        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir
        self._opt_level = opt_level
        self._max_runs = max_runs
        self._calculate_hessian = calculate_hessian
        self._num_cores = str(num_cores)
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self.incomplete = set()
        super().__init__(use_cache=use_cache)

    def _has_neg_frequencies(self, output_file):
        """
        Check for negative frequencies.

        Parameters
        ----------
        output_file : :class:`str`
            Name of output file with xTB results.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if a negative frequency is present.

        """
        xtbext = XTBExtractor(output_file=output_file)
        # Check for one negative frequency, excluding the first
        # 6 frequencies.
        return any(x < 0 for x in xtbext.frequencies[6:])

    def _is_complete(self, output_file):
        """
        Check if xTB optimization has completed and converged.

        Parameters
        ----------
        output_file : :class:`str`
            Name of xTB output file.

        Returns
        -------
        :class:`bool`
            Returns ``False`` if a negative frequency is present.

        Raises
        -------
        :class:`XTBOptimizerError`
            If the optimization failed.

        :class:`XTBConvergenceError`
            If the optimization did not converge.

        """
        if output_file is None:
            # No simulation has been run.
            return False
        # If convergence is achieved, then .xtboptok should exist.
        if os.path.exists('.xtboptok'):
            # Check for negative frequencies in output file if the
            # hessian was calculated.
            # Return True if there exists at least one.
            if self._calculate_hessian:
                return not self._has_neg_frequencies(output_file)
            else:
                return True
        elif os.path.exists('NOT_CONVERGED'):
            raise XTBConvergenceError('Optimization not converged.')
        else:
            raise XTBOptimizerError('Optimization failed to complete')

    def _run_xtb(self, xyz, out_file):
        """
        Run GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        # Set optimization level and type.
        if self._calculate_hessian:
            # Do optimization and check hessian.
            optimization = f'--ohess {self._opt_level}'
        else:
            # Do optimization.
            optimization = f'--opt {self._opt_level}'

        if self._solvent is not None:
            solvent = f'--gbsa {self._solvent} {self._solvent_grid}'
        else:
            solvent = ''

        cmd = (
            f'{memory} {self._xtb_path} {xyz} '
            f'--gfn {self._gfn_version} '
            f'{optimization} --parallel {self._num_cores} '
            f'--etemp {self._electronic_temperature} '
            f'{solvent} --chrg {self._charge} '
            f'--uhf {self._num_unpaired_electrons}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Uses the shell if unlimited_memory is True to run
                # multiple commands in one subprocess.
                shell=self._unlimited_memory
            )

    def _run_optimizations(self, mol):
        """
        Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """
        for run in range(self._max_runs):
            xyz = f'input_structure_{run+1}.xyz'
            out_file = f'optimization_{run+1}.output'
            mol.write(xyz)
            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            coord_file = 'xtbhess.coord'
            coord_exists = os.path.exists(coord_file)
            output_xyz = 'xtbopt.xyz'
            opt_complete = self._is_complete(out_file)
            if not opt_complete:
                if coord_exists:
                    # The calculation is incomplete.
                    # Update mol from xtbhess.coord and continue.
                    mol.update_from_file(coord_file)
                else:
                    # Update mol from xtbopt.xyz.
                    mol.update_from_file(output_xyz)
                    # If the negative frequencies are small, then GFN
                    # may not produce the restart file. If that is the
                    # case, exit optimization loop and warn.
                    self.incomplete.add(mol)
                    logging.warning(
                        f'Small negative frequencies present in {mol}.'
                    )
                    return False
            else:
                # Optimization is complete.
                # Update mol from xtbopt.xyz.
                mol.update_from_file(output_xyz)
                break
        return opt_complete

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

        # Remove mol from self.incomplete if present.
        if mol in self.incomplete:
            self.incomplete.remove(mol)

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            complete = self._run_optimizations(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            self.incomplete.add(mol)
            logging.warning(f'Optimization is incomplete for {mol}.')
