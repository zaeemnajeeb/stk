"""
Disulfide Factory
============

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from ..functional_groups import Disulfide


class DisulfideFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Disulfide` instances.

    Creates functional groups from substructures, which match the
    ``[S]([*])[S]([*])`` functional group string.

    Examples
    --------
    You want to create a building block which has :class:`.Disulfide`
    functional groups. You want one of the sulfur atoms in the functional
    group to be the *bonder* atoms, and the remaining SR group to be a 
    leaving group.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='CCCSSC',
            functional_groups=(stk.DisulfideFactory(), ),
        )

    You want to create a building block which has :class:`.Disulfide`
    functional groups. You want the sulfur atoms to be the *bonder*
    atoms and the hydrogen atoms to be the *deleter* atoms.

    .. code-block:: python

        import stk

        disulfide_factory = stk.DisulfideFactory(
            # The indices of the sulfur atoms in the functional
            # group string (see docstring) are 0 and 2.
            bonders=(0, 2),
            # The indices of the R groups in the functional
            # group string (see docstring) are 1 and 3.
            deleters=(1, 3),
        )
        building_block = stk.BuildingBlock(
            smiles='CCCSSC',
            functional_groups=(disulfide_factory, ),
        )


    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(
        self,
        bonders=(0,),
        deleters=(2, 3),
        placers=None,
    ):
        """
        Initialize a :class:`.DisulfideFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are *bonder* atoms.

        deleters : :class:`tuple` of :class:`int`
            The indices of atoms in the functional group string, which
            are *deleter* atoms.

        placers : :class:`tuple` of :class:`int`, optional
            The indices of atoms in the functional group string, which
            are *placer* atoms. If ``None``, `bonders` will be used.

        """

        self._bonders = bonders
        self._deleters = deleters
        self._placers = bonders if placers is None else placers

    def get_functional_groups(self, molecule):
        ids = _get_atom_ids('[S]([*])[S]([*])', molecule)
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield Disulfide(
                sulfur1=atoms[0],
                atom1=atoms[1],
                sulfur2=atoms[2],
                atom2=atoms[3],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
                placers=tuple(atoms[i] for i in self._placers),
            )
