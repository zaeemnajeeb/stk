"""
Disulfide Exchange Factory
==========================

"""

from .functional_group_factory import FunctionalGroupFactory
from .utilities import _get_atom_ids
from .utilities import _find_R_ids
from ..functional_groups import Disulfide


class DisulfideExchangeFactory(FunctionalGroupFactory):
    """
    Creates :class:`.Disulfide` instances.

    Creates functional groups from substructures, which match the
    ``[S]([*])[S]([*])`` functional group string.

    Examples
    --------
    You want to create a building block which has
    :class:`.Disulfide` functional groups, and you want this to
    behave like a disulfide in a exchange reaction.
    You want the sulfur atom in the functional group
    connected to the most atoms (largest R groups) to be the 
    *bonder* atom, and the other sulfur atom to be the 
    *deleter* atom along with any other atoms 
    connected to the *deleter*.

    .. code-block:: python

        import stk

        building_block = stk.BuildingBlock(
            smiles='COCCCOCCSSCCC',
            functional_groups=(stk.DisulfideExchangeFactory(), ),
        )

    You want to create a building block which has
    :class:`.Disulfide` functional groups, and you want this to
    behave like a disulfide in a metathesis reaction. 
    There are multiple disulfides and you want them all to 
    be correctly represented in the reaction. This situation
    can occur should there be more than 1 disulfide present.

    .. code-block:: python

        import stk

        disulfide_exchange_factory = stk.DisulfideExchangeFactory(
            # The indices of the sulfur atoms in the functional
            # group string (see docstring) are 0 and 2.
            # As there are 2 disulfides, specify the *bonder* atoms
            # in the order that the disulfides appear according to 
            # the atom_ids
            bonders=(2, 0),
            # As disulfides would be incorectly assigned based on
            # the default option, set largest_chain to false
            largest_chain = False,
        )
        building_block = stk.BuildingBlock(
            smiles='C(CC)CCCSSCCCSSC',
            functional_groups=(disulfide_exchange_factory, ),
        )

    See Also
    --------
    :class:`.GenericFunctionalGroup`
        Defines *bonders* and  *deleters*.

    """

    def __init__(
        self,
        bonders = None,
        deleters = None,
        placers = None,
        largest_chain = True
    ):
        """
        Initialize a :class:`.DisulfideFactory` instance.

        Parameters
        ----------
        bonders : :class:`tuple` of :class:`int`, optional
            The indices of atoms in the functional group string, which
            are *bonder* atoms. This should be either 0 or 2 for this
            reaction. Optional if `largest_chain` is True, otherwise 
            must be stated as a tuple of integers. Each value in 
            `bonders` refers to one alkene in the molecule. If there are
            three alkenes in the molecule, bonders will contain three 
            integers of either 0 or 3 e.g. bonders = (0, 0, 2)

        deleters : :class:`tuple` of :class:`int`, optional
            The indices of atoms in the functional group string, which
            are *deleter* atoms. Automatically assigned according
            to the *bonder* atoms.

        placers : :class:`tuple` of :class:`int`, optional
            The indices of atoms in the functional group string, which
            are *placer* atoms. If ``None``, `bonders` will be used.

        largest_chain : :bool:, optional
            Toggles whether the factory will assign `bonders` and 
            `deleters` for a given alkene based on the number of 
            atoms attached to a given carbon in the C=C functional
            group.  
        
        Returns
        -------
        :class:`.Disulfide`
            If `largest_chain` is ``True``, assigns the `bonders` 
            atom which results in the least number of *deleter*
            atoms. An example would be where `SC` is 
            assigned as the *deleter* atoms for 'CCCSSC'. This 
            refers to the `-SC` as it has the fewest atoms to
            be deleted relative to `CCCS-`. This count includes
            Hydrogens.


        :class:`.Disulfide`
            If `largest_chain` is ``False``, the *bonder* atoms
            must be explicitly stated via `bonders`.

        """

        if largest_chain is False:
            deleters = []
            if bonders is None:
                raise ValueError(
            'Bonders must either be (0,) or (2,) if largest_chain is False'
            )

            for i in bonders:
                if int(i) == 0: 
                    deleters.append(2)
                elif int(i) == 2:
                    deleters.append(0)
                else:
                    raise ValueError(
                'Bonders should be (0,) or (2,) if largest_chain is False'
                ) 
            deleters=tuple(deleters)


        self._bonders = bonders
        self._deleters = deleters
        self._largest_chain = largest_chain
        if self._bonders is not None:
            self._placers = bonders if placers is None else placers

    def get_functional_groups(self, molecule):
        functional_group_number = 0
        ids = _get_atom_ids('[S]([*])[S]([*])', molecule)
        for atom_ids in ids:
            atoms = tuple(molecule.get_atoms(atom_ids))

            if self._largest_chain is True:
                R_to_delete_ids_1 = _find_R_ids(
                molecule = molecule,
                deleters = set([atoms[0].get_id()]),
                bonder = atoms[2].get_id()
                )

                R_to_delete_ids_2 = _find_R_ids(
                molecule = molecule,
                deleters = set([atoms[2].get_id()]),
                bonder = atoms[0].get_id()
                )

                if len(R_to_delete_ids_1) < len(R_to_delete_ids_2):
                    R_deleters = tuple(molecule.get_atoms(list(R_to_delete_ids_1)))
                    _placers = _bonders = tuple((atoms[2], ))

                elif len(R_to_delete_ids_1) > len(R_to_delete_ids_2):
                    R_deleters = tuple(molecule.get_atoms(list(R_to_delete_ids_2)))
                    _placers = _bonders = tuple((atoms[0], ))

                else:
                    #if the R groups are identical in length, default to first
                    R_deleters = tuple(molecule.get_atoms(list(R_to_delete_ids_1)))
                    _placers = _bonders = tuple((atoms[2], ))

            if self._largest_chain is False:
                bonder_in_fg = self._bonders[functional_group_number]
                _bonders = tuple((atoms[bonder_in_fg], ))
                
                _placers = _bonders if self._placers is self._bonders else tuple(atoms[i] for i in self._placers)
                
                deleter_in_fg = self._deleters[functional_group_number]
                _deleters = tuple((atoms[deleter_in_fg], ))
                
                R_to_delete_ids = _find_R_ids(
                    molecule = molecule,
                    deleters = set([_deleters[0].get_id()]),
                    bonder = _bonders[0].get_id()
                )
                R_deleters = tuple(molecule.get_atoms(list(R_to_delete_ids)))

            yield Disulfide(
                sulfur1=atoms[0],
                atom1=atoms[1],
                sulfur2=atoms[2],
                atom2=atoms[3],
                bonders=_bonders,
                deleters=R_deleters,
                placers=_placers,
            )
            functional_group_number += 1
