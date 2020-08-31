"""
Acetal
======

"""

from .generic_functional_group import GenericFunctionalGroup


class Acetal(GenericFunctionalGroup):
    """
    Represents an acetal functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``([atom3])[carbon]([oxygen1][atom1])([oxygen2][atom2])([atom4])``.

    """

    def __init__(
        self,
        carbon,
        oxygen1,
        atom1,
        oxygen2,
        atom2,
        atom3,
        atom4,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Acetal` instance.

        Parameters
        ----------
        carbon : :class:`.C`
            The ``[carbon]`` atom.

        oxygen1 : :class:`.O`
            The ``[oxygen1]`` atom.

        atom1 : :class:`.Atom`
            The ``[atom1]`` atom.

        oxygen2 : :class:`.O`
            The ``[oxygen2]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom2]`` atom.

        atom3 : :class:`.Atom`
            The ``[atom3]`` atom.
        
        atom4 : :class:`.Atom`
            The ``[atom4]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            the bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._carbon = carbon
        self._oxygen1 = oxygen1
        self._atom1 = atom1
        self._oxygen2 = oxygen2
        self._atom2 = atom2
        self._atom3 = atom3
        self._atom4 = atom4
        atoms = (carbon, oxygen1, atom1, oxygen2, atom2, atom3, atom4)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_carbon(self):
        """
        Get the ``[carbon]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon]`` atom.

        """

        return self._carbon

    def get_oxygen1(self):
        """
        Get the ``[oxygen1]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen1]`` atom.

        """

        return self._oxygen1

    def get_atom1(self):
        """
        Get the ``[atom1]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_oxygen2(self):
        """
        Get the ``[oxygen2]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen2]`` atom.

        """

        return self._oxygen2

    def get_atom2(self):
        """
        Get the ``[atom2]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_atom3(self):
        """
        Get the ``[atom3]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom3]`` atom.

        """

        return self._atom3

    def get_atom4(self):
        """
        Get the ``[atom4]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom4]`` atom.

        """

        return self._atom4

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen1 = self._oxygen1
        clone._atom1 = self._atom1
        clone._oxygen2 = self._oxygen2
        clone._atom2 = self._atom2
        clone._atom3 = self._atom3
        clone._atom4 = self._atom4
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._atom3 = atom_map.get(
            self._atom3.get_id(),
            self._atom3,
        )
        clone._atom4 = atom_map.get(
            self._atom4.get_id(),
            self._atom4,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen1}, {self._atom1}, '
            f'{self._oxygen2}, {self._atom2}, {self._atom3}, '
            f'{self._atom4}), '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
