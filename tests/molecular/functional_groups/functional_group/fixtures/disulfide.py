import stk
import pytest

from ..case_data import GenericCaseData


@pytest.fixture
def disulfide(get_atom_ids):
    a, b, c, d = get_atom_ids(4)
    return _disulfide(stk.S(a), stk.C(b), stk.S(c), stk.C(d))


def _disulfide(sulfur1, atom1, sulfur2, atom2):
    bonders = (sulfur1, )
    deleters = (sulfur2, atom2)
    return GenericCaseData(
        functional_group=stk.Disulfide(
            sulfur1=sulfur1,
            atom1=atom1,
            sulfur2=sulfur2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(sulfur1, atom1, sulfur2, atom2),
        bonders=bonders,
        deleters=deleters,
    )