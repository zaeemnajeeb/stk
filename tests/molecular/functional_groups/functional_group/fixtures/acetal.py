import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def acetal(get_atom_ids):
    a, b, c, d, e, f, g = get_atom_ids(7)
    return _acetal(
        carbon=stk.C(a),
        oxygen1=stk.O(b),
        atom1=stk.C(c),
        oxygen2=stk.O(d),
        atom2=stk.C(e),
        atom3=stk.C(f),
        atom4=stk.C(g)
    )


def _acetal(carbon, oxygen1, atom1, oxygen2, atom2, atom3, atom4):
    bonders = (carbon)
    deleters = (oxygen1, atom1, oxygen2, atom2)
    return GenericCaseData(
        functional_group=stk.Acetal(
            carbon=carbon,
            oxygen1=oxygen1,
            atom1=atom1,
            oxygen2=oxygen2,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen1, atom1, oxygen2, atom2, atom3, atom4),
        bonders=bonders,
        deleters=deleters,
    )
