import pytest
import stk

from ..case_data import GenericCaseData


@pytest.fixture
def thioacetal(get_atom_ids):
    a, b, c, d, e, f, g = get_atom_ids(7)
    return _thioacetal(
        carbon=stk.C(a),
        sulfur1=stk.S(b),
        atom1=stk.C(c),
        sulfur2=stk.S(d),
        atom2=stk.C(e),
        atom3=stk.C(f),
        atom4=stk.C(g)
    )


def _thioacetal(carbon, sulfur1, atom1, sulfur2, atom2, atom3, atom4):
    bonders = (carbon, )
    deleters = (sulfur1, atom1, sulfur2, atom2)
    return GenericCaseData(
        functional_group=stk.Thioacetal(
            carbon=carbon,
            sulfur1=sulfur1,
            atom1=atom1,
            sulfur2=sulfur2,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, sulfur1, atom1, sulfur2, atom2, atom3, atom4),
        bonders=bonders,
        deleters=deleters,
    )
