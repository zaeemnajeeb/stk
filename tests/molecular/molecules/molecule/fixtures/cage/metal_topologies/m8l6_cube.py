import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Pd+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_metal_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _metal_atom.get_atoms(0)
_metal_atom = _metal_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(4))
)
linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=C([H])C(C3=C([H])C([H]'
        ')=NC([H])=C3[H])=C2[H])=C1[H]'
    ),
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#6]',
            bonders=(1, ),
            deleters=(),
        ),
    ]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M24L48(
                    building_blocks={
                        _metal_atom: range(6),
                        linker: range(6, 18)
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            ),
            smiles=(
                '[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C'
                '([H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3'
                '[H])C3C7=C([H])C([H])=C(C([H])=C7[H])N7->[Fe+2]89%10'
                '(<-N%11=C([H])C([H])=C([H])C([H])=C%11C([H])=N->8C8=C'
                '([H])C([H])=C(C([H])=C8[H])C(C8=C([H])C([H])=C(C([H]'
                ')=C8[H])N->4=C([H])C2=C1[H])C1C2=C([H])C([H])=C(C([H]'
                ')=C2[H])N2->[Fe+2]48%11(<-N%12=C([H])C([H])=C([H])C(['
                'H])=C%12C=2[H])<-N2=C([H])C([H])=C([H])C([H])=C2C([H]'
                ')=N->4C2=C([H])C([H])=C(C([H])=C2[H])C2C4=C([H])C([H]'
                ')=C(C([H])=C4[H])N4->[Fe+2]%12%13(<-N%14=C([H])C([H]'
                ')=C([H])C([H])=C%14C([H])=N->%12C%12=C([H])C([H])=C1C'
                '([H])=C%12[H])(<-N1=C([H])C([H])=C([H])C([H])=C1C=4[H'
                '])<-N1=C([H])C([H])=C([H])C([H])=C1C([H])=N->%13C1=C('
                '[H])C([H])=C(C([H])=C1[H])C(C1=C([H])C([H])=C(C([H])='
                'C1[H])N->9=C([H])C1=C([H])C([H])=C([H])C([H])=N->%101'
                ')C1C4=C([H])C([H])=C(C([H])=C4[H])N4->[Fe+2]9%10%12('
                '<-N%13=C([H])C([H])=C([H])C([H])=C%13C([H])=N->9C9=C('
                '[H])C([H])=C(C([H])=C9[H])C9C%13=C([H])C([H])=C(C([H]'
                ')=C%13[H])N%13->[Fe+2]%14%15(<-N%16=C([H])C([H])=C([H'
                '])C([H])=C%16C=%13[H])(<-N%13=C([H])C([H])=C([H])C([H'
                '])=C%13C([H])=N->%14C%13=C([H])C([H])=C(C([H])=C%13[H'
                '])C3C3=C([H])C([H])=C(C([H])=C3[H])N->%10=C([H])C3=C('
                '[H])C([H])=C([H])C([H])=N->%123)<-N3=C([H])C([H])=C(['
                'H])C([H])=C3C([H])=N->%15C3=C([H])C([H])=C(C([H])=C3['
                'H])C(C3=C([H])C([H])=C(C([H])=C3[H])N3->[Fe+2]%10%12%'
                '13(<-N%14=C([H])C([H])=C([H])C([H])=C%14C([H])=N->%10'
                'C%10=C([H])C([H])=C(C([H])=C%10[H])C9C9=C([H])C([H])='
                'C(C([H])=C9[H])N9->[Fe+2]%10%14(<-N%15=C([H])C([H])=C'
                '([H])C([H])=C%15C=9[H])(<-N9=C([H])C([H])=C([H])C([H]'
                ')=C9C([H])=N->%10C9=C([H])C([H])=C(C([H])=C9[H])C2C2='
                'C([H])C([H])=C(C([H])=C2[H])N->%12=C([H])C2=C([H])C(['
                'H])=C([H])C([H])=N->%132)<-N2=C([H])C([H])=C([H])C([H'
                '])=C2C([H])=N->%14C2=C([H])C([H])=C1C([H])=C2[H])<-N1'
                '=C([H])C([H])=C([H])C([H])=C1C=3[H])C(C1=C([H])C([H])'
                '=C(C([H])=C1[H])N->5=C([H])C1=C([H])C([H])=C([H])C([H'
                '])=N->61)C1=C([H])C([H])=C(C([H])=C1[H])N->8=C([H])C1'
                '=C([H])C([H])=C([H])C([H])=N->%111)<-N1=C([H])C([H])='
                'C([H])C([H])=C1C=4[H])<-N1=C([H])C([H])=C([H])C([H])='
                'C1C=7[H]'
            ),
        ),
    ),
)
def cage_six_plus_eight(request):
    return request.param
