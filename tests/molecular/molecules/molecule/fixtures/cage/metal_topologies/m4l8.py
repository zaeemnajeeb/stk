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
                stk.cage.M4L8(
                    building_blocks={
                        _metal_atom: range(4),
                        linker: range(4, 12)
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
                '[H]C1=C([H])C2=C([H])C(=C1[H])C1=C([H])C([H])=N(->['
                'Pd+2]34<-N5=C([H])C([H])=C(C([H])=C5[H])C5=C([H])C('
                '[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2]6(<-N'
                '7=C([H])C([H])=C(C([H])=C7[H])C7=C([H])C([H])=C([H]'
                ')C(=C7[H])C7=C([H])C([H])=N(->[Pd+2]8(<-N9=C([H])C('
                '[H])=C(C([H])=C9[H])C9=C([H])C([H])=C([H])C(=C9[H])'
                'C9=C([H])C([H])=N(->[Pd+2](<-N%10=C([H])C([H])=C(C('
                '[H])=C%10[H])C%10=C([H])C(=C([H])C([H])=C%10[H])C%1'
                '0=C([H])C([H])=N->3C([H])=C%10[H])(<-N3=C([H])C([H]'
                ')=C(C([H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])C3='
                'C([H])C([H])=N->4C([H])=C3[H])<-N3=C([H])C([H])=C(C'
                '([H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])C3=C([H]'
                ')C([H])=N->8C([H])=C3[H])C([H])=C9[H])<-N3=C([H])C('
                '[H])=C(C([H])=C3[H])C3=C([H])C(=C([H])C([H])=C3[H])'
                'C3=C([H])C([H])=N->6C([H])=C3[H])C([H])=C7[H])<-N3='
                'C([H])C([H])=C2C([H])=C3[H])C([H])=C5[H])C([H])=C1[H]'
            ),
        ),
    ),
)
def mcage_m4l8(request):
    return request.param
