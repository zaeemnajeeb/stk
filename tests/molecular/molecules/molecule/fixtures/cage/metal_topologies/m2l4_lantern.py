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
                stk.cage.M2L4Lantern(
                    building_blocks={
                        _metal_atom: range(2),
                        linker: range(2, 6)
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
                '[H])=C([H])C(=C5[H])C5=C([H])C([H])=N(->[Pd+2](<-N6'
                '=C([H])C([H])=C2C([H])=C6[H])(<-N2=C([H])C([H])=C(C'
                '([H])=C2[H])C2=C([H])C(=C([H])C([H])=C2[H])C2=C([H]'
                ')C([H])=N->3C([H])=C2[H])<-N2=C([H])C([H])=C(C([H])'
                '=C2[H])C2=C([H])C(=C([H])C([H])=C2[H])C2=C([H])C([H'
                '])=N->4C([H])=C2[H])C([H])=C5[H])C([H])=C1[H]'
            ),
        ),
    ),
)
def mcage_m2l4_lantern(request):
    return request.param
