import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Pd+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_palladium_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _palladium_atom.get_atoms(0)
_palladium_atom = _palladium_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(4))
)

_mo_1 = stk.BuildingBlock(
    smiles='c1cc2c(cn1)CCCCC2',
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
                stk.metal_complex.SquarePlanar(
                    metals={_palladium_atom: 0},
                    ligands={_mo_1: (0, 1, 2, 3)},
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
                '[H]c1c2c(c([H])n(->[Pd+2](<-n3c([H])c([H])c4c(c3[H])'
                'C([H])([H])C([H])([H])C([H])([H])C([H])([H])C4([H])['
                'H])(<-n3c([H])c([H])c4c(c3[H])C([H])([H])C([H])([H])'
                'C([H])([H])C([H])([H])C4([H])[H])<-n3c([H])c([H])c4c'
                '(c3[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C'
                '4([H])[H])c1[H])C([H])([H])C([H])([H])C([H])([H])C(['
                'H])([H])C2([H])[H]'
            ),
        ),
    ),
)
def metal_complex_squareplanar(request):
    return request.param
