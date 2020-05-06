import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Zn+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_zinc_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _zinc_atom.get_atoms(0)
_zinc_atom = _zinc_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(4))
)

_quad_1 = stk.BuildingBlock.init_from_file(
    smiles=(
        'Brc1ccc(-c2c3nc(c(-c4ccc(Br)cc4)c4ccc([nH]4)c(-c4ccc(Br)cc4)'
        'c4nc(c(-c5ccc(Br)cc5)c5ccc2[nH]5)C=C4)C=C3)cc1'
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
                stk.metal_complex.Porphyrin(
                    metals={_zinc_atom: 0},
                    ligands={_quad_1: 0},
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
                'Brc1ccc(C2=C3C=CC4=N3->[Zn+2]35<-N6=C(C=CC6=C(c6ccc'
                '(Br)cc6)c6ccc2[nH]->36)C(c2ccc(Br)cc2)=c2ccc([nH]->'
                '52)=C4c2ccc(Br)cc2)cc1'
            ),
        ),
    ),
)
def metal_complex_porphyrin(request):
    return request.param
