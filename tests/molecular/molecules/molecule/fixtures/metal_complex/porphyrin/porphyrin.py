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

_quad_1 = stk.BuildingBlock(
    smiles=(
        'Brc1ccc(C2=C3C=CC(=C(c4ccc(Br)cc4)C4=NC(=C(c5ccc(Br)cc5)C5=C'
        'C=C([N]5)C(c5ccc(Br)cc5)=C5C=CC2=N5)C=C4)[N]3)cc1'
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
                'Brc1ccc(C2=C3C=CC4=C(c5ccc(Br)cc5)C5=N6->[Zn+2]7('
                '<-[N]8C(=CC=C8C(c8ccc(Br)cc8)=C6C=C5)C(c5ccc(Br)cc5)'
                '=C5C=CC2=N->75)<-[N]34)cc1'
            ),
        ),
    ),
)
def metal_complex_porphyrin(request):
    return request.param
