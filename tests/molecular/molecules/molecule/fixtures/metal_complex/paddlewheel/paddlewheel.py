import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Cu+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_copper_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _copper_atom.get_atoms(0)
_copper_atom = _copper_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(4))
)

_bi_1 = stk.BuildingBlock.init_from_file(
    smiles='O=C(O)c1ccc(Br)cc1',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#8]~[#1]',
            bonders=(1, ),
            deleters=(2, ),
        ),
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#8]',
            bonders=(1, ),
            deleters=(),
        ),
    ]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.metal_complex.PaddleWheel(
                    metals={_copper_atom: (0, 1)},
                    ligands={_bi_1: (0, 1, 2, 3)},
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
                'Brc1ccc(C2=O->[Cu+2]34<-O=C(c5ccc(Br)cc5)O->[Cu+2]'
                '(<-O2)(<-OC(c2ccc(Br)cc2)=O->3)<-OC(c2ccc(Br)cc2)'
                '=O->4)cc1'
            ),
        ),
    ),
)
def metal_complex_paddlewheel(request):
    return request.param
