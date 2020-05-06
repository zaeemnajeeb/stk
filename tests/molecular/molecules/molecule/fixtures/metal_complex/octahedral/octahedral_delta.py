import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Fe+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_iron_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _iron_atom.get_atoms(0)
_iron_atom = _iron_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(6))
)

_iron_bi_1 = stk.BuildingBlock(
    smiles='BrN=Cc1ccccn1',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#35]',
            bonders=(1, ),
            deleters=(),
        ),
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
                stk.metal_complex.OctahedralDelta(
                    metals={_iron_atom: 0},
                    ligands={_iron_bi_1: (0, 1, 2)},
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
                '[H]C1=N(Br)->[Fe+2]23(<-N(Br)=C([H])c4c([H])c([H])c'
                '([H])c([H])n->24)(<-N(Br)=C([H])c2c([H])c([H])c([H])'
                'c([H])n->32)<-n2c([H])c([H])c([H])c([H])c21'
            ),
        ),
    ),
)
def metal_complex_octahedraldelta(request):
    return request.param
