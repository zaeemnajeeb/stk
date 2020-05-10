import pytest
import stk
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData


atom = rdkit.MolFromSmiles('[Fe+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_metal_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _metal_atom.get_atoms(0)
_metal_atom = _metal_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(6))
)

tritopic_linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=C([H])C(N(C2=C([H])C([H])=C(Br)C([H])=C2[H])C2=C([H])C('
        '[H])=C(Br)C([H])=C2[H])=C([H])C([H])=C1Br'
    ),
    functional_groups=[stk.BromoFactory()]
)
tetratopic_linker = stk.BuildingBlock(
    smiles=(
        '[H]C1=C([H])C(C(C2=C([H])C([H])=C(Br)C([H])=C2[H])C(C2=C([H])'
        'C([H])=C(Br)C([H])=C2[H])C2=C([H])C([H])=C(Br)C([H])=C2[H])=C'
        '([H])C([H])=C1Br'
    ),
    functional_groups=[stk.BromoFactory()]
)
complex_ligand = stk.BuildingBlock(
    smiles='[H]C1=NC(C([H])=NBr)=C([H])C([H])=C1[H]',
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
iron_complex = stk.ConstructedMolecule(
    stk.metal_complex.OctahedralDelta(
        metals={_metal_atom: 0},
        ligands={complex_ligand: (0, 1, 2)},
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
)
iron_complex = stk.BuildingBlock.init_from_molecule(
    molecule=iron_complex,
    functional_groups=[stk.BromoFactory()]
)


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.ConstructedMolecule(
                stk.cage.M4L4Tetrahedron(
                    building_blocks={
                        iron_complex: range(6),
                        tritopic_linker: range(6, 8),
                        tetratopic_linker: range(8, 11),
                    },
                )
            ),
            smiles=(
                '[H]C1=C([H])C([H])=N2->[Fe+2]3456<-N7=C([H])C([H])=C'
                '([H])C([H])=C7C([H])=N->3C3=C([H])C([H])=C(C([H])=C3'
                '[H])C3C7=C([H])C([H])=C(C([H])=C7[H])N7->[Fe+2]89%10'
                '(<-N%11=C([H])C([H])=C([H])C([H])=C%11C=7[H])<-N7=C('
                '[H])C([H])=C([H])C([H])=C7C([H])=N->8C7=C([H])C([H])'
                '=C(C([H])=C7[H])C7C8=C([H])C([H])=C(C([H])=C8[H])N8-'
                '>[Fe+2]%11%12%13(<-N%14=C([H])C([H])=C([H])C([H])=C%'
                '14C([H])=N->%11C%11=C([H])C([H])=C(C([H])=C%11[H])C3'
                'C3=C([H])C([H])=C(C([H])=C3[H])N3->[Fe+2]%11%14%15(<'
                '-N%16=C([H])C([H])=C([H])C([H])=C%16C([H])=N->%11C%1'
                '1=C([H])C([H])=C(C([H])=C%11[H])C%11C%16=C([H])C([H]'
                ')=C(C([H])=C%16[H])N%16->[Fe+2]%17%18(<-N%19=C([H])C'
                '([H])=C([H])C([H])=C%19C=%16[H])(<-N%16=C([H])C([H])'
                '=C([H])C([H])=C%16C([H])=N->%17C%16=C([H])C([H])=C(C'
                '([H])=C%16[H])C7C7=C([H])C([H])=C(C([H])=C7[H])N7->['
                'Fe+2]%16%17(<-N%19=C([H])C([H])=C([H])C([H])=C%19C(['
                'H])=N->%16C%16=C([H])C([H])=C(C([H])=C%16[H])C%11C%1'
                '1=C([H])C([H])=C(C([H])=C%11[H])N->4=C([H])C2=C1[H])'
                '(<-N1=C([H])C([H])=C([H])C([H])=C1C=7[H])<-N1=C([H])'
                'C([H])=C([H])C([H])=C1C([H])=N->%17C1=C([H])C([H])=C'
                '(C([H])=C1[H])N(C1=C([H])C([H])=C(C([H])=C1[H])N->5='
                'C([H])C1=C([H])C([H])=C([H])C([H])=N->61)C1=C([H])C('
                '[H])=C(C([H])=C1[H])N->9=C([H])C1=C([H])C([H])=C([H]'
                ')C([H])=N->%101)<-N1=C([H])C([H])=C([H])C([H])=C1C(['
                'H])=N->%18C1=C([H])C([H])=C(C([H])=C1[H])N(C1=C([H])'
                'C([H])=C(C([H])=C1[H])N->%14=C([H])C1=C([H])C([H])=C'
                '([H])C([H])=N->%151)C1=C([H])C([H])=C(C([H])=C1[H])N'
                '->%12=C([H])C1=C([H])C([H])=C([H])C([H])=N->%131)<-N'
                '1=C([H])C([H])=C([H])C([H])=C1C=3[H])<-N1=C([H])C([H'
                '])=C([H])C([H])=C1C=8[H]'
            ),
        ),
    ),
)
def mcage_m6l2l3_prism(request):
    return request.param
