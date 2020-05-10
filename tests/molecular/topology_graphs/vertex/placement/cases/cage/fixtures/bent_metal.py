import pytest
import numpy as np
import stk
from functools import partial
from scipy.spatial.distance import euclidean
from rdkit.Chem import AllChem as rdkit

from ....case_data import CaseData

vertices = stk.cage.vertices


atom = rdkit.MolFromSmiles('[Pd+2]')
atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))

_metal_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
atom_0, = _metal_atom.get_atoms(0)
_metal_atom = _metal_atom.with_functional_groups(
    (stk.SingleAtom(atom_0) for i in range(4))
)
palladium_bi_1 = stk.BuildingBlock(
    smiles='[H]N([H])C([H])([H])C([H])([H])N([H])[H]',
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts='[#7]~[#6]',
            bonders=(0, ),
            deleters=(),
        ),
    ]
)
palladium_cispbi_sqpl = stk.ConstructedMolecule(
    stk.metal_complex.CisProtectedSquarePlanar(
        metals={_metal_atom: 0},
        ligands={palladium_bi_1: 0},
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


@pytest.fixture
def bent_metal(position, bent_aligner_edge, bent_building_block):

    core_fg_direction = np.array([-1, 0, 0], dtype=np.float64)

    def get_core_fg_direction(building_block):
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        vector = core_centroid - position
        return (vector)/np.linalg.norm(vector)

    vertex = vertices._BentMetalComplexCageVertex(
        id=0,
        position=position,
        aligner_edge=bent_aligner_edge,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_bent_edges(vertex)),
        building_block=bent_building_block,
        position=position,
        alignment_tests={
            get_core_fg_direction: core_fg_direction,
        },
        functional_group_edges=(
            {0: 0, 1: 1} if bent_aligner_edge == 0 else {0: 1, 1: 0}
        ),
    )


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_bent_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [10, 10, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [10, -10, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)


@pytest.fixture(params=(0, 1))
def bent_aligner_edge(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock.init_from_molecule(
            molecule=palladium_cispbi_sqpl,
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[Pd]~[#7]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ]
        ),
    ),
)
def bent_building_block(request):
    return request.param


@pytest.fixture(
    params=(
        [1, 2, -20],
    ),
)
def position(request):
    """
    The `position` of a vertex.

    """

    return np.array(request.param, dtype=np.float64)
