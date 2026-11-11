"""Mathematical root results and specialist certificate boundaries."""

import pytest
import sympy as sp

import alggeom
from alggeom import decomposition
from alggeom.decomposition import (
    IncompleteDecompositionError,
    decompose,
    decompose_relation,
    decompositions,
)
from alggeom.decomposition.rational import search_rational

x, y, t, a = sp.symbols("x y t a", real=True)


def compose(maps, *inputs):
    value = maps[-1](*inputs)
    for component in reversed(maps[:-1]):
        value = component(value)
    return sp.cancel(value)


@pytest.mark.parametrize("variables", [x, (x,), [x]])
def test_lambda_chain(variables):
    maps = decompose(x**8, variables)
    assert len(maps) == 3
    assert all(isinstance(f, sp.Lambda) for f in maps)
    assert compose(maps, x) == x**8
    assert [sp.degree(f.expr, f.variables[0]) for f in maps] == [2, 2, 2]


def test_multivariate_lambdas():
    maps = decompose((x * x + y * y) ** 4, (v for v in (x, y)))
    assert len(maps) == 3
    assert maps[-1].variables == (x, y)
    assert sp.cancel(compose(maps, x, y) - (x * x + y * y) ** 4) == 0


@pytest.mark.parametrize("expression", [0, 3, x + 1, x**4 + x])
def test_single_lambda(expression):
    maps = decompose(expression, x)
    assert len(maps) == 1
    assert maps[0](x) == expression
    assert decompositions(expression, x) == ()


def test_all_pairs():
    pairs = decompositions(x**6, x, parameter=t)
    assert len(pairs) == 2
    assert {tuple(sp.degree(f.expr, f.variables[0]) for f in pair) for pair in pairs} == {
        (2, 3),
        (3, 2),
    }
    assert all(sp.cancel(outer(inner(x)) - x**6) == 0 for outer, inner in pairs)


@pytest.mark.parametrize(
    "call",
    [
        lambda: decompose(x**8, x, max_depth=1),
        lambda: decompositions(x**6, x, max_candidates=1),
        lambda: decompose_relation(y**8 - x, (y, x), max_depth=1),
    ],
)
def test_incomplete_results(call):
    with pytest.raises(IncompleteDecompositionError) as error:
        call()
    assert not error.value.result.complete
    assert error.value.result.verify()
    assert (
        len(
            error.value.result.components
            if hasattr(error.value.result, "components")
            else error.value.result.decompositions
        )
        >= 1
    )


def test_implicit_polynomials():
    components = decompose_relation(y**8 - x, (y, x))
    assert len(components) == 3
    result = components[-1]
    for outer in reversed(components[:-1]):
        result = sp.Poly(
            sp.resultant(outer.as_expr(), result.as_expr(), outer.gens[1]),
            outer.gens[0],
            result.gens[1],
            domain=sp.QQ,
        )
    assert result.monic() == sp.Poly(y**8 - x, y, x, domain=sp.QQ)


def test_export_surface():
    assert set(alggeom.__all__) == {
        "absolute_factorization",
        "groebner_cover",
        "groebner_certificates",
        "function_field_constants",
        "parameter_fields",
        "parameter_localization",
        "parameter_cover",
        "chevalley",
        "comprehensive_strata",
        "constructible",
        "families",
        "geometric",
        "linear_strata",
        "localized_linear_cover",
        "morphism_invariants",
        "nonlinear_strata",
        "projective_normalization",
        "fields",
        "ideals",
        "rings",
        "rational_functions",
        "schemes",
        "maps",
        "points",
        "properties",
        "local",
        "components",
        "noether",
        "normalization",
        "curves",
        "graded",
        "modules",
        "resolutions",
        "decomposition",
    }

    assert alggeom.decomposition is decomposition
    assert decompose is decomposition.decompose
    assert set(decomposition.__all__) == {
        "decompose",
        "decompositions",
        "decompose_relation",
        "decomposition_condition",
        "real_equivalent",
        "common_right_component",
        "DecompositionError",
        "DegreeLimitError",
        "IncompleteDecompositionError",
    }


def test_symbol_and_sequence_agree():
    direct = search_rational((x + 1 / x) ** 2, x, parameter=t)
    sequence = search_rational((x + 1 / x) ** 2, (x,), parameter=t)
    assert direct.decompositions == sequence.decompositions
    assert direct.candidate_count == sequence.candidate_count
    assert direct.complete and len(direct.decompositions) == 3
