"""Recursive composition, open leaves, budgets and proof-tree replay."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.decomposition import DecompositionError
from alggeom.decomposition.chains import algebraic_chain, rational_chain
from alggeom.decomposition.discovery import search_algebraic
from alggeom.decomposition.rational import RationalDecomposition

x, y, z = sp.symbols("x y z")
PAPER = x**4 - z * x * x * (x + 1) + z**3 * (x + 1) ** 2


def test_rational_chain():
    result = rational_chain(x**8, x)
    assert result.complete and result.status == "complete" and result.verify()
    assert len(result.components) == 3
    assert result.compose() == x**8
    assert all(leaf.stop_reason == "indecomposable" for leaf in result.leaves)
    assert all(leaf.search.complete and leaf.search.verify() for leaf in result.leaves)
    assert all(len(c.variables) == 1 for c in result.components)


def test_supplied_rational_split():
    t = sp.Symbol("t")
    initial = RationalDecomposition(1 / (1 + x**8), (x,), t, sp.QQ, t**4 / (t**4 + 1), 1 / x**2)
    result = rational_chain(initial.original_expression, x, initial=initial)
    assert result.root.supplied is initial and result.root.search is None
    assert len(result.components) == 3 and result.complete and result.verify()
    assert sp.cancel(result.compose() - initial.function) == 0
    assert not replace(result, root=replace(result.root, supplied=None)).verify()
    with pytest.raises(DecompositionError, match="initial certificate"):
        rational_chain(x**8, x, initial=initial)


def test_supplied_implicit_split():
    initial = search_algebraic(PAPER, (x, z)).decompositions[0]
    result = algebraic_chain(PAPER, (x, z), initial=initial, max_depth=1)
    assert result.root.supplied is initial and result.root.search is None and result.verify()
    assert len(result.components) == 2 and (not result.complete)
    with pytest.raises(DecompositionError, match="initial certificate"):
        algebraic_chain(PAPER, (x, z), initial=replace(initial, unit=7))


def test_multivariate_leaf():
    function = (x * x + y * y) ** 4
    result = rational_chain(function, (x, y))
    assert result.complete and result.verify() and (len(result.components) == 3)
    assert result.components[-1].variables == (x, y)
    assert all(len(c.variables) == 1 for c in result.components[:-1])
    assert sp.cancel(result.compose() - function) == 0


@pytest.mark.parametrize("function", [x**6, (x + 1 / x) ** 4, (x**2 + 1) ** 4])
def test_rational_identities(function):
    result = rational_chain(function, x)
    assert result.verify() and result.complete
    assert sp.cancel(result.compose() - function) == 0


@pytest.mark.parametrize("function", [0, 2, x + 1, (x + 1) / (x - 1), x**3 + x, x**4 + x])
def test_unsplit_leaves(function):
    result = rational_chain(function, x)
    assert result.complete and result.verify() and (len(result.components) == 1)
    assert result.root.split is None and result.root.resolved


@pytest.mark.parametrize(
    "options,reason",
    [
        ({"max_depth": 0}, "max_depth"),
        ({"max_depth": 1}, "max_depth"),
        ({"max_nodes": 1}, "max_nodes"),
        ({"max_nodes": 3}, "max_nodes"),
        ({"max_candidates": 0}, "search_incomplete"),
        ({"max_degree": 4}, "degree_limit"),
    ],
)
def test_open_limits(options, reason):
    result = rational_chain(x**8, x, **options)
    assert result.status == "incomplete" and (not result.complete) and result.verify()
    assert any(leaf.stop_reason == reason for leaf in result.leaves)
    assert result.compose() == x**8


def test_resolved_chain():
    result = rational_chain(x**8, x, max_candidates=1)
    assert result.root.search.decompositions and (not result.root.search.complete)
    assert result.complete and result.verify()


def test_number_fields_and_poles():
    K = sp.QQ.algebraic_field(sp.sqrt(2))
    result = rational_chain((x * x + sp.sqrt(2)) ** 4, x, constants=K)
    assert result.constants == K and result.complete and result.verify()
    original = sp.Mul(x**8 - 1, sp.Pow(x - 1, -1, evaluate=False), evaluate=False)
    guarded = rational_chain(original, x, max_depth=0)
    assert guarded.input_guard.subs(x, 1) is sp.false and guarded.verify()


@pytest.mark.parametrize("relation", [x - z**8, x**8 - z])
def test_graph_orientations(relation):
    result = algebraic_chain(relation, (x, z))
    assert result.complete and result.verify() and (len(result.components) == 3)
    assert result.exponent == 1
    assert (result.resultant() - result.relation * result.unit).is_zero
    assert result.root.rational_search.complete
    assert result.root.rational_search.verify()


def test_implicit_negative():
    relation = sp.Poly(x * x + z * x + z, x, z, domain=sp.QQ)
    result = algebraic_chain(relation, (x, z))
    assert result.complete and result.verify() and (len(result.components) == 1)
    assert result.root.search.complete and result.root.stop_reason == "indecomposable"


def test_resultant_power():
    result = algebraic_chain(x**3 - 3 * x - z, (x, z), max_depth=4)
    assert result.root.split.exponent == 2
    assert result.exponent >= 2 and result.verify() and (not result.complete)
    assert any(leaf.stop_reason == "repeated_component" for leaf in result.leaves)
    assert (result.resultant() - result.unit * result.relation**result.exponent).is_zero


def test_implicit_reconstruction():
    result = algebraic_chain(PAPER, (x, z), max_depth=3)
    assert result.root.split.relation == result.relation
    assert result.root.split.verify() and result.verify()
    assert len(result.components) >= 2
    assert all(
        leaf.resolved
        or leaf.stop_reason
        in {"max_depth", "repeated_component", "search_incomplete", "degree_limit", "unsupported"}
        for leaf in result.leaves
    )


@pytest.mark.parametrize(
    "options,reason",
    [
        ({"max_depth": 0}, "max_depth"),
        ({"max_nodes": 1}, "max_nodes"),
        ({"max_candidates": 0}, "search_incomplete"),
        ({"max_rational_degree": 3}, "degree_limit"),
    ],
)
def test_graph_limits(options, reason):
    result = algebraic_chain(x - z**8, (x, z), **options)
    assert not result.complete and result.verify()
    assert result.root.stop_reason == reason


def test_factor_degree_limit():
    result = algebraic_chain(x**7 + z * x + z**2, (x, z))
    assert result.root.stop_reason == "degree_limit" and result.verify() and (not result.complete)
    assert result.resultant() == result.relation


def test_altered_rational_tree():
    result = rational_chain(x**8, x)
    bad = replace(result.root, outer=replace(result.root.outer, function=x**2))
    assert not replace(result, root=bad).verify()
    bad = replace(result.root, split=replace(result.root.split, outer=x**2))
    assert not replace(result, root=bad).verify()
    limited = rational_chain(x**8, x, max_depth=0)
    assert not replace(
        limited, root=replace(limited.root, stop_reason="indecomposable", detail=None)
    ).verify()
    assert not replace(limited, max_depth=2).verify()
    assert not replace(result, max_nodes=1).verify()
    assert not replace(result, root=replace(result.root, inner=None)).verify()
    assert not replace(result, root=replace(result.root, search=None)).verify()


def test_altered_implicit_tree():
    result = algebraic_chain(x - z**4, (x, z))
    assert result.verify()
    assert not replace(result, root=replace(result.root, rational_search=None)).verify()
    bad = replace(result.root, split=replace(result.root.split, unit=7))
    assert not replace(result, root=bad).verify()
    bad = replace(result.root, inner=replace(result.root.inner, function=sp.Poly(y - z * z, y, z)))
    assert not replace(result, root=bad).verify()


@pytest.mark.parametrize(
    "options",
    [
        {"max_depth": -1},
        {"max_depth": True},
        {"max_nodes": 0},
        {"max_nodes": 1.5},
        {"max_candidates": -1},
        {"max_degree": 0},
    ],
)
@pytest.mark.parametrize("mode", ["rational", "algebraic"])
def test_invalid_budgets(options, mode):
    with pytest.raises(DecompositionError):
        if mode == "rational":
            rational_chain(x**8, x, **options)
        else:
            algebraic_chain(x - z**8, (x, z), **options)


@pytest.mark.parametrize(
    "function,variables", [(sp.sin(x), x), (x + y, x), (1.2 * x, x), (x, ()), (x, (x, x))]
)
def test_invalid_inputs(function, variables):
    with pytest.raises(DecompositionError):
        rational_chain(function, variables)


@pytest.mark.parametrize(
    "relation,variables", [(x * x - z * z, (x, z)), (x * x - z, (x,)), (x * x - z, (x, x))]
)
def test_invalid_relations(relation, variables):
    with pytest.raises(DecompositionError):
        algebraic_chain(relation, variables)
