"""Reject systematic changes to identities, embeddings, coverage and tree evidence."""

from dataclasses import replace
from functools import cache

import pytest
import sympy as sp
from algroots.quotient import QuotientAlgebra

from alggeom.decomposition.algebraic import reconstruct_algebraic
from alggeom.decomposition.chains import algebraic_chain, rational_chain
from alggeom.decomposition.discovery import search_algebraic
from alggeom.decomposition.factors import discover_factors
from alggeom.decomposition.rational import search_rational

x, z, a, t = sp.symbols("x z a t")
RELATION = x**4 - z * x**2 * (x + 1) + z**3 * (x + 1) ** 2


@pytest.fixture(scope="module")
def evidence():
    @cache
    def get(kind):
        if kind == "rational_search":
            result = search_rational(x**6, x, parameter=t)
        elif kind == "rational_pair":
            result = get("rational_search").decompositions[0]
        elif kind == "discovery":
            result = discover_factors(RELATION, (x, z), degrees=(2,))
        elif kind == "factor_degree":
            result = get("discovery").levels[0]
        elif kind == "factor_orbit":
            result = get("factor_degree").orbits[0]
        elif kind == "reconstruction":
            ambient = QuotientAlgebra.from_polynomials(
                [a**2 - z * a + z**3], (a,), domain=sp.QQ.frac_field(z)
            )
            result = reconstruct_algebraic(
                RELATION, (x, z), candidate=x**2 - a * (x + 1), ambient=ambient, parameter=t
            )
        elif kind == "algebraic_pair":
            result = get("reconstruction").decomposition
        elif kind == "algebraic_search":
            result = search_algebraic(x**3 - 3 * x - z, (x, z), parameter=t, exhaustive=True)
        elif kind == "rational_chain":
            result = rational_chain(x**8, x)
        elif kind == "algebraic_chain":
            result = algebraic_chain(x**8 - z, (x, z))
        else:
            raise ValueError(kind)
        assert result.verify(), kind
        return result

    return get


MUTATIONS = {
    "rational_pair": {
        "input": lambda c: replace(c, original_expression=c.original_expression + 1),
        "variables": lambda c: replace(c, variables=(z,)),
        "parameter": lambda c: replace(c, parameter=x),
        "constants": lambda c: replace(c, constants=sp.GF(5)),
        "outer": lambda c: replace(c, outer=c.outer + 1),
        "inner": lambda c: replace(c, inner=c.inner + 1),
    },
    "rational_search": {
        "input": lambda c: replace(c, original_expression=c.original_expression + 1),
        "coordinates": lambda c: replace(c, comparison_variables=(x,)),
        "missing_factor": lambda c: replace(c, factors=c.factors[1:]),
        "duplicate_factor": lambda c: replace(c, factors=(*c.factors, c.factors[0])),
        "factor_unit": lambda c: replace(c, factors=(2 * c.factors[0], *c.factors[1:])),
        "count": lambda c: replace(c, candidate_count=c.candidate_count + 1),
        "examined": lambda c: replace(c, examined=0),
        "boolean_count": lambda c: replace(c, candidate_count=True),
        "omitted_answer": lambda c: replace(c, decompositions=c.decompositions[1:]),
        "duplicated_answer": lambda c: replace(
            c, decompositions=(*c.decompositions, c.decompositions[0])
        ),
    },
    "factor_orbit": {
        "relation": lambda c: replace(c, relation=c.relation + 1),
        "degree": lambda c: replace(c, degree=c.degree + 1),
        "coefficient": lambda c: replace(
            c, coefficients=(c.coefficients[0] + 1, *c.coefficients[1:])
        ),
        "coefficient_count": lambda c: replace(c, coefficients=c.coefficients[1:]),
        "cofactor": lambda c: replace(c, cofactor=c.cofactor + 1),
    },
    "factor_degree": {
        "degree": lambda c: replace(c, degree=1),
        "parameter": lambda c: replace(c, parameter=x),
        "weights": lambda c: replace(c, weights=(0,) * c.degree),
        "defining_polynomial": lambda c: replace(c, defining_polynomial=c.defining_polynomial + 1),
        "coordinate": lambda c: replace(
            c,
            coordinate_polynomials=(c.coordinate_polynomials[0] + 1, *c.coordinate_polynomials[1:]),
        ),
        "coordinate_count": lambda c: replace(
            c, coordinate_polynomials=c.coordinate_polynomials[1:]
        ),
        "component": lambda c: replace(
            c, component_polynomials=(c.component_polynomials[0] + 1, *c.component_polynomials[1:])
        ),
        "component_order": lambda c: replace(
            c, component_polynomials=tuple(reversed(c.component_polynomials))
        ),
        "orbit": lambda c: replace(
            c, orbits=(replace(c.orbits[0], cofactor=c.orbits[0].cofactor + 1), *c.orbits[1:])
        ),
        "orbit_count": lambda c: replace(c, orbits=c.orbits[1:]),
        "omitted_component": lambda c: replace(
            c, orbits=c.orbits[1:], component_polynomials=c.component_polynomials[1:]
        ),
    },
    "discovery": {
        "requested_degree": lambda c: replace(c, requested_degrees=(2, 3)),
        "missing_level": lambda c: replace(c, levels=()),
        "duplicate_level": lambda c: replace(c, levels=(*c.levels, c.levels[0])),
        "dimension_budget": lambda c: replace(c, max_factor_dimension=5),
        "separator_budget": lambda c: replace(c, max_separating_attempts=0),
    },
    "reconstruction": {
        "candidate": lambda c: replace(c, candidate=c.candidate + 1),
        "parameter": lambda c: replace(c, parameter=z),
        "coefficient_generators": lambda c: replace(
            c, coefficient_field=replace(c.coefficient_field, generators=(a, a))
        ),
        "rationality": lambda c: replace(
            c, rationality=replace(c.rationality, status="positive_genus")
        ),
        "parametrization": lambda c: replace(c, parametrization=None),
        "missing_answer": lambda c: replace(c, decomposition=None),
        "rejection": lambda c: replace(c, reason="positive_genus"),
        "nested_unit": lambda c: replace(c, decomposition=replace(c.decomposition, unit=2)),
    },
    "algebraic_pair": {
        "relation": lambda c: replace(c, relation=c.relation + 1),
        "outer": lambda c: replace(c, outer=c.outer + 1),
        "inner": lambda c: replace(c, inner=c.inner + 1),
        "power": lambda c: replace(c, exponent=c.exponent + 1),
        "zero_power": lambda c: replace(c, exponent=0),
        "unit": lambda c: replace(c, unit=2),
        "variable_unit": lambda c: replace(c, unit=z),
        "inexact_unit": lambda c: replace(c, unit=sp.Float(1)),
        "coordinate_order": lambda c: replace(
            c, inner=sp.Poly(c.inner.as_expr(), z, t, domain=sp.QQ)
        ),
    },
    "algebraic_search": {
        "parameter": lambda c: replace(c, parameter=x),
        "attempt_order": lambda c: replace(c, attempts=tuple(reversed(c.attempts))),
        "attempt_index": lambda c: replace(
            c, attempts=(replace(c.attempts[0], index=1), *c.attempts[1:])
        ),
        "duplicate_attempt": lambda c: replace(c, attempts=(*c.attempts, c.attempts[0])),
        "missing_diagnostic": lambda c: replace(
            c, attempts=(replace(c.attempts[0], result=None), *c.attempts[1:])
        ),
        "answer_and_error": lambda c: replace(
            c, attempts=(replace(c.attempts[0], error="unresolved"), *c.attempts[1:])
        ),
    },
}


@pytest.mark.parametrize(
    "kind,mutate",
    [
        pytest.param(kind, mutate, id=f"{kind}-{name}")
        for kind, changes in MUTATIONS.items()
        for name, mutate in changes.items()
    ],
)
def test_certificate_changes(evidence, kind, mutate):
    assert mutate(evidence(kind)).verify() is False


TREE_MUTATIONS = {
    "missing_search": lambda c: replace(c, root=replace(c.root, search=None, rational_search=None)),
    "missing_child": lambda c: replace(c, root=replace(c.root, inner=None)),
    "swapped_children": lambda c: replace(
        c, root=replace(c.root, outer=c.root.inner, inner=c.root.outer)
    ),
    "child_function": lambda c: replace(
        c, root=replace(c.root, outer=replace(c.root.outer, function=c.root.outer.function + 1))
    ),
    "child_constants": lambda c: replace(
        c,
        root=replace(
            c.root, inner=replace(c.root.inner, constants=sp.QQ.algebraic_field(sp.sqrt(2)))
        ),
    ),
    "split_and_stop": lambda c: replace(c, root=replace(c.root, stop_reason="indecomposable")),
    "split_and_detail": lambda c: replace(c, root=replace(c.root, detail="unresolved")),
    "node_budget": lambda c: replace(c, max_nodes=1),
    "depth_budget": lambda c: replace(c, max_depth=0),
}


@pytest.mark.parametrize("kind", ["rational_chain", "algebraic_chain"])
@pytest.mark.parametrize("mutate", [pytest.param(f, id=n) for n, f in TREE_MUTATIONS.items()])
def test_tree_changes(evidence, kind, mutate):
    assert mutate(evidence(kind)).verify() is False


@pytest.mark.parametrize(
    "options,reason",
    [
        ({"max_depth": 0}, "max_depth"),
        ({"max_nodes": 1}, "max_nodes"),
        ({"max_degree": 4}, "degree_limit"),
        ({"max_candidates": 0}, "search_incomplete"),
    ],
)
def test_open_leaf_changes(options, reason):
    chain = rational_chain(x**8, x, **options)
    assert chain.verify() and not chain.complete and chain.root.stop_reason == reason
    for change in (
        {"stop_reason": "indecomposable", "detail": None},
        {"detail": ""},
        {"stop_reason": "repeated_component"},
    ):
        assert replace(chain, root=replace(chain.root, **change)).verify() is False


def test_diagnostic_scope(evidence):
    chain = rational_chain(x**8, x, max_depth=0)
    diagnostic = replace(
        chain,
        root=replace(chain.root, stop_reason="unsupported", detail="exact operation unavailable"),
    )
    assert diagnostic.verify() and not diagnostic.complete and diagnostic.status == "incomplete"
    assert (
        replace(
            diagnostic, root=replace(diagnostic.root, stop_reason="indecomposable", detail=None)
        ).verify()
        is False
    )
    search = evidence("algebraic_search")
    attempt = replace(search.attempts[0], result=None, error="exact operation unavailable")
    unresolved = replace(search, attempts=(attempt, *search.attempts[1:]))
    assert unresolved.verify() and not unresolved.complete and unresolved.status == "incomplete"


def test_repeated_leaf_changes():
    chain = algebraic_chain(x**3 - 3 * x - z, (x, z), max_depth=4)
    assert chain.verify() and not chain.complete

    def remove_cycle(node):
        if node.stop_reason == "repeated_component":
            return replace(node, stop_reason="indecomposable", detail=None)
        if node.split is None:
            return node
        return replace(node, outer=remove_cycle(node.outer), inner=remove_cycle(node.inner))

    assert any(leaf.stop_reason == "repeated_component" for leaf in chain.leaves)
    assert replace(chain, root=remove_cycle(chain.root)).verify() is False
