"""Automatic factor coefficient extensions and exhaustive orbit accounting."""

from dataclasses import replace

import pytest
import sympy as sp
from algroots.subfields import field_normal_form

from alggeom.decomposition import DecompositionError, DegreeLimitError
from alggeom.decomposition.discovery import CandidateAttempt, search_algebraic
from alggeom.decomposition.factors import FactorIssue, discover_factors

x, z, t = sp.symbols("x z t")
PAPER = x**4 - z * x * x * (x + 1) + z**3 * (x + 1) ** 2


def test_factor_orbits():
    discovery = discover_factors(PAPER, (x, z), degrees=(2,))
    assert discovery.complete and discovery.verify()
    (level,) = discovery.levels
    assert level.geometric_factor_count == 6
    assert [o.extension_degree for o in level.orbits] == [2, 4]
    assert sum(o.extension_degree for o in level.orbits) == 6
    for orbit in level.orbits:
        assert orbit.verify()
        difference = sp.Poly(orbit.factor * orbit.cofactor - PAPER, x)
        assert all(field_normal_form(orbit.ambient, c) == 0 for c in difference.all_coeffs())
    assert not replace(
        level, orbits=level.orbits[:1], component_polynomials=level.component_polynomials[:1]
    ).verify()
    bad = level.coordinate_polynomials[0] + 1
    assert not replace(
        level, coordinate_polynomials=(bad, *level.coordinate_polynomials[1:])
    ).verify()
    assert not replace(level.orbits[0], cofactor=level.orbits[0].cofactor + 1).verify()
    assert not replace(discovery, requested_degrees=(2, 3)).verify()


def test_automatic_reconstruction():
    result = search_algebraic(PAPER, (x, z), parameter=t)
    assert result.status == "decomposable" and result.verify()
    assert result.examined == 1 and result.candidate_count == 4
    assert result.discovery.complete and (not result.complete)
    (D,) = result.decompositions
    assert D.verify() and D.exponent == 1
    assert not result.attempts[0].result.coefficient_field.contains(z)
    bad = replace(result.attempts[0].result, decomposition=replace(D, unit=2))
    assert not replace(result, attempts=(CandidateAttempt(0, bad),)).verify()
    assert not replace(result, attempts=(replace(result.attempts[0], index=1),)).verify()


def test_coefficient_extension():
    discovery = discover_factors(x**4 + x + z, (x, z), degrees=(2,))
    assert discovery.complete and discovery.verify()
    (orbit,) = discovery.orbits
    assert orbit.extension_degree == 6 and orbit.verify()
    assert 4 % orbit.extension_degree != 0


def test_resultant_power():
    result = search_algebraic(x**3 - 3 * x - z, (x, z), parameter=t, exhaustive=True)
    assert result.status == "decomposable" and result.complete and result.verify()
    assert result.examined == result.candidate_count == 2
    assert result.decompositions[0].exponent == 2
    assert result.attempts[1].result.reason == "trivial_component"


@pytest.mark.parametrize("polynomial", [x * x + z * x + z, x**3 + z * x + z**3 + z])
def test_exhaustive_rejections(polynomial):
    result = search_algebraic(polynomial, (x, z), parameter=t, exhaustive=True)
    assert result.complete and result.status == "indecomposable" and result.verify()
    assert all(a.result.reason in {"positive_genus", "trivial_component"} for a in result.attempts)
    assert not replace(result.discovery, levels=()).verify()


def test_partial_coverage():
    F = x**3 + z * x + z**3 + z
    limited = search_algebraic(F, (x, z), max_candidates=0)
    assert limited.status == "incomplete" and (not limited.complete) and limited.verify()
    limited = search_algebraic(F, (x, z), max_candidates=1, exhaustive=True)
    assert limited.status == "incomplete" and limited.verify()
    scoped = search_algebraic(F, (x, z), factor_degrees=(2,), exhaustive=True)
    assert scoped.discovery.complete and scoped.status == "incomplete" and scoped.verify()
    assert not scoped.complete


def test_dimension_limit():
    discovery = discover_factors(PAPER, (x, z), degrees=(2, 3, 4), max_factor_dimension=5)
    assert not discovery.complete and discovery.verify()
    assert discovery.issues == (FactorIssue(2, "factor_dimension_limit"),)
    assert [level.degree for level in discovery.levels] == [3, 4]
    assert not replace(discovery, max_factor_dimension=6).verify()
    result = search_algebraic(PAPER, (x, z), max_factor_dimension=1, max_candidates=0)
    assert result.status == "incomplete" and result.verify()
    with pytest.raises(DegreeLimitError):
        discover_factors(x**7 + z, (x, z))


def test_separating_form():
    F = PAPER.subs(x, -x)
    failed = discover_factors(F, (x, z), degrees=(2,), max_separating_attempts=1)
    assert not failed.complete and failed.verify()
    assert failed.issues == (FactorIssue(2, "separator_limit"),)
    retried = discover_factors(F, (x, z), degrees=(2,), max_separating_attempts=2)
    assert retried.complete and retried.verify()


def test_nonmonic_input():
    F = z * x * x + x + 1
    result = discover_factors(F, (x, z))
    assert result.requested_degrees == (1, 2) and result.complete and result.verify()
    for orbit in result.orbits:
        assert all(
            field_normal_form(orbit.ambient, c) == 0
            for c in sp.Poly(z * orbit.factor * orbit.cofactor - F, x).all_coeffs()
        )


def test_fields_and_symbols():
    K = sp.QQ.algebraic_field(sp.sqrt(2))
    result = discover_factors(x * x + sp.sqrt(2) * z * x + z, (x, z), constants=K)
    assert result.complete and result.verify()
    assert all(o.ambient.domain == K.frac_field(z) for o in result.orbits)
    u = sp.Dummy("factor_extension", positive=True)
    v = sp.Symbol("factor_coefficient_0", real=True)
    result = discover_factors(u * u + v * u + v, (u, v))
    assert result.complete and result.verify()


def test_linear_output():
    result = search_algebraic(x - z**4, (x, z))
    assert result.status == "degree_below_two" and result.complete and result.verify()


@pytest.mark.parametrize("degrees", [(0,), (2, 1), (2, 2), (5,), (True,)])
def test_invalid_degrees(degrees):
    with pytest.raises(DecompositionError):
        discover_factors(PAPER, (x, z), degrees=degrees)


@pytest.mark.parametrize(
    "options", [{"max_factor_dimension": 0}, {"max_separating_attempts": 0}, {"max_degree": True}]
)
def test_invalid_budgets(options):
    with pytest.raises(DecompositionError):
        discover_factors(PAPER, (x, z), **options)


@pytest.mark.parametrize(
    "options",
    [{"max_candidates": -1}, {"max_candidates": True}, {"exhaustive": 1}, {"factor_degrees": (1,)}],
)
def test_invalid_options(options):
    with pytest.raises(DecompositionError):
        search_algebraic(PAPER, (x, z), **options)
