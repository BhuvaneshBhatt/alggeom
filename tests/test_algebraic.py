"""Supplied factor reconstruction and exact resultant identities."""

from dataclasses import replace

import pytest
import sympy as sp
from algroots.quotient import QuotientAlgebra

from alggeom.decomposition import DecompositionError
from alggeom.decomposition.algebraic import AlgebraicDecomposition, reconstruct_algebraic
from alggeom.decomposition.rational import RationalDecomposition

x, y, z, a = sp.symbols("x y z a")


def test_coefficient_subfield():
    F = x**4 - z * x * x * (x + 1) + z**3 * (x + 1) ** 2
    A = QuotientAlgebra.from_polynomials([a * a - z * a + z**3], (a,), domain=sp.QQ.frac_field(z))
    result = reconstruct_algebraic(F, (x, z), candidate=x * x - a * (x + 1), ambient=A, parameter=y)
    assert result.status == "decomposable" and (not result.complete) and result.verify()
    assert not result.coefficient_field.contains(z)
    D = result.decomposition
    assert D.verify() and D.exponent == 1
    assert D.outer == sp.Poly(x * x + y * (x + 1), x, y, domain=sp.QQ)
    assert D.inner == sp.Poly(y * y + z * y + z**3, y, z, domain=sp.QQ)
    assert not replace(D, inner=sp.Poly(y * y + z * y + z**3 + 1, y, z)).verify()
    assert not replace(D, unit=z).verify()
    assert not replace(result, candidate=x * x + a * (x + 1)).verify()
    assert not replace(result, reason="positive_genus").verify()
    assert not replace(result, decomposition=replace(D, exponent=2)).verify()


def test_resultant_power():
    A = QuotientAlgebra.from_polynomials([a**3 - 3 * a - z], (a,), domain=sp.QQ.frac_field(z))
    result = reconstruct_algebraic(
        x**3 - 3 * x - z, (x, z), candidate=x * x + a * x + a * a - 3, ambient=A, parameter=y
    )
    assert result.verify() and result.decomposition.exponent == 2
    assert not replace(result.decomposition, exponent=1).verify()
    assert not replace(result.decomposition, unit=2).verify()


def test_full_factor():
    A = QuotientAlgebra.from_polynomials([a - z], (a,), domain=sp.QQ.frac_field(z))
    F = x * x - z**4 - 1
    result = reconstruct_algebraic(F, (x, z), candidate=F, ambient=A, parameter=y)
    assert result.verify() and result.decomposition.inner.degree(z) == 4


@pytest.mark.parametrize(
    "modulus,reason",
    [(a * a - z**3 + z, "positive_genus"), (a * a + z * z - 3, "no_ground_field_point")],
)
def test_candidate_rejection(modulus, reason):
    A = QuotientAlgebra.from_polynomials([modulus], (a,), domain=sp.QQ.frac_field(z))
    g = x * x + a * x + z
    F = sp.resultant(g, modulus, a)
    result = reconstruct_algebraic(F, (x, z), candidate=g, ambient=A, parameter=y)
    assert result.status == "candidate_rejected" and result.reason == reason
    assert not result.complete and (not result.decompositions) and result.verify()
    assert not replace(result, reason="other").verify()


def test_trivial_components():
    A = QuotientAlgebra.from_polynomials([a**3 - z], (a,), domain=sp.QQ.frac_field(z))
    result = reconstruct_algebraic(x**3 - z, (x, z), candidate=x**3 - a**3, ambient=A, parameter=y)
    assert result.reason == "trivial_component" and result.verify()
    linear = reconstruct_algebraic(x**3 - z, (x, z), candidate=x - a, ambient=A, parameter=y)
    assert linear.reason == "trivial_component" and linear.verify()


@pytest.mark.parametrize(
    "options", [{"candidate": x * x + a}, {"candidate": 1}, {"parameter": x}, {"ambient": None}]
)
def test_invalid_candidates(options):
    A = QuotientAlgebra.from_polynomials([a * a - z], (a,), domain=sp.QQ.frac_field(z))
    kw = {"candidate": x - a, "ambient": A, "parameter": y}
    kw.update(options)
    with pytest.raises(DecompositionError):
        reconstruct_algebraic(x * x - z, (x, z), **kw)


def test_invalid_resultant():
    D = AlgebraicDecomposition(
        sp.Poly(x * x - z, x, z), sp.Poly(x * x - y, x, y), sp.Poly(y - z, y, z), 1, sp.Integer(1)
    )
    assert not D.verify()
    bad = replace(D, inner=sp.Poly(y * y - z, y, z))
    assert not bad.verify()
    assert not RationalDecomposition(x**4, (x,), y, sp.QQ, y, x**4).verify()
    assert not replace(bad, unit=sp.Float(1)).verify()


def test_reducible_component():
    A = QuotientAlgebra.from_polynomials([a**4 - z], (a,), domain=sp.QQ.frac_field(z))
    result = reconstruct_algebraic(
        x**4 - z, (x, z), candidate=x**3 + a * x * x + a * a * x + a**3, ambient=A, parameter=y
    )
    assert result.reason == "reducible_candidate" and result.verify()


def test_factor_reduction():
    modulus = a * a - z * a + z**3
    A = QuotientAlgebra.from_polynomials([modulus], (a,), domain=sp.QQ.frac_field(z))
    F = x**4 - z * x * x * (x + 1) + z**3 * (x + 1) ** 2
    g = x * x - a * (x + 1)
    result = reconstruct_algebraic(
        F, (x, z), candidate=modulus * x**3 + 2 * a * g, ambient=A, parameter=y
    )
    assert result.candidate == g.expand() and result.verify()


def test_number_field():
    K = sp.QQ.algebraic_field(sp.sqrt(2))
    A = QuotientAlgebra.from_polynomials([a - z], (a,), domain=K.frac_field(z))
    F = x * x - z**4 - sp.sqrt(2)
    result = reconstruct_algebraic(F, (x, z), candidate=F, ambient=A, parameter=y, constants=K)
    assert result.verify() and result.decomposition.constants == K
