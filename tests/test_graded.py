"""Hilbert invariants and saturated projective closure with independent expectations."""

from itertools import product

import pytest
import sympy as sp

from alggeom.graded import (
    degree,
    hilbert_function,
    hilbert_polynomial,
    hilbert_series,
    is_homogeneous,
    krull_dimension,
    projective_closure,
    projective_dimension,
)
from alggeom.ideals import PolynomialRing

x, y, z, w, h, q, n = sp.symbols("x y z w h q n")


def test_twisted_cubic():
    ring = PolynomialRing((x, y, z, w))
    ideal = ring.ideal(x * z - y * y, x * w - y * z, y * w - z * z)
    assert sp.cancel(hilbert_series(ideal, parameter=q) - (1 + 2 * q) / (1 - q) ** 2) == 0
    assert hilbert_polynomial(ideal, parameter=n) == sp.Poly(3 * n + 1, n, domain=sp.QQ)
    assert [hilbert_function(ideal, i) for i in range(5)] == [1, 4, 7, 10, 13]
    assert krull_dimension(ideal) == 2 and projective_dimension(ideal) == 1 and degree(ideal) == 3


def test_artinian():
    ring = PolynomialRing((x, y))
    ideal = ring.ideal(x * x, y**3)
    assert sp.expand(hilbert_series(ideal, parameter=q) - (1 + q) * (1 + q + q * q)) == 0
    assert hilbert_polynomial(ideal, parameter=n).is_zero
    assert [hilbert_function(ideal, i) for i in range(5)] == [1, 2, 2, 1, 0]
    assert degree(ideal) == 6 and krull_dimension(ideal) == 0 and projective_dimension(ideal) == -1
    assert ring.quotient(*ideal.generators).vector_space_dimension == 6


def test_projective_closure():
    ring = PolynomialRing((x, y, z))
    affine = ring.ideal(y - x * x, z - x**3)
    closure = projective_closure(affine, homogenizer=h)
    assert closure.ring.variables == (x, y, z, h)
    assert closure == closure.ring.ideal(y * h - x * x, z * h - x * y, x * z - y * y)
    assert degree(closure) == 3 and projective_dimension(closure) == 1
    assert closure.quotient(closure.ring.ideal(h)) == closure
    assert is_homogeneous(closure) and closure.saturate(h) == closure


@pytest.mark.parametrize(
    "generators", [(x * x, y**3), (x * y, x * z, y * z), (x**2 * y, y * y * z, z**3)]
)
def test_monomial_counts(generators):
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal(*generators)
    for d in range(7):
        monomials = [
            sp.prod(v**e for v, e in zip(ring.variables, exponents, strict=True))
            for exponents in product(range(d + 1), repeat=3)
            if sum(exponents) == d
        ]
        expected = sum(
            not any(sp.div(m, g, *ring.variables)[1] == 0 for g in generators) for m in monomials
        )
        assert hilbert_function(ideal, d) == expected


def test_boundaries():
    ring = PolynomialRing((x, y))
    assert hilbert_series(ring.ideal(1), parameter=q) == 0 and degree(ring.ideal(1)) == 0
    assert hilbert_polynomial(ring.ideal(1), parameter=n).is_zero
    assert krull_dimension(ring.ideal(x * y - 1)) == 1
    with pytest.raises(ValueError, match="homogeneous"):
        hilbert_series(ring.ideal(x * y - 1))
    assert is_homogeneous(ring.ideal(x + 1, x))
    scalar = PolynomialRing(())
    assert hilbert_series(scalar, parameter=q) == 1 and hilbert_function(scalar, 0) == 1
    assert hilbert_function(scalar, 1) == 0 and degree(scalar) == 1
    assert projective_closure(scalar, homogenizer=h).ring.variables == (h,)
    with pytest.raises(ValueError):
        hilbert_function(ring, -1)


def test_nonreduced_closure():
    ring = PolynomialRing((x, y), sp.QQ.algebraic_field(sp.sqrt(-3)))
    closure = projective_closure(ring.ideal(x * x, y), homogenizer=h)
    assert degree(closure) == 2 and projective_dimension(closure) == 0
    assert hilbert_polynomial(closure, parameter=n) == sp.Poly(2, n, domain=sp.QQ)
    assert closure == closure.ring.ideal(x * x, y)
