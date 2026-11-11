"""Exact algebraic constant fields and geometric-primality certificates."""

import pytest
import sympy as sp

from alggeom.function_field_constants import (
    algebraic_constant_field,
    certify_geometric_prime,
)
from alggeom.geometric import geometric_components, is_geometrically_integral
from alggeom.ideals import PolynomialRing


def test_nonlinear_smooth_prime():
    x, y, z = sp.symbols("x y z")
    prime = PolynomialRing((x, y, z), sp.QQ, "lex").ideal(y**2 - x**3 + x, z - x * y)
    assert certify_geometric_prime(prime)
    constants = algebraic_constant_field(prime)
    assert constants.field == sp.QQ
    assert constants.geometric_primes == (prime,)
    assert is_geometrically_integral(prime)


def test_quadratic_constant_field():
    x, y, z = sp.symbols("x y z")
    prime = PolynomialRing((x, y, z), sp.QQ, "lex").ideal(x**2 + y**2, z - y**2)
    constants = algebraic_constant_field(prime, max_candidates=12)
    assert sp.Poly(constants.minimal_polynomial).all_coeffs() == [1, 0, 1]
    assert constants.numerator == x
    assert constants.denominator == y
    assert constants.field.mod.degree() == 2
    assert len(constants.geometric_primes) == 2
    for component in constants.geometric_primes:
        assert certify_geometric_prime(component)
    assert not is_geometrically_integral(prime)
    assert len(geometric_components(prime).components) == 2


def test_nonnormal_constants():
    x, y, z = sp.symbols("x y z")
    prime = PolynomialRing((x, y, z), sp.QQ, "lex").ideal(x**3 - 2 * y**3, z - y**2)
    constants = algebraic_constant_field(prime, max_candidates=12, max_degree=3)
    assert sp.Poly(constants.minimal_polynomial).all_coeffs() == [1, 0, 0, -2]
    assert constants.field.mod.degree() == 3
    assert constants.splitting_field.mod.degree() == 6
    assert len(constants.geometric_primes) == 3


def test_nonprime_rejected():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    assert not certify_geometric_prime(ring.ideal(x * y))
    with pytest.raises(ValueError, match="not certified prime"):
        algebraic_constant_field(ring.ideal(x * y))


def test_constant_search_inconclusive():
    x, y, z = sp.symbols("x y z")
    prime = PolynomialRing((x, y, z), sp.QQ, "lex").ideal(x**2 + y**2, z - y**2)
    with pytest.raises(NotImplementedError, match="could not be certified"):
        algebraic_constant_field(prime, max_candidates=1)


def test_smooth_projective_conic():
    x, y = sp.symbols("x y")
    prime = PolynomialRing((x, y), sp.QQ, "lex").ideal(x**2 + y**2 + 1)
    assert certify_geometric_prime(prime)
    assert algebraic_constant_field(prime).field == sp.QQ


def test_constant_field_primary_splitting():
    from alggeom.geometric import _constant_field_primary_decomposition

    x, y, z = sp.symbols("x y z")
    prime = PolynomialRing((x, y, z), sp.QQ, "lex").ideal(x**2 + y**2, z - y**2)
    result = _constant_field_primary_decomposition(prime, max_components=8)
    assert len(result.components) == 2
    assert all(certify_geometric_prime(component) for component in result.components)
    intersection = result.components[0].intersection(result.components[1])
    assert intersection == prime.change_ring(result.components[0].ring)
