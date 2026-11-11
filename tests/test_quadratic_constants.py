"""Reproducible completeness and adversarial checks for quadratic extensions."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.function_field_constants import algebraic_constant_field, certify_geometric_prime
from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing


@pytest.mark.parametrize("exponent", [1, 4, 7])
def test_unbounded_primitive_degree(exponent):
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    g = x**exponent + 1
    prime = ring.ideal(g**2 * y**2 - 2)
    result = algebraic_constant_field(prime, numerator_degree=1, max_candidates=1)
    assert result.field.mod.degree() == 2
    assert sp.cancel(result.numerator / result.denominator - 2 * g * y) == 0
    assert len(result.geometric_primes) == 2
    assert result.verify(prime)
    assert not certify_geometric_prime(prime)
    assert all(certify_geometric_prime(component) for component in result.geometric_primes)
    assert not replace(result, denominator=x).verify(prime)
    assert not replace(result, geometric_primes=result.geometric_primes[::-1][:1] * 2).verify(prime)


@pytest.mark.parametrize(
    "expression",
    [
        lambda x, y: y**2 - x**11 - 3,
        lambda x, y: y**2 + x**8 + 1,
        lambda x, y: (x**5 + 1) * y**2 + x * y + 3,
    ],
)
def test_discriminant_geometric_completeness(expression):
    x, y = sp.symbols("x y")
    prime = PolynomialRing((x, y)).ideal(expression(x, y))
    result = algebraic_constant_field(prime, max_candidates=1)
    assert result.field == sp.QQ
    assert result.verify(prime)
    assert certify_geometric_prime(prime)
    assert not replace(result, numerator=x).verify(prime)
    assert not replace(result, minimal_polynomial=sp.Symbol("z") + 1).verify(prime)


def test_embedded_primary_replay():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    f = y**2 - 2 * x**8
    ideal = ring.ideal(f * x, f * y)
    decomposition = absolute_primary_decomposition(ideal)
    assert len(decomposition.components) == 3
    assert sorted(component.krull_dimension for component in decomposition.components) == [0, 1, 1]
    assert decomposition.verify(ideal)
    assert not replace(decomposition, components=decomposition.components[:2]).verify(ideal)
    assert not replace(decomposition, components=decomposition.components * 2).verify(ideal)
