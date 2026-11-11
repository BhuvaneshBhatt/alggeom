"""Splitting nontrivial constant fields in positive-dimensional associated primes."""

import sympy as sp

from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing


def _reconstruct(decomposition):
    ring = decomposition.components[0].ring
    intersection = ring.ideal(1)
    for primary in decomposition.components:
        intersection = intersection.intersection(primary)
    return intersection


def test_positive_dimensional_absolute_split():
    x, y, z = sp.symbols("x y z")
    ideal = PolynomialRing((x, y, z)).ideal(x * x + 1, y * z)
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 4
    assert split.field.to_sympy(split.field.from_sympy(sp.I)) == sp.I
    assert _reconstruct(split) == ideal.change_ring(split.components[0].ring)
    assert all(primary.krull_dimension == 1 for primary in split.components)


def test_embedded_primary_split():
    x, y, z = sp.symbols("x y z")
    ideal = PolynomialRing((x, y, z)).ideal(x * x + 1, y * y, y * z)
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 4
    assert _reconstruct(split) == ideal.change_ring(split.components[0].ring)
    assert sorted(primary.krull_dimension for primary in split.components) == [0, 0, 1, 1]


def test_geometric_integrality():
    x, y = sp.symbols("x y")
    ideal = PolynomialRing((x, y)).ideal(x * x - y)
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 1
    assert split.components[0] == ideal
