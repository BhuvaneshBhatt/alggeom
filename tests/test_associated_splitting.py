"""Certified algebraic splitting of coupled positive-dimensional associated primes."""

import sympy as sp

from alggeom.geometric import absolute_primary_decomposition, geometric_components
from alggeom.ideals import PolynomialRing


def test_nonprincipal_associated_splitting():
    x, y, z = sp.symbols("x y z")
    ideal = PolynomialRing((x, y, z)).ideal(x**2 + y**2, x * z, y * z)
    decomposition = absolute_primary_decomposition(ideal)
    assert len(decomposition.components) == 3
    assert all(component.krull_dimension == 1 for component in decomposition.components)
    restored = decomposition.components[0].ring.ideal(1)
    for component in decomposition.components:
        restored = restored.intersection(component)
    assert restored == ideal.change_ring(restored.ring)
    assert len(geometric_components(ideal).components) == 3


def test_embedded_associated_splitting():
    x, y, z = sp.symbols("x y z")
    ideal = PolynomialRing((x, y, z)).ideal(x**2 + y**2, x * z, y * z, z**2)
    decomposition = absolute_primary_decomposition(ideal)
    assert len(decomposition.components) == 3
    assert sorted(comp.krull_dimension for comp in decomposition.components) == [0, 1, 1]
    restored = decomposition.components[0].ring.ideal(1)
    for component in decomposition.components:
        restored = restored.intersection(component)
    assert restored == ideal.change_ring(restored.ring)
    assert len(geometric_components(ideal).components) == 2
