"""Absolute-component certificates on nonzero-dimensional algebraic sets."""

import sympy as sp

from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing


def test_nonreduced_coordinate_subspace():
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal(x * x, x * y)
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 2
    intersection = ring.ideal(1)
    for component in split.components:
        intersection = intersection.intersection(component)
    assert intersection == ideal
    assert any(component.krull_dimension == 2 for component in split.components)


def test_smooth_ternary_hypersurface():
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal((x * x + y * y + z * z) ** 2)
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 1
    assert split.components[0] == ideal


def test_two_smooth_factors():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    ideal = ring.ideal((x * x + y * y + 1) * (x + y - 1))
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 2
    assert split.components[0].intersection(split.components[1]) == ideal


def test_linear_support_nonlinear_primary_factors():
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal(x * (x - 1), y * y)
    split = absolute_primary_decomposition(ideal)
    assert len(split.components) == 2
    intersection = ring.ideal(1)
    for component in split.components:
        intersection = intersection.intersection(component)
    assert intersection == ideal


def test_homogeneous_cubic_splitting():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    original = ring.ideal(x**3 + y**3)
    decomposition = absolute_primary_decomposition(original)
    assert len(decomposition.components) == 3
    combined = decomposition.components[0].ring.ideal(1)
    for primary in decomposition.components:
        combined = combined.intersection(primary)
    assert combined == original.change_ring(decomposition.components[0].ring)


def test_binary_form_at_infinity():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    original = ring.ideal(y * (x - y) ** 2)
    decomposition = absolute_primary_decomposition(original)
    assert len(decomposition.components) == 2
    assert decomposition.components[0].intersection(decomposition.components[1]) == original


def test_high_degree_geometric_components():
    from alggeom.geometric import geometric_components, is_geometrically_integral

    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    assert len(geometric_components(ring.ideal(x**3 + y**3)).components) == 3
    assert is_geometrically_integral(ring.ideal(x * x + y * y + 1))
