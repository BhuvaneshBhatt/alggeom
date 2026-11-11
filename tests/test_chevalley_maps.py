"""Integration of exact Chevalley projection with polynomial and rational maps."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.maps import RationalMap
from alggeom.schemes import AffineScheme


def test_polynomial_nonlinear_map():
    s, t, x, y = sp.symbols("s t x y")
    source = AffineScheme(PolynomialRing((s, t)))
    target = AffineScheme(PolynomialRing((x, y)))
    image = source.morphism(target, (s**2, s * t)).constructible_image.as_constructible()
    assert image.contains((0, 0))
    assert image.contains((4, 3))
    assert not image.contains((0, 3))
    assert image.closure_ideal.is_zero


def test_rational_map_exact_image():
    s, u = sp.symbols("s u")
    image = RationalMap(PolynomialRing((s,)), (1 / s,), (u,)).constructible_image.as_constructible()
    assert image.contains((1,))
    assert image.contains((-2,))
    assert not image.contains((0,))


def test_identical_source_target_symbols():
    x, y = sp.symbols("x y")
    source = AffineScheme(PolynomialRing((x, y)))
    target = AffineScheme(PolynomialRing((x,)))
    image = source.morphism(target, (x * x,)).constructible_image.as_constructible()
    assert image.contains((2,))
    assert image.contains((0,))
