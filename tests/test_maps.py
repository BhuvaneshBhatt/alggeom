"""Image closures, retained domains and tangent/singular scheme structure."""

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.maps import (
    RationalMap,
    graph_ideal,
    implicitize,
    singular_locus,
    tangent_ideal,
    tangent_space,
)

t, x, y, u, v = sp.symbols("t x y u v")


def test_cusp_image():
    image = implicitize((t * t, t**3), (t,), targets=(x, y))
    assert image == image.ring.ideal(y * y - x**3)
    graph = graph_ideal((t * t, t**3), (t,), targets=(x, y))
    assert graph.contains(y * y - x**3) and not graph.contains(t)
    assert image.krull_dimension == 1


def test_rational_domains():
    ring = PolynomialRing((t,))
    inverse = RationalMap(ring, (1 / t,), (x,))
    assert inverse.image_closure.is_zero and inverse.graph_ideal == inverse.graph_ideal.ring.ideal(
        t * x - 1
    )
    assert inverse.evaluate((2,)) == (sp.Rational(1, 2),)
    with pytest.raises(ValueError, match="denominator"):
        inverse.evaluate((0,))
    original = sp.Mul(t * t - 1, sp.Pow(t - 1, -1, evaluate=False), evaluate=False)
    punctured = RationalMap(ring, (original,), (x,))
    assert punctured.functions == (t + 1,) and punctured.domain_guard.subs(t, 1) is sp.false
    assert punctured.image_closure.is_zero
    with pytest.raises(ValueError, match="denominator"):
        punctured.evaluate((1,))


def test_source_and_composition():
    source = PolynomialRing((x, y)).quotient(x * y - 1)
    parametrization = RationalMap(PolynomialRing((t,)), (t, 1 / t), (u, v))
    outer = RationalMap(source, (x + y,), (t,))
    composed = outer.compose(parametrization)
    assert composed.evaluate((2,)) == (sp.Rational(5, 2),)
    assert composed.domain_guard.subs(t, 0) is sp.false
    invalid = RationalMap(PolynomialRing((t,)), (t, t), (u, v))
    with pytest.raises(ValueError, match="does not lie"):
        outer.compose(invalid)
    empty = RationalMap(PolynomialRing((t,)).quotient(t), (1 / t,), (x,))
    assert empty.graph_ideal.is_unit and empty.image_closure.is_unit


def test_singularities():
    ring = PolynomialRing((x, y))
    cusp = ring.ideal(y * y - x**3)
    assert singular_locus(cusp) == ring.ideal(x * x, y)
    assert tangent_space(cusp, (0, 0)) == sp.eye(2)
    assert tangent_ideal(cusp, (0, 0)).is_zero
    assert tangent_space(cusp, (1, 1)).shape == (2, 1)
    assert tangent_ideal(cusp, (1, 1)) == ring.ideal(-3 * (x - 1) + 2 * (y - 1))
    assert singular_locus(ring.ideal(y - x * x)).is_unit
    assert singular_locus(ring.ideal(x * x, y)).contains(x)
    with pytest.raises(ValueError, match="satisfy"):
        tangent_space(cusp, (0, 1))


def test_mixed_dimension_boundary():
    ring = PolynomialRing((x, y, t))
    ideal = ring.ideal(x * y, x * t)
    with pytest.raises(ValueError, match="mixed-dimensional"):
        singular_locus(ideal)
    assert singular_locus(ideal, reduced=True) == ring.ideal(x, y, t)
    assert tangent_space(ideal, (0, 0, 0)) == sp.eye(3)


def test_number_field_map():
    alpha = sp.sqrt(-3)
    field = sp.QQ.algebraic_field(alpha).frac_field(t)
    source = PolynomialRing((u,), field)
    image = implicitize((u**2 + alpha, u**3), u, targets=(x, y), source=source)
    assert image.ring.domain == field
    assert image == image.ring.ideal(y * y - (x - alpha) ** 3)
    assert tangent_space(image, (alpha, 0)) == sp.eye(2)
    assert singular_locus(image) == image.ring.ideal((x - alpha) ** 2, y)


def test_composition_into_pole():
    outer = RationalMap(PolynomialRing((y,)), (1 / y,), (u,))
    inner = RationalMap(PolynomialRing((x,)), (sp.Integer(0),), (y,))
    with pytest.raises(ValueError, match="empty domain"):
        outer.compose(inner)


def test_composition_with_nonempty_pole_avoidance():
    outer = RationalMap(PolynomialRing((y,)), (1 / y,), (u,))
    inner = RationalMap(PolynomialRing((x,)), (x,), (y,))
    composed = outer.compose(inner)
    assert composed.evaluate((2,)) == (sp.Rational(1, 2),)
    with pytest.raises(ValueError, match="excluded"):
        composed.evaluate((0,))
