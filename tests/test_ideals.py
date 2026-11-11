"""Ideal identities, elimination contractions and coordinate-algebra dimensions."""

import pytest
import sympy as sp

from alggeom.ideals import CoordinateRing, PolynomialRing

x, y, z, t = sp.symbols("x y z t")


@pytest.mark.parametrize("order", ["lex", "grlex", "grevlex"])
def test_membership(order):
    ring = PolynomialRing((x, y), order=order)
    ideal = ring.ideal(x * y - 1, y * y - y)
    f = x**2 - 1
    coefficients = ideal.membership(f)
    assert (
        sp.expand(sum(a * g for a, g in zip(coefficients, ideal.generators, strict=True)) - f) == 0
    )
    assert ideal.normal_form(x) == 1 and ideal.normal_form(y) == 1
    assert ideal.membership(x + 1) is None
    assert ideal == ideal.change_ring(ring.with_order("lex"))


def test_ideal_operations():
    ring = PolynomialRing((x, y))
    left, right = ring.ideal(x), ring.ideal(y)
    assert ring.ideal(x, y) == left + right
    assert left * right == ring.ideal(x * y) == left.intersection(right)
    assert ring.ideal(x * y, y * y).quotient(right) == ring.ideal(x, y)
    assert ring.ideal(x * y).saturate(x) == right
    assert ring.ideal(x * y).saturate(ring.ideal(x, y)) == ring.ideal(x * y)
    assert ring.ideal(x * y).saturate(x * y).is_unit
    assert left.quotient(ring.ideal()).is_unit and left.saturate(ring.ideal()).is_unit


def test_elimination():
    ring = PolynomialRing((t, x, y))
    implicit = ring.ideal(x - t * t, y - t**3).eliminate(t)
    assert implicit.ring.variables == (x, y)
    assert implicit == implicit.ring.ideal(y * y - x**3)
    point = PolynomialRing((x,)).ideal(x - 1).eliminate(x)
    assert point.ring.variables == () and point.is_zero
    impossible = PolynomialRing((x,)).ideal(1).eliminate(x)
    assert impossible.is_unit


def test_dimensions():
    ring = PolynomialRing((x, y, z))
    assert ring.ideal().krull_dimension == 3
    assert ring.ideal(x * y, x * z).krull_dimension == 2
    assert ring.ideal(1).krull_dimension == -1
    affine = CoordinateRing(ring.ideal(x * y))
    assert affine.krull_dimension == 2 and affine.vector_space_dimension == sp.oo
    with pytest.raises(ValueError, match="positive-dimensional"):
        affine.finite_algebra()
    finite = CoordinateRing(ring.ideal(x * x, y * y, z))
    assert finite.krull_dimension == 0 and finite.vector_space_dimension == 4
    assert finite.finite_algebra().dimension == 4
    scalar = PolynomialRing(()).quotient()
    assert scalar.krull_dimension == 0 and scalar.vector_space_dimension == 1
    assert scalar.finite_algebra().dimension == 1
    zero = PolynomialRing(()).quotient(1)
    assert zero.krull_dimension == -1 and zero.vector_space_dimension == 0


@pytest.mark.parametrize(
    "constants",
    [
        sp.QQ,
        sp.QQ.algebraic_field(sp.sqrt(-3)),
        sp.QQ.frac_field(t),
        sp.QQ.algebraic_field(sp.sqrt(-3)).frac_field(t),
    ],
)
def test_coefficient_fields(constants):
    alpha = (
        sp.sqrt(-3)
        if constants.is_AlgebraicField
        or (constants.is_FractionField and constants.domain.is_AlgebraicField)
        else 1
    )
    parameter = t if constants.is_FractionField else 1
    ring = PolynomialRing((x, y), domain=constants)
    ideal = ring.ideal(x - alpha * y, y * y - parameter)
    assert ideal.contains(x * x - alpha**2 * parameter)
    assert ideal.krull_dimension == 0
    witness = ideal.membership(x * x - alpha**2 * parameter)
    assert (
        ring.convert(
            sum(a * g for a, g in zip(witness, ideal.generators, strict=True))
            - x * x
            + alpha**2 * parameter
        )
        == 0
    )


def test_validation():
    ring = PolynomialRing((x, y))
    with pytest.raises(ValueError):
        PolynomialRing((x, x))
    with pytest.raises(ValueError):
        PolynomialRing((x,), domain=sp.QQ.frac_field(x))
    with pytest.raises(ValueError):
        ring.ideal(sp.Float(0.5) * x)
    with pytest.raises(sp.PolynomialError):
        ring.ideal(1 / x)
    with pytest.raises(ValueError):
        ring.ideal(x).eliminate(z)
