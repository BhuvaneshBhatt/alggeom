"""Coordinate arithmetic and contravariant polynomial algebra maps."""

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.rings import RingHomomorphism

x, y, t = sp.symbols("x y t")


def test_elements():
    ring = PolynomialRing((x, y)).quotient(x * y - 1)
    a, b = ring.gens
    assert a * b == 1 and a.inverse() == b and a**-2 == b * b
    assert (a + 2) * (b - 3) == ring(1 + 2 * y - 3 * x - 6)
    assert x + a == 2 * a and 1 / a == b and a / b == a * a
    assert ring.quotient(a - 1).gens == (ring.quotient(x - 1)(1), ring.quotient(x - 1)(1))
    with pytest.raises(TypeError):
        hash(a)
    with pytest.raises(TypeError):
        a + PolynomialRing((x,)).quotient()(x)
    zero = PolynomialRing((x,)).quotient(1)
    assert zero.zero == zero.one and zero.zero.inverse() == zero.zero


def test_nonunits():
    ring = PolynomialRing((x, y)).quotient(x * y)
    a, b = ring.gens
    assert a * b == 0 and not a.is_unit and (1 + ring(0)).is_unit
    with pytest.raises(ZeroDivisionError):
        a.inverse()
    with pytest.raises(TypeError):
        a ** sp.Rational(1, 2)


def test_maps():
    source = PolynomialRing((x, y)).quotient(y * y - x**3)
    target = PolynomialRing((t,)).quotient()
    hom = source.hom((t * t, t**3), target)
    assert hom(source(x + y)) == target(t * t + t**3)
    assert hom.kernel == source.ideal and hom.is_injective and not hom.is_surjective
    assert hom.lift(t * t) == source(x)
    with pytest.raises(ValueError, match="image"):
        hom.lift(t)
    assert hom.preimage(target.ambient.ideal(t)) == source.ambient.ideal(x, y)
    assert hom.extend(source.ambient.ideal(x, y)) == target.ambient.ideal(t * t)
    with pytest.raises(ValueError, match="relations"):
        source.hom((t, t), target)


def test_inverse():
    ring = PolynomialRing((x, y)).quotient()
    hom = RingHomomorphism(ring, ring, (x + y, y))
    inverse = hom.inverse()
    assert inverse.images == (ring(x - y), ring(y))
    assert hom.compose(inverse).images == ring.gens
    assert inverse.compose(hom).images == ring.gens
    quotient = ring.quotient(x * x, y)
    projection = ring.hom(quotient.gens, quotient)
    assert projection.kernel == quotient.ideal and projection.is_surjective
    assert not projection.is_injective


def test_coefficients():
    ring = PolynomialRing((x,)).quotient(x * x - 2)
    extended = ring.base_change(sp.QQ.algebraic_field(sp.sqrt(2)))
    assert extended(x * x) == 2 and extended.domain == sp.QQ.algebraic_field(sp.sqrt(2))
    hom = extended.hom((sp.sqrt(2),), PolynomialRing((), extended.domain))
    assert hom(x) == hom.target(sp.sqrt(2))
    assert hom.kernel == extended.ambient.ideal(x - sp.sqrt(2))
    with pytest.raises(ValueError, match="coefficient field"):
        ring.hom((sp.sqrt(2),), hom.target)
    with pytest.raises(sp.polys.polyerrors.CoercionFailed):
        extended.base_change(sp.QQ)


@pytest.mark.parametrize("field", [sp.QQ, sp.QQ.algebraic_field(sp.I), sp.QQ.frac_field(t)])
def test_simultaneous_coordinates(field):
    source = PolynomialRing((x, y), field).quotient(x * y - 1)
    target = PolynomialRing((x, y), field).quotient(x * y - 1)
    swap = source.hom((y, x), target)
    assert swap(x + 2 * y) == target(y + 2 * x)
    assert swap.compose(swap).images == source.gens
    with pytest.raises(ValueError, match="source relations"):
        source.hom((y + 1, x), target)
