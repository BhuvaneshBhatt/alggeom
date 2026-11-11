"""Localization identities in reduced, reducible, and nonreduced quotient rings."""

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.parameter_localization import ParameterLocalization


def test_zero_divisor_is_killed_by_localization():
    a, b = sp.symbols("a b")
    ring = PolynomialRing((a, b))
    domain = ParameterLocalization(ring.ideal(a * b))
    assert domain.element(a * b).is_zero()
    assert domain.element(a) != 0
    assert domain.element(b) != 0
    opened = domain.localize_at(a)
    assert opened.element(b).is_zero()  # V(ab) intersect D(a) = V(b) intersect D(a)
    inv = opened.inverse_generator(0)
    assert opened.element(a) * inv == 1
    assert (opened.element(a + 2) * inv) == (opened.element(1) + 2 * inv)
    assert opened.restrict(domain.element(a)) * inv == 1
    assert opened.element(a) + opened.element(b) == opened.element(a)
    with pytest.raises(ZeroDivisionError, match="empty"):
        opened.localize_at(b)


def test_nilpotent_zero_and_fraction_arithmetic():
    a, b = sp.symbols("a b")
    ring = PolynomialRing((a, b))
    base = ParameterLocalization(ring.ideal(a**2 * b), (a,))
    inv = base.inverse_generator(0)
    assert base.element(b).is_zero()
    assert (base.element(a) * inv) == 1
    assert (inv + 1) * base.element(a) == base.element(1 + a)
    assert inv**2 * base.element(a**2) == 1
    assert base.element(a) != 0
    with pytest.raises(ValueError, match="nonnegative"):
        base.element(1, (-1,))
    with pytest.raises(TypeError, match="different localizations"):
        base.coerce(ParameterLocalization(ring.ideal(a**2 * b)).element(a))


def test_reduction_on_nonlinear_implicit_quotient():
    a, b = sp.symbols("a b")
    ring = PolynomialRing((a, b))
    opened = ParameterLocalization(ring.ideal(b**2 - a**3), (a,))
    assert opened.element(b**2) == opened.element(a**3)
    assert opened.element(b) ** 2 * opened.inverse_generator(0) ** 3 == 1
    further = opened.localize_at(b)
    assert further.restrict(opened.inverse_generator(0)) * further.element(a) == 1
    with pytest.raises(TypeError, match="prefix open"):
        opened.restrict(further.inverse_generator(1))


def test_nonradical_branching():
    from alggeom.cover_proofs import check_parameter_partition

    a, b = sp.symbols("a b")
    ring = PolynomialRing((a, b))
    ambient = ParameterLocalization(ring.ideal())
    opened, closed = ambient.split_ideal(ring.ideal(a * b, a**2))
    assert len(opened) == 2
    assert check_parameter_partition((*opened, closed), (a, b))
    # The surviving closed component has a=0, including the nonreduced scheme.
    assert closed.ideal.contains(a**2)
    assert closed.ideal.contains(a * b)
    assert not closed.ideal.contains(a)
    # On a=0, an exceptional ideal (a,b) can be branched again without
    # needing a prime decomposition of the scheme or a field of fractions.
    further_open, further_closed = closed.split_ideal(ring.ideal(a, b))
    assert len(further_open) == 1
    opened_b, vanished_b = further_open[0], further_closed
    assert opened_b.condition.has(sp.Ne(b, 0))
    assert vanished_b.ideal == ring.ideal(a, b)


def test_vanishing_guard_returns_empty_open():
    a, b = sp.symbols("a b")
    ring = PolynomialRing((a, b))
    parent = ParameterLocalization(ring.ideal(a))
    opened, closed = parent.split_at(a)
    assert opened is None
    assert closed.ideal == parent.ideal
