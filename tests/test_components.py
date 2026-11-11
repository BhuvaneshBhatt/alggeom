from functools import reduce

import pytest
import sympy as sp

from alggeom.components import (
    DecompositionLimitError,
    _saturation_exponent,
    associated_primes,
    is_primary,
    is_prime,
    minimal_primes,
    primary_decomposition,
    radical,
)
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, z, t = sp.symbols("x y z t")


def same_ideals(actual, expected):
    assert len(actual) == len(expected)
    assert all(any(a == b for a in actual) for b in expected)


@pytest.mark.parametrize(
    "generators,components,primes",
    [
        ((x**2, x * y), None, ((x,), (x, y))),
        ((x * y, x * z), ((x,), (y, z)), ((x,), (y, z))),
        ((x**2, x * y, y**2), ((x**2, x * y, y**2),), ((x, y),)),
        ((x**2 * y, x * y**2), None, ((x,), (y,), (x, y))),
        (
            (x * z - y**2, y * z - x**3),
            ((x * z - y**2, y * z - x**3, z**2 - x**2 * y), (x, y)),
            ((x * z - y**2, y * z - x**3, z**2 - x**2 * y), (x, y)),
        ),
    ],
)
def test_primary_components(generators, components, primes):
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal(*generators)
    factors = primary_decomposition(ideal)
    if components is not None:
        same_ideals(factors, tuple(ring.ideal(*g) for g in components))
    else:
        assert len(factors) == len(primes)
    same_ideals(associated_primes(ideal), tuple(ring.ideal(*g) for g in primes))
    assert reduce(lambda a, b: a.intersection(b), factors) == ideal
    for component in factors:
        assert is_primary(component)
        assert is_prime(radical(component))


def test_radicals():
    ring = PolynomialRing((x, y, z))
    ideal = ring.ideal(x**2, x * y)
    assert ideal.radical() == ring.ideal(x)
    same_ideals(ideal.minimal_primes(), (ring.ideal(x),))
    assert not is_prime(ideal)
    assert not is_primary(ideal)
    assert ideal.radical_contains(x)
    assert not ideal.radical_contains(y)
    assert is_prime(ring.ideal(x * y - z**2))
    assert is_primary(ring.ideal((x * y - z**2) ** 3))
    scheme = AffineScheme(ideal)
    assert scheme.reduced().ideal == ring.ideal(x)
    assert tuple(s.ideal for s in scheme.irreducible_components()) == (ring.ideal(x),)
    projective = ProjectiveScheme(ring.ideal(x**2 * y))
    same_ideals(
        tuple(s.ideal for s in projective.irreducible_components()), (ring.ideal(x), ring.ideal(y))
    )


def test_field_relative_components():
    rational = PolynomialRing((x, y))
    ideal = rational.ideal((x**2 + 1) ** 2, y**2)
    assert is_primary(ideal)
    assert radical(ideal) == rational.ideal(x**2 + 1, y)
    extended = PolynomialRing((x, y), sp.QQ.algebraic_field(sp.I))
    changed = ideal.change_ring(extended)
    factors = primary_decomposition(changed)
    same_ideals(
        factors, (extended.ideal((x - sp.I) ** 2, y**2), extended.ideal((x + sp.I) ** 2, y**2))
    )
    assert reduce(lambda a, b: a.intersection(b), factors) == changed
    parameterized = PolynomialRing((x, y), sp.QQ.algebraic_field(sp.sqrt(-3)).frac_field(t))
    ideal = parameterized.ideal((x - t) ** 2, (x - t) * y)
    same_ideals(
        primary_decomposition(ideal),
        (parameterized.ideal(x - t), parameterized.ideal((x - t) ** 2, y)),
    )


def test_decomposition_limits():
    ring = PolynomialRing((x, y, z))
    with pytest.raises(DecompositionLimitError, match="max_nodes"):
        primary_decomposition(ring.ideal(x**2, x * y), max_nodes=1)
    with pytest.raises(DecompositionLimitError, match="max_dimension"):
        radical(ring.ideal(x * x, y * y, z * z), max_dimension=4)
    with pytest.raises(DecompositionLimitError, match="max_exponent"):
        _saturation_exponent(ring.ideal(x**4, x * y**3), y, ring.ideal(x), 1)
    with pytest.raises(ValueError, match="max_nodes"):
        radical(ring.ideal(), max_nodes=0)
    assert radical(ring.ideal(1)).is_unit
    assert primary_decomposition(ring.ideal(1)) == ()
    assert minimal_primes(ring.ideal()) == (ring.ideal(),)
    assert not is_primary(ring.ideal(1))
    scalar = PolynomialRing(())
    assert is_prime(scalar.ideal())
    with pytest.raises(TypeError, match="Ideal"):
        is_primary(x)


def test_coordinate_changes():
    ring = PolynomialRing((x, y, z), sp.QQ.algebraic_field(sp.sqrt(-3)).frac_field(t))
    form = x + y + z
    ideal = ring.ideal(form**2, form * (2 * x + y))
    factors = primary_decomposition(ideal)
    assert len(factors) == 2
    assert factors[0].intersection(factors[1]) == ideal
    same_ideals(associated_primes(ideal), (ring.ideal(form), ring.ideal(form, 2 * x + y)))
