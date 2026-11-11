import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.rational_functions import FractionField, local_ring

x, y, z = sp.symbols("x y z")


def test_quotient_function_field():
    ring = PolynomialRing((x, y)).quotient(y * y - x**3)
    field = ring.fraction_field()
    parameter = field(y / x)
    assert parameter**2 == field(x) and parameter**3 == field(y)
    assert parameter.inverse() * parameter == 1
    assert field((y * y) / x) == field(x * x)
    assert field(x).to_coordinate() == ring(x)
    assert ring(x) == field(x) and field(x) == ring(x)
    assert ring(x) + parameter == parameter + field(x)
    assert ring(x) * parameter == field(y)
    assert ring(x) - parameter == field(x) - parameter
    with pytest.raises(ValueError, match="base coordinate ring"):
        parameter.to_coordinate()
    with pytest.raises(ZeroDivisionError):
        field.fraction(1, y * y - x**3)
    with pytest.raises(ValueError, match="integral"):
        FractionField(PolynomialRing((x, y)).quotient(x * y))


def test_principal_localization():
    ring = PolynomialRing((x, y)).quotient(x * y)
    localized = ring.localize(x)
    assert localized(y) == 0 and localized(x).is_unit
    assert localized(1 / x) * localized(x) == 1
    assert localized.natural_map(y) == 0
    assert localized.kernel == ring.ambient.ideal(y)
    zero = ring.localize(x * y)
    assert zero.zero == zero.one


def test_local_regularity():
    ring = PolynomialRing((x, y, z)).quotient(x * y - z)
    local = local_ring(ring, (0, 0, 0))
    assert local(z / x) == local(y)
    assert local(1 / (1 + x)).is_unit
    assert not local(x).is_unit
    with pytest.raises(ValueError, match="not regular"):
        local(1 / x)
    cusp = PolynomialRing((x, y)).quotient(y * y - x**3)
    with pytest.raises(ValueError, match="not regular"):
        local_ring(cusp, (0, 0))(y / x)


def test_local_component_removal():
    ring = PolynomialRing((x, y)).quotient(x * y)
    local = local_ring(ring, (1, 0))
    assert local.kernel == ring.ambient.ideal(y)
    assert local(y) == 0 and local(x).is_unit
    node = local_ring(ring, (0, 0))
    with pytest.raises(ZeroDivisionError, match="zero-divisor"):
        node.fraction(x, x)


def test_fraction_parent_coercion():
    ambient = PolynomialRing((x,))
    base = ambient.quotient()
    field = base.fraction_field()
    unrelated = ambient.quotient(x * x - x)
    idempotent = unrelated(x)
    assert idempotent**2 == idempotent
    for operation in (
        lambda: field(idempotent),
        lambda: field.fraction(idempotent, x),
        lambda: field.fraction(1, idempotent),
        lambda: field.one + idempotent,
        lambda: idempotent + field.one,
        lambda: field.one * idempotent,
    ):
        with pytest.raises(TypeError, match="base coordinate ring"):
            operation()
    assert field(base(x)) == field(x)
    assert field.fraction(field(x), base(x)) == 1
    other_field = ambient.quotient(x * x - 2).fraction_field()
    with pytest.raises(TypeError, match="base coordinate ring"):
        field.fraction(1, other_field(x))
    ordered = ambient.with_order("lex").quotient()
    assert field(ordered(x)) ** 2 == field(x * x)
