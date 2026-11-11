from dataclasses import replace

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.normalization import integral_closure, is_normal
from alggeom.schemes import AffineScheme

x, y, z = sp.symbols("x y z")


@pytest.mark.parametrize("equation", [y * y - x**3, y * y - x * x * (x + 1), y * y - x**5])
def test_curve_integral_closure(equation):
    ring = PolynomialRing((x, y))
    original = ring.quotient(equation)
    data = integral_closure(original)
    assert data.verify() and data.mapping.is_injective
    assert data.ring.is_normal and not original.is_normal
    assert data.conductor.radical() == ring.ideal(x, y)
    assert data.module_generators[0] == 1
    bad = replace(data, steps=())
    assert not bad.verify()
    assert not replace(data, module_generators=(data.module_generators[0],)).verify()
    broken_step = replace(data.steps[0], denominator=sp.S.Zero)
    assert not broken_step.verify()
    assert not replace(data, steps=(broken_step, *data.steps[1:])).verify()
    scheme = AffineScheme(original)
    morphism = scheme.normalization_map()
    assert morphism.target == scheme
    assert morphism.source.is_normal


def test_positive_genus_and_surfaces():
    curve = PolynomialRing((x, y)).quotient(y * y - x**3 * (x * x + 1))
    normalized = integral_closure(curve)
    assert normalized.verify() and normalized.ring.is_normal
    surface = PolynomialRing((x, y, z)).quotient(x * x - y * y * z)
    data = integral_closure(surface)
    assert data.verify() and data.ring.krull_dimension == 2
    assert data.conductor == surface.ambient.ideal(x, y)
    normal = PolynomialRing((x, y, z)).quotient(x * y - z * z)
    assert is_normal(normal) and integral_closure(normal).ring == normal


def test_reducible_and_empty():
    ring = PolynomialRing((x, y))
    node = ring.quotient(x * y)
    result = integral_closure(node)
    assert len(result.components) == 2 and result.verify()
    assert result.mapping.is_injective
    assert not replace(result, components=result.components[:1]).verify()
    with pytest.raises(RuntimeError, match="max_components"):
        integral_closure(node, max_components=1)
    assert not node.is_normal
    assert ring.quotient(x * (x - 1), y).is_normal
    empty = integral_closure(ring.quotient(1))
    assert empty.ring.ideal.is_unit and empty.verify()


def test_fields_and_limits():
    ring = PolynomialRing((x, y), sp.QQ.algebraic_field(sp.sqrt(-3)))
    data = integral_closure(ring.quotient(y * y - x**3))
    assert data.verify() and data.ring.domain == ring.domain
    with pytest.raises(RuntimeError, match="max_steps"):
        integral_closure(ring.quotient(y * y - x**3), max_steps=0)
    with pytest.raises(RuntimeError, match="max_module_generators"):
        integral_closure(ring.quotient(y * y - x**3), max_module_generators=1)
    with pytest.raises(ValueError, match="reduced"):
        integral_closure(ring.quotient(x * x))


def test_step_regular_denominator():
    from alggeom.normalization import NormalizationStep

    ambient = PolynomialRing((x, y))
    source = ambient.quotient(x * y)
    removed_component = NormalizationStep(source, ambient.quotient(y), ambient.ideal(x, y), x, ())
    assert source(x) != 0 and source(x) * source(y) == 0
    assert not removed_component.verify()


def test_step_constants():
    step = integral_closure(PolynomialRing((x, y)).quotient(y * y - x**3)).steps[0]
    extended = step.target.base_change(sp.QQ.frac_field(z))
    assert not replace(step, target=extended).verify()
