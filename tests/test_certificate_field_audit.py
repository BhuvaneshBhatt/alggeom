"""Explicit binding and omission tests for finite projections, normalizations and resolutions."""

import sympy as sp
from _field_audit import audit_fields

from alggeom.ideals import CoordinateRing, PolynomialRing
from alggeom.modules import PresentedModule
from alggeom.noether import noether_normalization
from alggeom.normalization import integral_closure
from alggeom.resolutions import free_resolution

x, y = sp.symbols("x y")


def test_finite_projection_fields():
    S = PolynomialRing((x, y))
    data = noether_normalization(S.ideal(x * y - 1))
    forward = data.coordinate_map.source.hom(
        tuple(-f.as_expr() for f in data.coordinate_map.images), data.ring
    )
    backward = data.ring.hom(
        tuple(-f.as_expr() for f in data.inverse_map.images), data.inverse_map.target
    )
    base = data.base.hom((0,), data.ring)
    report = audit_fields(
        data,
        {
            "coordinate_map": forward,
            "inverse_map": backward,
            "base_embedding": base,
            "integral_relations": (sp.S.Zero,),
            "module_generators": (),
        },
        lambda value: value.verify(),
    )
    assert len(report["fields"]) == 5


def test_normalization_fields():
    S = PolynomialRing((x, y))
    A = S.quotient(y * y - x**3)
    data = integral_closure(A)
    step = data.steps[0]
    point = CoordinateRing(S.ideal(x, y))
    audit_fields(
        step,
        {
            "source": point,
            "target": step.target.ambient.quotient(1),
            "test_ideal": S.ideal(1),
            "denominator": sp.S.Zero,
            "numerators": tuple(sp.S.Zero for c in step.numerators),
        },
        lambda value: value.verify(),
    )
    zero_map = A.hom((0, 0), data.ring)
    point_closure = integral_closure(point)
    audit_fields(
        data,
        {
            "mapping": zero_map,
            "steps": (),
            "inverse_functions": tuple(f.parent.zero for f in data.inverse_functions),
            "module_generators": (),
            "components": (point_closure,),
        },
        lambda value: value.verify(),
    )


def test_resolution_fields():
    S = PolynomialRing((x, y))
    data = free_resolution(S.ideal(x, y), max_length=1)
    assert not data.complete
    audit_fields(
        data,
        {
            "module": PresentedModule(S, [[x]]),
            "generator_degrees": (1,),
            "differentials": tuple(sp.zeros(*d.shape) for d in data.differentials),
            "shifts": (),
            "augmentation": sp.zeros(*data.augmentation.shape),
            "complete": True,
        },
        lambda value: value.verify(),
    )
