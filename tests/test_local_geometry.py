import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.local import (
    embedding_dimension,
    intersection_multiplicity,
    local_dimension,
    module_length,
    multiplicity,
    tangent_cone,
)
from alggeom.maps import singular_locus
from alggeom.modules import PresentedModule
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, z = sp.symbols("x y z")


def test_cusp_local_data():
    ring = PolynomialRing((x, y))
    curve = AffineScheme(ring.ideal(y * y - x**3))
    assert tangent_cone(curve, (0, 0)) == ring.ideal(y * y)
    assert curve.local_dimension((0, 0)) == 1
    assert curve.embedding_dimension((0, 0)) == 2
    assert curve.multiplicity((0, 0)) == 2
    assert curve.is_reduced and curve.is_integral and curve.is_irreducible
    assert curve.codimension == 1 and curve.is_equidimensional and curve.is_unmixed
    assert not curve.is_smooth() and not curve.is_smooth((0, 0))
    assert curve.is_smooth((1, 1))
    assert curve.singular_locus().ideal == ring.ideal(x * x, y)
    assert curve.singular_locus(reduced=True).ideal == ring.ideal(x, y)


def test_initial_ideal_cancellations():
    ring = PolynomialRing((x, y))
    ideal = ring.ideal(x * x + y**3, x * x - y**3)
    assert tangent_cone(ideal, (0, 0)) == ring.ideal(x * x, y**3)
    assert local_dimension(ideal, (0, 0)) == 0
    assert embedding_dimension(ideal, (0, 0)) == 2
    assert multiplicity(ideal, (0, 0)) == 6


def test_mixed_components():
    ring = PolynomialRing((x, y, z))
    mixed = AffineScheme(ring.ideal(x * y, x * z))
    assert not mixed.is_equidimensional
    assert not mixed.is_smooth((0, 0, 0))
    assert singular_locus(mixed, reduced=True).ideal == ring.ideal(x, y, z)
    with pytest.raises(ValueError, match="mixed-dimensional"):
        mixed.singular_locus()
    disjoint = AffineScheme(ring.ideal(x).intersection(ring.ideal(x - 1, y, z)))
    assert not disjoint.is_equidimensional and disjoint.is_smooth()
    embedded = AffineScheme(ring.ideal(x * x, x * y))
    assert embedded.is_equidimensional and not embedded.is_unmixed
    assert not embedded.is_reduced and not embedded.is_integral
    fat = AffineScheme(
        PolynomialRing((x, y)).ideal(x * x).intersection(PolynomialRing((x, y)).ideal(x - 1, y))
    )
    assert not fat.is_smooth()
    assert fat.singular_locus(reduced=True).ideal == fat.ideal.ring.ideal(x)


def test_projective_local_data():
    ring = PolynomialRing((x, y, z))
    curve = ProjectiveScheme(ring.ideal(y * y * z - x**3))
    point = curve.point((0, 0, 1))
    assert curve.local_dimension(point) == 1
    assert curve.embedding_dimension(point) == 2
    assert curve.multiplicity(point) == 2
    assert not curve.is_smooth(point)
    assert not curve.is_normal


def test_intersection_lengths():
    ring = PolynomialRing((x, y))
    cusp = ring.ideal(y * y - x**3)
    assert intersection_multiplicity(cusp, ring.ideal(y), (0, 0)) == 3
    assert intersection_multiplicity(cusp, ring.ideal(x), (0, 0)) == 2
    assert intersection_multiplicity(ring.ideal(y - x * x), ring.ideal(y), (0, 0)) == 2
    with pytest.raises(ValueError, match="proper and finite"):
        intersection_multiplicity(cusp, cusp, (0, 0))
    module = PresentedModule(ring, [[x * (x - 1), y * y]])
    assert module_length(module) == 4
    assert module_length(module, (0, 0)) == 2
    assert module_length(module, (2, 0)) == 0


def test_higher_tor_correction():
    ring = PolynomialRing((x, y, z))
    line = ring.ideal(x, y)
    plane_with_embedded_point = ring.ideal(x * z, y * z, z * z)
    assert (line + plane_with_embedded_point).krull_dimension == 0
    assert ring.quotient(*(line + plane_with_embedded_point).generators).vector_space_dimension == 2
    assert intersection_multiplicity(line, plane_with_embedded_point, (0, 0, 0)) == 1
