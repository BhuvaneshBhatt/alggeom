"""Scheme structure, field-valued points, projective charts and exact fibers."""

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, z, t, s, u, v = sp.symbols("x y z t s u v")


def test_points_and_fields():
    scheme = AffineScheme(PolynomialRing((x,)).ideal(x * x + 1))
    assert not scheme.is_empty and scheme.dimension == 0
    assert scheme.rational_points() == ()
    field = sp.QQ.algebraic_field(sp.I)
    points = scheme.rational_points(field=field)
    assert {p.coordinates for p in points} == {(-sp.I,), (sp.I,)}
    assert all(p.scheme == scheme and p.field == field and p.evaluate(x * x) == -1 for p in points)
    with pytest.raises(ValueError, match="equations"):
        scheme.point((0,))
    assert AffineScheme(PolynomialRing(())).rational_points()[0].coordinates == ()


def test_closed_operations():
    ring = PolynomialRing((x, y))
    first = AffineScheme(ring.ideal(x * x))
    second = AffineScheme(ring.ideal(y))
    intersection = first.intersection(second)
    assert intersection.ideal == ring.ideal(x * x, y)
    assert intersection.coordinate_ring.vector_space_dimension == 2
    assert first.union(second).ideal == ring.ideal(x * x * y)
    assert first.principal_open(x).is_empty
    with pytest.raises(ValueError, match="zero-dimensional"):
        first.rational_points()


def test_affine_fibers():
    source = AffineScheme(PolynomialRing((t,)))
    target = AffineScheme(PolynomialRing((x,)))
    morphism = source.morphism(target, (t * t,))
    fiber = morphism.fiber(target.point((0,)))
    assert fiber.ideal == source.ambient.ideal(t * t)
    assert fiber.coordinate_ring.vector_space_dimension == 2
    assert fiber.rational_points()[0].coordinates == (0,)
    assert morphism(source.point((2,))).coordinates == (4,)
    assert morphism.image_closure == target
    field = sp.QQ.algebraic_field(sp.sqrt(2))
    extended = morphism.fiber(target.point((2,), field=field))
    assert {p.coordinates for p in extended.rational_points()} == {(sp.sqrt(2),), (-sp.sqrt(2),)}
    point = target.subscheme(x - 1)
    assert morphism.preimage(point).ideal == source.ambient.ideal(t * t - 1)
    product = morphism.fiber_product(morphism, variables=(u, v))
    assert product.ideal == product.ambient.ideal(u * u - v * v)


def test_affine_isomorphism():
    scheme = AffineScheme(PolynomialRing((x, y)))
    morphism = scheme.morphism(scheme, (x + y, y))
    assert morphism.inverse().coordinates == (x - y, y)
    assert morphism.compose(morphism.inverse())(scheme.point((3, 4))) == scheme.point((3, 4))
    with pytest.raises(ValueError, match="relations"):
        AffineScheme(PolynomialRing((x,))).morphism(
            AffineScheme(PolynomialRing((u, v)).ideal(u * v - 1)), (x, x)
        )


def test_projective_charts():
    curve = ProjectiveScheme(PolynomialRing((x, y, z)).ideal(x * z - y * y))
    assert curve.dimension == 1 and curve.degree == 2
    first, last = curve.affine_chart(0), curve.affine_chart(2)
    assert first.scheme.ideal == first.scheme.ambient.ideal(z - y * y)
    assert first.transition_coordinates(last) == (1 / z, y / z)
    assert first.overlap_guard(last) == sp.Ne(z, 0)
    point = first.scheme.point((2, 4))
    projective = first.to_projective(point)
    assert projective == curve.point((2, 4, 8))
    assert last.from_projective(projective).coordinates == (sp.Rational(1, 4), sp.Rational(1, 2))
    assert first.from_projective(projective) == point
    assert first.embedding.image_closure == curve
    with pytest.raises(ValueError, match="all be zero"):
        curve.point((0, 0, 0))
    with pytest.raises(ValueError, match="chart"):
        first.from_projective(curve.point((0, 0, 1)))


def test_projective_maps():
    line = ProjectiveScheme(PolynomialRing((s, t)).ideal())
    plane = ProjectiveScheme(PolynomialRing((x, y, z)).ideal())
    veronese = line.morphism(plane, (s * s, s * t, t * t))
    assert veronese.image_closure.ideal == plane.ambient.ideal(x * z - y * y)
    assert veronese(line.point((1, 2))) == plane.point((1, 2, 4))
    cover = line.morphism(line, (s * s, t * t))
    fiber = cover.fiber(line.point((1, 0)))
    assert fiber.ideal == line.ambient.ideal(t * t) and fiber.degree == 2
    assert tuple(p.coordinates for p in fiber.rational_points()) == ((1, 0),)
    assert cover.compose(cover).coordinates == (s**4, t**4)
    with pytest.raises(ValueError, match="base locus"):
        line.morphism(line, (s * s, s * t))
    with pytest.raises(ValueError, match="same degree"):
        line.morphism(line, (s * s, t))


def test_projective_saturation():
    ring = PolynomialRing((x, y))
    assert ProjectiveScheme(ring.ideal(x * x, x * y)).ideal == ring.ideal(x)
    assert ProjectiveScheme(ring.ideal(x * x, y * y)).is_empty
    with pytest.raises(ValueError, match="homogeneous"):
        ProjectiveScheme(ring.ideal(x - 1))
