"""Independent local, fiber, component and morphism invariants."""

import pytest
import sympy as sp

from alggeom.constructible import ConstructibleSet, LocallyClosedSet
from alggeom.families import FiberInvariants, fiber_stratification
from alggeom.geometric import (
    absolute_primary_decomposition,
    geometric_components,
    is_geometrically_integral,
)
from alggeom.ideals import PolynomialRing
from alggeom.local import intersection_multiplicity
from alggeom.maps import RationalMap
from alggeom.morphism_invariants import (
    branch_locus,
    is_birational,
    is_finite,
    morphism_degree,
    ramification_locus,
)
from alggeom.projective_normalization import projective_normalization_atlas
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, z, a, b, s, t, u, v = sp.symbols("x y z a b s t u v")


def test_absolute_univariate_primary_multiplicity():
    ring = PolynomialRing((x,))
    original = ring.ideal((x**2 + 1) ** 2)
    split = absolute_primary_decomposition(original)
    assert len(split.components) == 2
    assert {component.krull_dimension for component in split.components} == {0}
    expected = original.change_ring(split.components[0].ring)
    assert split.components[0].intersection(split.components[1]) == expected
    assert not is_geometrically_integral(original)


def test_absolute_finite_algebra_components():
    ring = PolynomialRing((x, y))
    original = ring.ideal(x**2 - 1, y**2 - 1)
    split = absolute_primary_decomposition(original)
    assert len(split.components) == 4
    reconstructed = split.components[0].ring.ideal(1)
    for component in split.components:
        reconstructed = reconstructed.intersection(component)
    assert reconstructed == original.change_ring(split.components[0].ring)
    assert len(geometric_components(original).components) == 4
    assert not is_geometrically_integral(original)


def test_nonmonogenic_local_algebra_primary():
    ring = PolynomialRing((x, y))
    original = ring.ideal(x**2, y**2)
    factors = absolute_primary_decomposition(original)
    assert len(factors.components) == 1
    assert factors.components[0] == original
    assert factors.components[0].radical() == ring.ideal(x, y)


def test_line_image_missing_origin():
    source = AffineScheme(PolynomialRing((s, t)).ideal(s * t - 1))
    target = AffineScheme(PolynomialRing((u,)))
    stratified = source.morphism(target, (s,)).constructible_image.as_constructible()
    assert isinstance(stratified, ConstructibleSet)
    assert stratified.contains((1,))
    assert not stratified.contains((0,))
    assert stratified.closure_ideal.is_zero


def test_line_image_closed_and_surjective():
    target = AffineScheme(PolynomialRing((u,)))
    source = AffineScheme(PolynomialRing((s,)))
    assert source.morphism(target, (s**2,)).constructible_image.as_constructible().contains((-1,))
    closed = AffineScheme(PolynomialRing((s,)).ideal(s**2 + 1))
    image = closed.morphism(target, (s,)).constructible_image.as_constructible()
    assert image.closure_ideal == target.ambient.ideal(u**2 + 1)
    assert not image.contains((0,))


def test_higher_dimensional_linear_family():
    family = fiber_stratification((a * x + b * y - 1, a * y), (x, y), (a, b))
    assert family.specialize((1, 0)) == FiberInvariants(0, 1)
    assert family.specialize((0, 1)) == FiberInvariants(1, sp.oo)
    assert family.specialize((0, 0)) == FiberInvariants(-1, 0)
    for value in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        matches = [piece for piece in family.pieces if piece.locus.contains(value)]
        assert len(matches) == 1
        actual = family.specialize(value)
        assert matches[0].dimension == actual.dimension
        assert matches[0].geometric_points == actual.geometric_points


def test_multivariate_polynomial_image_strata():
    source = AffineScheme(PolynomialRing((s, t)))
    target = AffineScheme(PolynomialRing((x, y)))
    image = source.morphism(target, (s, s * t)).constructible_image
    pieces = image.as_constructible()
    assert len(pieces.pieces) == 2
    for point in [(0, 0), (1, 2), (0, 3)]:
        assert pieces.contains(point) == image.contains(point)


def test_affine_plane_morphism_degree_and_branch():
    source = AffineScheme(PolynomialRing((s, t)))
    target = AffineScheme(PolynomialRing((x, y)))
    square = source.morphism(target, (s**2, t))
    assert is_finite(square) and morphism_degree(square) == 2
    assert ramification_locus(square).ideal == source.ambient.ideal(s)
    assert branch_locus(square).ideal == target.ambient.ideal(x)
    open_map = source.morphism(target, (s, s * t))
    assert not is_finite(open_map)
    assert morphism_degree(open_map) == 1 and is_birational(open_map)


def test_projective_cusp_normalization_cocycle():
    curve = ProjectiveScheme(PolynomialRing((x, y, z)).ideal(x * x * z - y**3))
    atlas = projective_normalization_atlas(curve)
    assert len(atlas.patches) == 3
    assert atlas.verify_cocycle()
    transition = atlas.transition(2, 1)
    assert transition.guard == y
    assert len(transition.coordinates) == len(atlas._patch(1).normalized_scheme.variables)


def test_constructible_closure_is_reduced():
    ring = PolynomialRing((x,))
    thickened = LocallyClosedSet(ring.ideal(x**2))
    assert thickened.closure_ideal == ring.ideal(x)
    source = AffineScheme(PolynomialRing((s,)).ideal(s**2))
    target = AffineScheme(PolynomialRing((u,)))
    image = source.morphism(target, (s,)).constructible_image
    assert image.closure_ideal == target.ambient.ideal(u)
    assert image.as_constructible().closure_ideal == target.ambient.ideal(u)


def test_ramification_on_singular_source():
    source = AffineScheme(PolynomialRing((s, t)).ideal(t**2 - s**3))
    target = AffineScheme(PolynomialRing((x,)))
    projection = source.morphism(target, (s,))
    assert is_finite(projection)
    assert ramification_locus(projection).ideal == source.ambient.ideal(t, t**2 - s**3)
    assert branch_locus(projection).ideal == target.ambient.ideal(x**3)


def test_finiteness_over_nonreduced_target():
    source = AffineScheme(PolynomialRing((s, t)).ideal(t**2 - s, s**2))
    target = AffineScheme(PolynomialRing((x,)).ideal(x**2))
    mapping = source.morphism(target, (s,))
    assert is_finite(mapping)


def test_local_component_selection():
    ring = PolynomialRing((x, y))
    first = ring.ideal(y * (x - 1))
    second = ring.ideal((y - x) * (x - 1))
    assert (first + second).krull_dimension == 1
    assert intersection_multiplicity(first, second, (0, 0)) == 1
    with pytest.raises(ValueError, match="proper and finite"):
        intersection_multiplicity(first, second, (1, 0))


def test_geometric_integrality_and_splitting():
    ring = PolynomialRing((x, y))
    assert is_geometrically_integral(ring.ideal(x * x + y * y + 1))
    reducible = ring.ideal(x * x + y * y)
    assert not is_geometrically_integral(reducible)
    decomposition = geometric_components(reducible)
    assert decomposition.field == sp.QQ.algebraic_field(sp.I)
    assert len(decomposition.components) == 2
    assert all(component.krull_dimension == 1 for component in decomposition.components)
    assert decomposition.components[0].intersection(
        decomposition.components[1]
    ) == reducible.change_ring(decomposition.components[0].ring)
    assert is_geometrically_integral(ring.ideal(y * y - x**3))


def test_constructible_boolean_operations():
    ring = PolynomialRing((x, y))
    plane = ConstructibleSet((LocallyClosedSet(ring.ideal()),))
    axis = ConstructibleSet((LocallyClosedSet(ring.ideal(x)),))
    complement = plane - axis
    assert complement.contains((1, 2))
    assert not complement.contains((0, 2))
    assert complement.closure_ideal.is_zero
    point = ConstructibleSet((LocallyClosedSet(ring.ideal(x, y)),))
    union = complement | point
    assert union.contains((0, 0))
    assert not union.contains((0, 1))


def test_exact_image_membership():
    source = AffineScheme(PolynomialRing((s, t)))
    target = AffineScheme(PolynomialRing((u, v)))
    image = source.morphism(target, (s, s * t)).constructible_image
    assert image.closure_ideal.is_zero
    assert image.contains((2, 4))
    assert image.contains((0, 0))
    assert not image.contains((0, 1))
    fraction = RationalMap(PolynomialRing((s,)), (1 / s,), (u,)).constructible_image
    assert fraction.contains((2,))
    assert not fraction.contains((0,))


def test_univariate_family_strata():
    family = fiber_stratification(x * x - a, x, a)
    assert family.generic == FiberInvariants(0, 2)
    assert sp.simplify(family.exceptional_polynomial.subs(a, 0)) == 0
    assert family.specialize(0) == FiberInvariants(0, 1)
    assert family.specialize(9) == FiberInvariants(0, 2)
    assert len(family.strata()) == 2
    assert family.strata()[1].invariants == FiberInvariants(0, 1)
    linear = fiber_stratification(a * x - 1, x, a)
    assert linear.generic == FiberInvariants(0, 1)
    assert linear.specialize(0) == FiberInvariants(-1, 0)
    varying_dimension = fiber_stratification(a * x, x, a)
    assert varying_dimension.specialize(0) == FiberInvariants(1, sp.oo)
    doubled = fiber_stratification((x - a) ** 2, x, a)
    assert doubled.generic == FiberInvariants(0, 1)
    cubic_parameters = fiber_stratification(x**2 - (a**3 - a + 1), x, a)
    with pytest.raises(NotImplementedError, match="root degree"):
        cubic_parameters.strata()


def test_affine_line_morphism_invariants():
    source = AffineScheme(PolynomialRing((t,)))
    target = AffineScheme(PolynomialRing((u,)))
    square = source.morphism(target, (t * t - 1,))
    assert square.is_dominant and square.is_finite
    assert morphism_degree(square) == square.generic_degree == 2
    assert not is_birational(square)
    assert ramification_locus(square).ideal == source.ambient.ideal(t)
    assert branch_locus(square).ideal == target.ambient.ideal(u + 1)
    linear = source.morphism(target, (2 * t + 1,))
    assert is_birational(linear) and linear.ramification_locus().is_empty
    constant = source.morphism(target, (1,))
    assert not is_finite(constant) and not constant.is_dominant
    thick_target = AffineScheme(PolynomialRing((u,)).ideal(u**2))
    reduced_source = AffineScheme(PolynomialRing((t,)).ideal(t))
    assert reduced_source.morphism(thick_target, (0,)).is_dominant
    with pytest.raises(ValueError, match="no generically finite"):
        morphism_degree(constant)


@pytest.mark.slow
def test_projective_normalization_patch_cover():
    from alggeom.projective_normalization import projective_normalization_atlas

    curve = ProjectiveScheme(PolynomialRing((x, y, z)).ideal(x * z - y**2))
    atlas = projective_normalization_atlas(curve)
    assert {patch.chart.index for patch in atlas.patches} == {0, 1, 2}
    assert all(patch.normalized_scheme.dimension == 1 for patch in atlas.patches)
    assert atlas.patches[0].normalized_scheme.ideal.contains(atlas.overlap_guard(0, 2) - y**2)
