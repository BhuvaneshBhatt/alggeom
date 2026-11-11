"""Generated coordinate changes with independently prescribed normalized models."""

import pytest
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from alggeom.ideals import PolynomialRing
from alggeom.local import local_dimension, multiplicity, tangent_cone
from alggeom.noether import noether_normalization
from alggeom.normalization import integral_closure
from alggeom.rational_functions import Localization, LocalRing
from alggeom.schemes import ProjectiveScheme

x, y, z, t = sp.symbols("x y z t")


@given(st.integers(-2, 2))
@settings(max_examples=5, deadline=None)
def test_transported_cusp(shear):
    S = PolynomialRing((x, y))
    A = S.quotient((y - shear * x) ** 2 - x**3)
    assert local_dimension(A, (0, 0)) == 1 and multiplicity(A, (0, 0)) == 2
    assert tangent_cone(A, (0, 0)) == S.ideal((y - shear * x) ** 2)
    B = integral_closure(A)
    T = PolynomialRing((t,)).quotient()
    parametrization = A.hom((t * t, t**3 + shear * t * t), T)
    substitute = {x: t * t, y: t**3 + shear * t * t}
    images = tuple(
        sp.cancel(f.as_expr().subs(substitute, simultaneous=True)) for f in B.inverse_functions
    )
    comparison = B.ring.hom(images, T)
    assert comparison.is_injective and comparison.is_surjective
    assert comparison.compose(B.mapping).images == parametrization.images
    assert B.conductor == S.ideal(x, y)
    projected = noether_normalization(A)
    base_coordinate = projected.inverse_map(projected.base_embedding.images[0])
    polynomial = parametrization(base_coordinate).as_expr()
    assert projected.generic_rank == sp.degree(polynomial, t)


def test_nonmaximal_prime():
    S = PolynomialRing((x, y))
    A = S.quotient()
    local = LocalRing(A, S.ideal(x))
    assert local(y).is_unit and not local(x).is_unit
    assert local((x * y) / y) == local(x)
    assert local.residue_field.base.ideal == S.ideal(x)
    with pytest.raises(ValueError, match="regular"):
        local(1 / x)


def test_localization_budget():
    A = PolynomialRing((x,)).quotient()
    localized = Localization(A, (x,), max_power=1)
    # x^2 is invertible in this localization; a requested witness for (x^2)
    # exceeds the explicit power budget even though existence is certified.
    with pytest.raises(RuntimeError, match="max_power"):
        localized._unit_in(A.ambient.ideal(x * x))


def test_projective_chart_agreement():
    C = ProjectiveScheme(PolynomialRing((x, y, z)).ideal(y * y * z - x**3))
    point = C.point((1, 1, 1))
    for index in (0, 1, 2):
        chart = C.affine_chart(index)
        p = chart.from_projective(point)
        assert chart.scheme.local_dimension(p) == 1
        assert chart.scheme.embedding_dimension(p) == 1
        assert chart.scheme.multiplicity(p) == 1


@given(st.integers(-3, 3))
@settings(max_examples=7, deadline=None)
def test_transported_components(shear):
    from alggeom.components import primary_decomposition, radical
    from alggeom.graded import degree, hilbert_function, krull_dimension

    S = PolynomialRing((x, y))
    u, v = x + shear * y, y
    ideal = S.ideal(u * u, u * v)
    components = primary_decomposition(ideal)
    assert len(components) == 2
    assert components[0].intersection(components[1]) == ideal
    assert radical(ideal) == S.ideal(u)
    assert any(S.ideal(u) == Q for Q in components)
    assert any(S.ideal(u * u, v) == Q for Q in components)
    assert krull_dimension(ideal) == 1 and degree(ideal) == 1
    assert [hilbert_function(ideal, d) for d in range(5)] == [1, 2, 1, 1, 1]
