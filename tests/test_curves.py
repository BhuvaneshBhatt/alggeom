"""Public curve models connect exact curve certificates to geometric maps."""

import pytest
import sympy as sp
from algroots.curves import CurveError

from alggeom.curves import PlaneCurve, genus, normalization, parametrize
from alggeom.ideals import PolynomialRing

x, y, t, u, v = sp.symbols("x y t u v")


def test_cusp():
    curve = PlaneCurve(y * y - x**3, (x, y))
    assert curve.degree == 3 and genus(curve) == 0
    analysis = normalization(curve)
    assert (analysis.arithmetic_genus, analysis.delta_sum) == (1, 1)
    assert analysis.verify()
    mapping = parametrize(curve, t)
    assert mapping.functions == (t * t, t**3)
    assert mapping.image_closure == curve.affine_scheme.ideal
    result = curve.rationality(t)
    assert result.parametrization.inverse == y / x and result.verify()
    morphism = curve.normalization_map(parameters=(u, v))
    assert morphism.coordinates == (u * u * v, u**3, v**3)
    assert morphism.image_closure.ideal == curve.projective_scheme.ideal
    point = morphism(morphism.source.point((1, 0)))
    assert point.coordinates == (0, 1, 0)


def test_rationality_obstructions():
    elliptic = PlaneCurve(y * y - x**3 + x, (x, y))
    assert elliptic.genus == 1
    assert elliptic.rationality(t).status == "positive_genus"
    with pytest.raises(ValueError, match="positive_genus"):
        elliptic.parametrize(t)
    conic = PlaneCurve(x * x + y * y - 3, (x, y))
    assert conic.genus == 0
    assert conic.rationality(t).status == "no_ground_field_point"
    with pytest.raises(ValueError, match="no_ground_field_point"):
        conic.normalization_map()
    extended = conic.parametrize(t, extend_constants=True)
    assert extended.source.domain == sp.QQ.algebraic_field(sp.sqrt(3))
    equation = conic.equation.subs(
        dict(zip((x, y), extended.functions, strict=True)), simultaneous=True
    )
    assert sp.cancel(equation) == 0
    assert conic.normalization_map(extend_constants=True).target.domain == extended.source.domain


@pytest.mark.parametrize("alpha", [sp.sqrt(2), sp.sqrt(-3), sp.real_root(2, 3)])
def test_number_fields(alpha):
    field = sp.QQ.algebraic_field(alpha)
    curve = PlaneCurve(y - x * x + alpha * x, (x, y), field)
    assert curve.genus == 0
    result = curve.rationality(t, point=(0, 0))
    assert result.rational and result.verify()
    mapping = curve.parametrize(t, point=(0, 0))
    assert mapping.image_closure == curve.affine_scheme.ideal
    assert mapping.source.domain == field


def test_supported_models():
    with pytest.raises(CurveError, match="reduced and irreducible"):
        PlaneCurve((y - x * x) ** 2, (x, y))
    with pytest.raises(CurveError, match="geometrically integral"):
        genus(PlaneCurve(x * x + 1, (x, y)))
    with pytest.raises(CurveError, match="constants"):
        PlaneCurve(y - x * x, (x, y), sp.QQ.frac_field(t))
    curve = PlaneCurve(y * y - x**3, (x, y))
    assert curve.affine_scheme.ambient == PolynomialRing((x, y))
    with pytest.raises(ValueError, match="two"):
        curve.normalization_map(parameters=(u,))
