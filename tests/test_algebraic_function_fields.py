"""Ideals and quotient actions over parameter-dependent algebraic fields."""

import pytest
import sympy as sp
from algroots.native_f4 import compiled_available

from alggeom.fields import AlgebraicFunctionField
from alggeom.ideals import PolynomialRing

x, y, t, a = sp.symbols("x y t a")


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_field_ideal(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    K = AlgebraicFunctionField(sp.Poly(a * a - t, a, domain=sp.QQ.frac_field(t)))
    ring = PolynomialRing((x, y), K)
    ideal = ring.ideal(x * x - a, y - x, groebner_backend=backend)
    assert ideal.normal_form(x * x - a) == 0
    assert ideal.normal_form(x - y) == 0
    assert ring.convert(a * a - t) == 0


def test_field_coordinate_collision():
    K = AlgebraicFunctionField(sp.Poly(a * a - t, a, domain=sp.QQ.frac_field(t)))
    for variable in (t, a):
        with pytest.raises(ValueError, match="distinct"):
            PolynomialRing((variable,), K)
