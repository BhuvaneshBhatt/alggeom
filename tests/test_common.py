"""Common intermediate fields and exact simultanous outer maps."""

from dataclasses import replace

import pytest
import sympy as sp

from alggeom.decomposition import common_right_component
from alggeom.decomposition.common import common_components

x, y, t = sp.symbols("x y t")


@pytest.mark.parametrize("alpha", [0, sp.sqrt(2), sp.sqrt(-3)])
def test_shared_generator(alpha):
    constants = sp.QQ if alpha == 0 else sp.QQ.algebraic_field(alpha)
    h = x**2 + alpha * x + 1
    functions = (h**2, h**3)
    result = common_components(functions, x, constants=constants, parameter=t)
    assert result.status == "rational" and result.verify()
    assert all(
        sp.Poly(
            sp.together(outer(result.inner(x)) - f).as_numer_denom()[0], x, domain=constants
        ).is_zero
        for outer, f in zip(result.outers, functions, strict=True)
    )
    assert common_right_component(functions, x, constants=constants)(x) == result.inner(x)
    assert not replace(result, functions=(functions[0] + 1, functions[1])).verify()


def test_multivariate():
    h = x * x + y * y
    result = common_components([h * h, h**3], (x, y), parameter=t)
    assert result.verify() and result.status == "rational"
    assert all(
        sp.cancel(g(result.inner(x, y)) - f) == 0
        for g, f in zip(result.outers, (h * h, h**3), strict=True)
    )
    independent = common_components([x, y], (x, y))
    assert independent.status == "higher_transcendence_degree" and independent.verify()
    assert independent.inner is None and independent.outers == ()


def test_constant_and_primitive():
    result = common_components([2, 3], x)
    assert result.status == "constant" and result.verify() and result.outers == ()
    assert common_right_component([2, 3], x) is None
    result = common_components([x * x, x**3], x)
    assert result.verify() and sp.degree(result.inner(x), x) == 1
    with pytest.raises(ValueError, match="at least one"):
        common_components([], x)


def test_nonmonic_poles():
    alpha = sp.sqrt(2)
    constants = sp.QQ.algebraic_field(alpha)
    h = (alpha * x + y) / (x - 2 * y)
    functions = (h * h, h**3, alpha)
    result = common_components(functions, (x, y), constants=constants, parameter=t)
    assert result.verify() and result.status == "rational"
    assert result.input_guard.subs(x, 2 * y) is sp.false
    assert all(
        sp.cancel(g(result.inner(x, y)) - f) == 0
        for g, f in zip(result.outers, functions, strict=True)
    )
