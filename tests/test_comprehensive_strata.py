"""Exact dimension and point-cardinality strata of coupled polynomial fibers."""

import pytest
import sympy as sp

from alggeom.comprehensive_strata import (
    StratificationLimitError,
    comprehensive_fiber_stratification,
)
from alggeom.families import FiberInvariants


def test_coupled_nonlinear_dimension():
    x, y, a = sp.symbols("x y a")
    family = comprehensive_fiber_stratification((x * y - a, y), (x, y), (a,))
    assert family.specialize((0,)) == FiberInvariants(1, sp.oo)
    assert family.specialize((1,)) == FiberInvariants(-1, 0)
    assert family.specialize((-1,)) == FiberInvariants(-1, 0)


def test_coupled_zero_dimensional_system():
    x, y, a, b = sp.symbols("x y a b")
    family = comprehensive_fiber_stratification((x * y - a, x + y - b), (x, y), (a, b))
    assert family.specialize((0, 0)) == FiberInvariants(0, None)
    assert family.specialize((1, 0)) == FiberInvariants(0, None)
    assert family.specialize((2, 1)) == FiberInvariants(0, None)


def test_geometric_distinct_points():
    x, a = sp.symbols("x a")
    family = comprehensive_fiber_stratification(
        (x * x - a,), (x,), (a,), count_points=True, max_points=2
    )
    assert family.specialize((0,)) == FiberInvariants(0, 1)
    assert family.specialize((1,)) == FiberInvariants(0, 2)
    assert family.specialize((-1,)) == FiberInvariants(0, 2)


def test_dimension_jump_and_inequation():
    x, y, a = sp.symbols("x y a")
    family = comprehensive_fiber_stratification((a * x, y), (x, y), (a,), nonzero=(x,))
    assert family.specialize((0,)) == FiberInvariants(1, sp.oo)
    assert family.specialize((1,)) == FiberInvariants(-1, 0)


def test_point_bound_exhaustion_is_explicit():
    x, a = sp.symbols("x a")
    with pytest.raises(StratificationLimitError):
        comprehensive_fiber_stratification((x**5 - a,), (x,), (a,), count_points=True, max_points=4)


def test_coupled_exact_point_count():
    x, y, a = sp.symbols("x y a")
    family = comprehensive_fiber_stratification(
        (x * x - a, y - x), (x, y), (a,), count_points=True, max_points=4
    )
    assert family.specialize((0,)) == FiberInvariants(0, 1)
    assert family.specialize((1,)) == FiberInvariants(0, 2)
    assert family.specialize((-1,)) == FiberInvariants(0, 2)


def test_coupled_fiber_dispatch():
    from alggeom.families import fiber_stratification

    x, y, a = sp.symbols("x y a")
    family = fiber_stratification((x * y - a, y), (x, y), (a,))
    assert family.specialize((0,)) == FiberInvariants(1, sp.oo)
    assert family.specialize((1,)) == FiberInvariants(-1, 0)


def test_constant_family_without_parameters():
    x = sp.Symbol("x")
    family = comprehensive_fiber_stratification((x * x - 1,), (x,), (), count_points=True)
    assert family.specialize(()) == FiberInvariants(0, 2)


def test_algebraic_specialization():
    x, a = sp.symbols("x a")
    family = comprehensive_fiber_stratification((x * x - a,), (x,), (a,), count_points=True)
    assert family.specialize((sp.sqrt(2),)) == FiberInvariants(0, 2)
