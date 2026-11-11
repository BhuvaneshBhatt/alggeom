"""Rank-certified nonlinear polynomial fiber partitions."""

import sympy as sp

from alggeom.families import FiberInvariants, fiber_stratification


def assert_unique(family, samples):
    for sample in samples:
        matching = [p for p in family.pieces if p.locus.contains(sample)]
        assert len(matching) == 1, (sample, matching)
        assert matching[0].invariants == family.specialize(sample)


def test_multiple_parameter_quadratic():
    x, a, b = sp.symbols("x a b")
    family = fiber_stratification((a * x * x + b * x + 1,), (x,), (a, b))
    assert family.specialize((1, 0)) == FiberInvariants(0, 2)
    assert family.specialize((1, 2)) == FiberInvariants(0, 1)
    assert family.specialize((0, 2)) == FiberInvariants(0, 1)
    assert family.specialize((0, 0)) == FiberInvariants(-1, 0)
    assert_unique(family, [(1, 0), (1, 2), (0, 2), (0, 0), (2, 3)])


def test_repeated_roots_and_zero_fiber():
    x, a, b = sp.symbols("x a b")
    family = fiber_stratification((a * (x - b) ** 2,), (x,), (a, b))
    assert family.specialize((1, 0)) == FiberInvariants(0, 1)
    assert family.specialize((0, 1)) == FiberInvariants(1, sp.oo)
    assert_unique(family, [(1, 0), (1, 2), (0, 1), (-2, 7)])


def test_multidimensional_hypersurface_dimension():
    x, y, a, b = sp.symbols("x y a b")
    family = fiber_stratification((a * x * x + b * y * y,), (x, y), (a, b))
    assert family.specialize((1, 2)) == FiberInvariants(1, sp.oo)
    assert family.specialize((0, 0)) == FiberInvariants(2, sp.oo)
    assert_unique(family, [(1, 2), (0, 0), (0, 1), (1, 0)])


def test_separable_nonlinear_system():
    x, y, a, b = sp.symbols("x y a b")
    family = fiber_stratification((x * x - a, y * y - b), (x, y), (a, b))
    assert family.specialize((1, 1)) == FiberInvariants(0, 4)
    assert family.specialize((0, 1)) == FiberInvariants(0, 2)
    assert family.specialize((0, 0)) == FiberInvariants(0, 1)
    assert_unique(family, [(0, 0), (0, 1), (1, 0), (1, 1)])
