"""Exact module kernels, witnesses, presentation invariants, and differentials."""

import pytest
import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.maps import singular_locus
from alggeom.modules import (
    PresentedModule,
    Submodule,
    differential,
    differentials,
    fitting_ideal,
    syzygies,
)

x, y, t = sp.symbols("x y t")


@pytest.mark.parametrize("field", [sp.QQ, sp.QQ.algebraic_field(sp.sqrt(-3)), sp.QQ.frac_field(t)])
def test_kernel(field):
    ring = PolynomialRing((x, y), field)
    matrix = sp.ImmutableMatrix([[x, y, 0], [0, 0, x]])
    kernel = syzygies(matrix, ring)
    assert matrix * kernel == sp.zeros(2, kernel.cols)
    module = Submodule(ring, kernel)
    assert module.contains([-y, x, 0])
    assert not module.contains([1, 0, 0])
    target = matrix * sp.Matrix([x + 1, y, 2])
    columns = Submodule(ring, matrix)
    witness = columns.coordinates(target)
    assert sp.expand(matrix * sp.Matrix(witness) - target) == sp.zeros(2, 1)
    assert columns.coordinates([1, 0]) is None
    assert columns.normal_form(target) == sp.zeros(2, 1)


def test_quotient_kernel():
    ring = PolynomialRing((x,))
    coordinate = ring.quotient(x * x)
    kernel = syzygies([[x]], coordinate)
    assert Submodule(coordinate, kernel).contains([x])
    assert not Submodule(coordinate, kernel).contains([1])
    assert kernel.applyfunc(lambda value: coordinate.normal_form(x * value)) == sp.zeros(
        1, kernel.cols
    )
    assert syzygies([[x]], ring).shape == (1, 0)
    columns = Submodule(coordinate, [[x]])
    witness = columns.coordinates([x + x * x])
    assert coordinate.normal_form(x * witness[0] - x) == 0
    assert columns.normal_form([1]) == sp.Matrix([1])


def test_fitting():
    ring = PolynomialRing((x, y))
    matrix = sp.Matrix.diag(x, y)
    module = PresentedModule(ring, matrix)
    assert fitting_ideal(module, 0) == ring.ideal(x * y)
    assert fitting_ideal(module, 1) == ring.ideal(x, y)
    assert fitting_ideal(module, 2) == ring.ideal(1)
    changed = PresentedModule(
        ring, sp.Matrix([[1, x], [0, 1]]) * matrix * sp.Matrix([[1, 0], [y, 1]])
    )
    assert all(changed.fitting_ideal(i) == module.fitting_ideal(i) for i in range(3))
    coordinate = ring.quotient(x * x)
    assert PresentedModule(coordinate, [[y]]).fitting_ideal(0) == ring.ideal(x * x, y)
    assert PresentedModule(ring, sp.eye(2)).is_zero
    assert not module.is_zero


def test_differentials():
    ring = PolynomialRing((x, y))
    curve = ring.quotient(y * y - x**3)
    omega = differentials(curve)
    assert omega.presentation == sp.Matrix([-3 * x * x, 2 * y])
    assert omega.fitting_ideal(1) == singular_locus(curve.ideal)
    assert omega.normal_form(differential(y * y, curve) - differential(x**3, curve)) == sp.zeros(
        2, 1
    )
    f, g = x + y, x * y
    assert omega.normal_form(
        differential(f * g, curve) - f * differential(g, curve) - g * differential(f, curve)
    ) == sp.zeros(2, 1)
    assert differentials(ring).presentation.shape == (2, 0)
    parametric = PolynomialRing((x,), sp.QQ.frac_field(t)).quotient()
    assert differential(t * x, parametric) == sp.Matrix([t])


def test_empty_presentations():
    ring = PolynomialRing((x,))
    assert syzygies(sp.zeros(0, 2), ring) == sp.eye(2)
    assert syzygies(sp.zeros(2, 0), ring).shape == (0, 0)
    assert PresentedModule(ring, sp.zeros(0, 2)).is_zero
    free = PresentedModule(ring, sp.zeros(2, 0))
    assert free.fitting_ideal(1).is_zero and free.fitting_ideal(2).is_unit
    assert free.normal_form([x, 1]) == sp.Matrix([x, 1])
    scalar = PolynomialRing(())
    assert Submodule(scalar, [[2], [3]]).coordinates([4, 6]) == (2,)
    assert Submodule(scalar, syzygies([[2, 3]], scalar)).contains([3, -2])
    zero = ring.quotient(1)
    assert PresentedModule(zero, sp.zeros(2, 0)).is_zero
    assert syzygies([[0]], zero) == sp.zeros(1, 0)
    with pytest.raises(ValueError):
        free.normal_form([1])
    with pytest.raises(ValueError):
        free.fitting_ideal(-1)


def test_algebraic_coefficients():
    alpha = sp.sqrt(-3)
    field = sp.QQ.algebraic_field(alpha).frac_field(t)
    ring = PolynomialRing((x, y), field)
    matrix = sp.Matrix([[x - alpha * y, y - t]])
    kernel = syzygies(matrix, ring)
    assert (matrix * kernel).applyfunc(ring.convert) == sp.zeros(1, kernel.cols)
    assert Submodule(ring, kernel).contains([y - t, -x + alpha * y])
    assert PresentedModule(ring, matrix).fitting_ideal(0) == ring.ideal(x - alpha * y, y - t)
