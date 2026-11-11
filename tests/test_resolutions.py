from dataclasses import replace

import pytest
import sympy as sp

from alggeom.graded import hilbert_series
from alggeom.ideals import PolynomialRing
from alggeom.modules import PresentedModule
from alggeom.resolutions import free_resolution

x, y, z, w, t = sp.symbols("x y z w t")


@pytest.mark.parametrize(
    "powers,shifts", [((1, 1), ((0,), (1, 1), (2,))), ((2, 3), ((0,), (2, 3), (5,)))]
)
def test_complete_intersections(powers, shifts):
    ring = PolynomialRing((x, y))
    ideal = ring.ideal(x ** powers[0], y ** powers[1])
    resolution = free_resolution(ideal)
    assert resolution.shifts == shifts
    assert resolution.complete and resolution.verify()
    assert resolution.projective_dimension == 2
    assert resolution.regularity == sum(powers) - 2
    assert (
        sp.cancel(resolution.hilbert_series(parameter=t) - hilbert_series(ideal, parameter=t)) == 0
    )
    assert resolution.betti_numbers == {
        (i, j): shifts[i].count(j) for i in range(3) for j in shifts[i]
    }
    assert all(
        matrix.rows == len(shifts[i]) and matrix.cols == len(shifts[i + 1])
        for i, matrix in enumerate(resolution.differentials)
    )


def test_twisted_cubic():
    ring = PolynomialRing((x, y, z, w))
    ideal = ring.ideal(x * z - y**2, x * w - y * z, y * w - z**2)
    resolution = free_resolution(ideal)
    assert resolution.shifts == ((0,), (2, 2, 2), (3, 3))
    assert resolution.betti_table == {0: (1, 0, 0), 1: (0, 3, 2)}
    assert resolution.projective_dimension == 2
    assert resolution.regularity == 1
    assert resolution.verify()
    assert sp.cancel(resolution.hilbert_series(parameter=t) - (1 + 2 * t) / (1 - t) ** 2) == 0


def test_presentations():
    ring = PolynomialRing((x, y))
    module = PresentedModule(ring, [[1, x, 0, 0], [x, x * x, y, y * y]])
    # This presentation is homogeneous for shifts (1,0), and represents S/(y).
    with pytest.raises(ValueError, match="homogeneous"):
        module.free_resolution()
    resolution = module.free_resolution(generator_degrees=(1, 0))
    assert resolution.shifts == ((0,), (1,))
    assert resolution.verify()
    assert resolution.augmentation == sp.Matrix([0, 1])
    free = free_resolution(PresentedModule(ring, sp.zeros(2, 0)), generator_degrees=(-2, 3))
    assert free.shifts == ((-2, 3),)
    assert free.projective_dimension == 0
    assert free.regularity == 3
    assert free.verify()
    zero = free_resolution(PresentedModule(ring, sp.eye(2)))
    assert zero.shifts == ((),)
    assert zero.projective_dimension == -1
    assert zero.regularity == -sp.oo
    assert zero.hilbert_series() == 0
    assert zero.verify()


def test_quotient_scalars():
    ring = PolynomialRing((x,))
    module = PresentedModule(ring.quotient(x * x), sp.zeros(1, 0))
    with pytest.raises(ValueError, match="restrict_scalars"):
        free_resolution(module)
    resolution = free_resolution(module.restrict_scalars())
    assert resolution.shifts == ((0,), (2,))
    assert resolution.verify()
    assert resolution.hilbert_series(parameter=t) == t + 1
    scalar = PolynomialRing(())
    result = free_resolution(PresentedModule(scalar, [[2], [3]]))
    assert result.shifts == ((0,),)
    assert result.verify()


def test_truncation_and_mutations():
    ring = PolynomialRing((x, y))
    ideal = ring.ideal(x, y)
    truncated = free_resolution(ideal, max_length=1)
    assert truncated.shifts == ((0,), (1, 1))
    assert not truncated.complete and truncated.verify()
    for invariant in (
        lambda: truncated.projective_dimension,
        lambda: truncated.regularity,
        truncated.hilbert_series,
    ):
        with pytest.raises(ValueError, match="complete"):
            invariant()
    assert not replace(truncated, complete=True).verify()
    result = free_resolution(ideal)
    assert not replace(result, differentials=(result.differentials[0], sp.Matrix([x, y]))).verify()
    assert not replace(result, augmentation=sp.zeros(1, 1)).verify()
    assert not replace(result, differentials=(sp.zeros(1, 2), result.differentials[1])).verify()
    assert not replace(result, shifts=((1,), (2, 2), (3,))).verify()
    assert not replace(
        result, differentials=(result.differentials[0], sp.Matrix([y + x * x, -x]))
    ).verify()
    assert free_resolution(ideal, max_length=0).verify()
    with pytest.raises(ValueError, match="degree"):
        free_resolution(ideal, generator_degrees=(sp.Rational(1, 2),))


@pytest.mark.parametrize(
    "field", [sp.QQ.algebraic_field(sp.sqrt(-3)), sp.QQ.algebraic_field(sp.I).frac_field(t)]
)
def test_coefficient_fields(field):
    ring = PolynomialRing((x, y), field)
    scalar = field.to_sympy(field.one + field.one)
    resolution = free_resolution(ring.ideal(x * x, scalar * x * y, y * y))
    assert resolution.shifts == ((0,), (2, 2, 2), (3, 3))
    assert resolution.verify()
