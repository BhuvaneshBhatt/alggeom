"""Recover a space curve, close it projectively, and inspect a cusp's differentials."""

import sympy as sp

from alggeom.graded import degree, hilbert_polynomial, projective_closure, projective_dimension
from alggeom.ideals import CoordinateRing
from alggeom.maps import implicitize, singular_locus, tangent_space
from alggeom.modules import differentials

t, x, y, z, h, n = sp.symbols("t x y z h n")
curve = implicitize((t, t**2, t**3), t, targets=(x, y, z))
assert curve == curve.ring.ideal(y - x**2, z - x**3)
closure = projective_closure(curve, homogenizer=h)
assert projective_dimension(closure) == 1 and degree(closure) == 3
assert hilbert_polynomial(closure, parameter=n) == sp.Poly(3 * n + 1, n, domain=sp.QQ)

cusp = implicitize((t**2, t**3), t, targets=(x, y))
omega = differentials(CoordinateRing(cusp))
assert tangent_space(cusp, (0, 0)) == sp.eye(2)
assert omega.fitting_ideal(1) == singular_locus(cusp)
print("Space curve degree:", degree(closure))
print("Projective Hilbert polynomial:", hilbert_polynomial(closure, parameter=n).as_expr())
print("Cusp singular subscheme:", singular_locus(cusp).groebner_basis)
