"""Enumerate complete field-valued point sets without confusing field order and characteristic."""

import sympy as sp

from alggeom.fields import field_order, finite_field
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, z = sp.symbols("x y z")
for field in (finite_field(3), finite_field(3, 2)):
    q = field_order(field)
    parabola = AffineScheme(PolynomialRing((x, y), field).ideal(y - x**2))
    assert len(parabola.rational_points()) == q
    line = ProjectiveScheme(PolynomialRing((x, y), field))
    assert len(line.rational_points()) == q + 1
    plane = ProjectiveScheme(PolynomialRing((x, y, z), field))
    assert len(plane.rational_points()) == q * q + q + 1
    print(f"GF({q}): {q} points on the affine parabola, {q + 1} on P^1, {q * q + q + 1} on P^2.")
