"""Locate rank drops of a polynomial system through Fitting ideals."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.modules import PresentedModule, Submodule, syzygies

x, y = sp.symbols("x y")
ring = PolynomialRing((x, y))
matrix = sp.diag(x, y)
model = PresentedModule(ring, matrix)
assert model.fitting_ideal(0) == ring.ideal(x * y)
assert model.fitting_ideal(1) == ring.ideal(x, y)
assert model.fitting_ideal(2).is_unit
# The generic matrix has rank two; one equation disappears on each axis,
# and both disappear at the origin.
assert matrix.subs({x: 1, y: 1}).rank() == 2
assert matrix.subs({x: 0, y: 1}).rank() == 1
assert matrix.subs({x: 0, y: 0}).rank() == 0

constraint = sp.Matrix([[x, y]])
relations = syzygies(constraint, ring)
assert Submodule(ring, relations).contains([-y, x])
print("Rank below two:", model.fitting_ideal(0).groebner_basis)
print("Rank zero:", model.fitting_ideal(1).groebner_basis)
