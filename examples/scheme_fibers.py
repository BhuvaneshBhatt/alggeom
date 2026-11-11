"""Ramification multiplicity survives in scheme-theoretic fibers."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme

t, x = sp.symbols("t x")
line = AffineScheme(PolynomialRing((t,)))
target = AffineScheme(PolynomialRing((x,)))
square = line.morphism(target, (t * t,))
fiber = square.fiber(target.point((0,)))
assert fiber.coordinate_ring.vector_space_dimension == 2
assert fiber.coordinate_ring(t) != 0
assert fiber.coordinate_ring(t) ** 2 == 0
assert len(fiber.rational_points()) == 1
field = sp.QQ.algebraic_field(sp.sqrt(2))
split = square.fiber(target.point((2,), field=field))
assert {point.coordinates for point in split.rational_points()} == {(-sp.sqrt(2),), (sp.sqrt(2),)}
print(
    "Fiber over zero:",
    fiber.ideal.generators,
    "length:",
    fiber.coordinate_ring.vector_space_dimension,
)
print("Fiber over two over QQ(sqrt(2)):", [point.coordinates for point in split.rational_points()])
