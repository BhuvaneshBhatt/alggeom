"""Compare finite coordinates with adding integral functions on curves and a surface."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.noether import noether_normalization
from alggeom.normalization import integral_closure

x, y, z = sp.symbols("x y z")
S = PolynomialRing((x, y))
projection = noether_normalization(S.ideal(x * y - 1))
assert projection.verify() and projection.generic_rank == 2

for equation in (y**2 - x**3, y**2 - x**3 * (x * x + 1)):
    source = S.quotient(equation)
    normalized = integral_closure(source)
    assert normalized.verify() and normalized.mapping.is_injective
    assert normalized.ring.is_normal and not source.is_normal

surface = PolynomialRing((x, y, z)).quotient(x * x - y * y * z)
normalized = integral_closure(surface)
assert normalized.verify() and normalized.ring.krull_dimension == 2
assert normalized.conductor == surface.ambient.ideal(x, y)
print(
    "A finite projection has degree 2; curve and surface integral closures add verified finite birational coordinates."
)
