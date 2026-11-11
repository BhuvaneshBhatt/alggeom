"""Discover factor coefficient extensions and verify an implicit composition."""

import sympy as sp

from alggeom.decomposition.discovery import search_algebraic
from alggeom.decomposition.factors import discover_factors

x, z, t = sp.symbols("x z t")
relation = x**4 - z * x * x * (x + 1) + z**3 * (x + 1) ** 2
factors = discover_factors(relation, (x, z), degrees=(2,))
assert factors.verify() and factors.complete
assert sum(orbit.extension_degree for orbit in factors.orbits) == 6
search = search_algebraic(relation, (x, z), parameter=t)
pair = search.decompositions[0]
assert search.verify() and pair.verify()
resultant = sp.Poly(sp.resultant(pair.outer.as_expr(), pair.inner.as_expr(), t), x, z, domain=sp.QQ)
assert resultant == pair.unit * pair.relation**pair.exponent
print("Outer and inner relations:", pair.outer.as_expr(), pair.inner.as_expr())
