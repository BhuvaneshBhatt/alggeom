"""Find the half-response contour of an illustrative radial response model."""

import sympy as sp
from semialg.decision import equivalent, is_satisfiable

from alggeom.decomposition.rational import search_rational

x, y, t = sp.symbols("x y t", real=True)
radius_squared = x**2 + y**2
response = sp.cancel(radius_squared**2 / (1 + radius_squared**2))
measurement = sp.Rational(1, 2)
search = search_rational(response, (x, y), parameter=t)
(pair,) = search.decompositions

# The scalar outer equation has two candidates; test which inner level sets exist.
scalar_candidates = sp.solve(pair.outer - measurement, t)
admissible_levels = tuple(
    value for value in scalar_candidates if is_satisfiable(sp.Eq(pair.inner, value), (x, y))
)
contour = sp.Or(*(sp.Eq(pair.inner, value) for value in admissible_levels))
assert search.complete and pair.verify()
assert equivalent(contour, sp.Eq(radius_squared, 1), (x, y))
print("Recovered intermediate:", pair.inner)
print("Scalar candidates:", scalar_candidates)
print("Admissible levels:", admissible_levels)
print("Half-response contour:", contour)
