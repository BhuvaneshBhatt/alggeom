"""Remove the generic double covering in the cusp parametrization (t**4,t**6)."""

import sympy as sp
from algroots.function_fields import luroth_generator

from alggeom.decomposition.rational import search_rational

t, s, u, v = sp.symbols("t s u v")
original_coordinates = (t**4, t**6)
search = search_rational(original_coordinates[0], t, parameter=s)
(pair,) = search.decompositions

# Field coordinates transport the other curve coordinate through the same substitution.
field = luroth_generator([pair.inner], variables=(t,), parameter=s)
new_coordinates = tuple(field.express(value).expression for value in original_coordinates)
new_x, new_y = new_coordinates
inverse_candidates = sp.solve(sp.together(new_y / new_x - v / u).as_numer_denom()[0], s)
(inverse_parameter,) = inverse_candidates

assert search.complete and field.verify()
assert all(
    sp.cancel(field.to_ambient(new) - old) == 0
    for new, old in zip(new_coordinates, original_coordinates, strict=True)
)
assert sp.cancel(new_y**2 - new_x**3) == 0
assert sp.cancel(inverse_parameter.subs({u: new_x, v: new_y}) - s) == 0
print("Original parametrization:", original_coordinates)
print("Recovered substitution:", field.generator)
print("New parametrization:", new_coordinates)
print("Rational inverse away from u=0:", inverse_parameter)
