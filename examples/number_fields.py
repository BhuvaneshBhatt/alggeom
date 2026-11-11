"""Compare rational decompositions over QQ and a quadratic extension."""

import sympy as sp

from alggeom.decomposition.rational import search_rational

x, t = sp.symbols("x t")
function = x**2 + 2 / x**2
constants = sp.QQ.algebraic_field(sp.sqrt(2))
rational = search_rational(function, x, parameter=t)
extended = search_rational(function, x, parameter=t, constants=constants)
assert rational.complete and rational.verify() and len(rational.decompositions) == 1
assert extended.complete and extended.verify() and len(extended.decompositions) == 3

# The two additional substitutions require sqrt(2) in the coefficient field.
for sign in (-1, 1):
    inner = x + sign * sp.sqrt(2) / x
    outer = t**2 - sign * 2 * sp.sqrt(2)
    assert sp.cancel(outer.subs(t, inner) - function) == 0

print("Over QQ:", [(d.outer, d.inner) for d in rational.decompositions])
print("Over QQ(sqrt(2)):", [(d.outer, d.inner) for d in extended.decompositions])
