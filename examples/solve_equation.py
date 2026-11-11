"""Solve an expanded sextic by recovering a cubic after a quadratic."""

import sympy as sp

from alggeom.decomposition import decompose

x = sp.Symbol("x")
polynomial = sp.expand((x**2 + x) ** 3 - (x**2 + x))
components = decompose(polynomial, x)

# Invert the outer map first, then solve each inner preimage equation.
preimages = [sp.S.Zero]
for component in components:
    argument = component.variables[0]
    preimages = [root for value in preimages for root in sp.solve(component.expr - value, argument)]
solutions = frozenset(preimages)
real_solutions = frozenset(root for root in solutions if root.is_real is True)

assert len(solutions) == 6
assert all(sp.simplify(polynomial.subs(x, root)) == 0 for root in solutions)
assert sp.Poly(sp.prod(x - root for root in solutions).expand(), x) == sp.Poly(polynomial, x)
print("Expanded equation:", sp.Eq(polynomial, 0))
print("Component maps:", components)
print("Real roots:", sorted(real_solutions, key=sp.default_sort_key))
