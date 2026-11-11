"""Decide when a quartic can be evaluated as two quadratic stages."""

import sympy as sp
from semialg.decision import equivalent

from alggeom.decomposition.parameters import polynomial_pair

x, t, a, b, c = sp.symbols("x t a b c", real=True)
polynomial = x**4 + a * x**3 + b * x**2 + c * x
pair = polynomial_pair(polynomial, x, degrees=(2, 2), parameters=(a, b, c), component_variable=t)
expected_condition = sp.Eq(a**3 - 4 * a * b + 8 * c, 0)
assert equivalent(pair.condition, expected_condition, (a, b, c))

assignment = {a: 2, b: 3, c: 2}
outer = pair.outer.subs(assignment)
inner = pair.inner.subs(assignment)
assert pair.condition.subs(assignment) is sp.true
assert sp.expand(outer.subs(t, inner) - polynomial.subs(assignment)) == 0
assert pair.condition.subs({a: 2, b: 3, c: 1}) is sp.false
print("Coefficient condition:", pair.condition)
print("Two stages at the selected coefficients:", outer, inner)
