"""Find special parameter values and keep degree drops in the result."""

import sympy as sp

from alggeom.decomposition import decomposition_condition
from alggeom.decomposition.parameters import polynomial_family

x, a, b = sp.symbols("x a b", real=True)
condition = decomposition_condition(x**4 + a * x, x, degrees=(2, 2), parameters=(a,))
assert condition.subs(a, 0) is sp.true
assert condition.subs(a, 1) is sp.false
family = polynomial_family(a * x**4 + b * x * x, x, parameters=(a, b))
assert family.verify()
assert family.select({a: 1, b: 2}).status == "decomposable"
assert family.select({a: 0, b: 2}).degree == 2
assert family.select({a: 0, b: 0}).degree == -1
print("Quartic decomposition condition:", condition)
