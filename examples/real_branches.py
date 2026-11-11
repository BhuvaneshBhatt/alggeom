"""Select the principal fourth root and distinguish an incorrect branch."""

import sympy as sp

from alggeom.decomposition.algebraic import AlgebraicDecomposition
from alggeom.decomposition.real import verify_decomposition

x, t, z = sp.symbols("x t z", real=True)
pair = AlgebraicDecomposition(
    sp.Poly(x**4 - z, x, z, domain=sp.QQ),
    sp.Poly(x * x - t, x, t, domain=sp.QQ),
    sp.Poly(t * t - z, t, z, domain=sp.QQ),
    1,
    sp.Integer(1),
)
principal = verify_decomposition(
    pair, domain=z >= 0, target_selector=x >= 0, outer_selector=x >= 0, inner_selector=t >= 0
)
assert principal.equivalent and principal.certified and principal.function_verified
wrong = verify_decomposition(
    pair,
    domain=z > 0,
    target_selector=x >= 0,
    outer_selector=x >= 0,
    inner_selector=t < 0,
    check_functions=False,
)
assert not wrong.equivalent
print("Principal branch verified:", principal.function_verified)
