"""Calculate exact maps, kernels, tensor products and Hom of torsion modules."""

import sympy as sp

from alggeom.ideals import PolynomialRing
from alggeom.modules import PresentedModule, tensor_product

x, y = sp.symbols("x y")
S = PolynomialRing((x, y))
M = PresentedModule(S, [[x**2]])
N = PresentedModule(S, [[x]])
projection = M.hom([[1]], N)
assert projection.kernel.annihilator() == S.ideal(x)
assert projection.cokernel.is_zero
H = N.hom_module(M)
f = N.hom([[x]], M)
assert H.homomorphism(H.from_homomorphism(f))(N.gens[0]) == M([x])
assert tensor_product(M, N).annihilator() == S.ideal(x)
print(
    "S/(x^2) maps onto S/(x); its kernel, Hom in the reverse direction and tensor product have annihilator (x)."
)
