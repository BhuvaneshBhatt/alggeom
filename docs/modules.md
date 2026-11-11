# Module arithmetic

Use a presented module to retain equations among vector-valued polynomial quantities. A matrix represents relations by columns; its rows correspond to module generators. This convention also makes Jacobian differentials and Fitting ideals agree with the geometry interfaces.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.modules import PresentedModule, direct_sum, tensor_product

x, y = sp.symbols("x y")
S = PolynomialRing((x, y))
M = PresentedModule(S, [[x**2]])
N = PresentedModule(S, [[x]])
projection = M.hom([[1]], N)
assert projection.is_surjective and not projection.is_injective
assert projection.cokernel.is_zero
assert projection.kernel.annihilator() == S.ideal(x)
assert projection.kernel_map(projection.kernel.gens[0]) == M([x])
assert projection(M([x])) == N.zero
```

`module(vector)` constructs a parent-aware element and reduces it modulo the presentation. Elements support addition, subtraction and scalar multiplication. `source.hom(matrix, target)` checks that every source relation maps to zero. Its `kernel`, `image` and `cokernel` are modules with explicit `kernel_map`, `image_map` and `cokernel_map`. `lift` expresses an image element in source generators; `inverse` requires an isomorphism. These presentations need not be minimal or unique.

```python
H = N.hom_module(M)
f = N.hom([[x]], M)
h = H.from_homomorphism(f)
assert H.homomorphism(h)(N.gens[0]) == M([x])
assert H.annihilator() == S.ideal(x)
assert tensor_product(M, N).annihilator() == S.ideal(x)
assert direct_sum(N, PresentedModule(S, [[y]])).annihilator() == S.ideal(x * y)
```

`HomModule` is a presented module whose elements convert to and from module homomorphisms. Tensor products are over the common base ring. For a module over S/I, annihilators and Fitting ideals are returned as inverse-image ideals in S, so they include I. Module maps, tensor and Hom respect the quotient relations; `restrict_scalars()` separately views the same object as a polynomial-ring module.

The implementation uses module standard bases and syzygies, so large presentations can be expensive. It supports exact fields and finite presentations over polynomial or coordinate rings. General Ext and Tor functors and local-ring module standard bases are not provided. Standard-graded polynomial resolutions, shifts and Betti tables are described in [resolutions](resolutions.md).
