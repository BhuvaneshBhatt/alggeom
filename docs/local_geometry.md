# Rational functions and local geometry

Local calculations describe what happens near a chosen point. A cusp has dimension one, but needs two tangent directions at its singular point. Its tangent cone retains a double tangent line, which determines multiplicity.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme
from alggeom.local import intersection_multiplicity

x, y, z = sp.symbols("x y z")
S = PolynomialRing((x, y))
C = AffineScheme(S.ideal(y**2 - x**3))
assert C.local_dimension((0, 0)) == 1
assert C.embedding_dimension((0, 0)) == 2
assert C.tangent_cone((0, 0)) == S.ideal(y**2)
assert C.multiplicity((0, 0)) == 2
assert not C.is_smooth((0, 0)) and C.is_smooth((1, 1))
assert intersection_multiplicity(C, S.ideal(y), (0, 0)) == 3
```

`tangent_cone` translates the point, introduces a scale parameter, saturates and specializes to the lowest-degree initial ideal. Taking only the initial forms of the supplied generators can miss cancellations; the saturation handles the entire ideal. Hilbert degree of that graded cone is the Hilbert–Samuel multiplicity. Projective inputs use an affine chart containing the point. Points may carry an explicit coefficient-field extension.

`is_reduced`, `is_irreducible`, `is_integral`, `is_equidimensional`, `is_unmixed`, `codimension` and `is_smooth` are available as ring/scheme members and through `alggeom.properties`. Irreducibility is relative to the declared field. Equidimensionality compares minimal components; unmixedness also checks associated primes, including embedded ones. The empty scheme is reduced, normal, smooth, equidimensional and unmixed; it is not integral or irreducible, and its embedded codimension is infinity.

## Function fields and regular fractions

```python
A = C.coordinate_ring
K = A.fraction_field()
t = K(y / x)
assert t**2 == K(x) and t**3 == K(y)
assert t.inverse() * t == 1

B = PolynomialRing((x, y, z)).quotient(x * y - z)
local = AffineScheme(B).local_ring((0, 0, 0))
assert local(z / x) == local(y)
assert local(1 / (1 + x)).is_unit
assert not local(x).is_unit
```

A `FractionField` requires a domain. `A.localize(f)` inverts powers of f, records the kernel I:f^∞, and supplies an explicit affine model A[u]/(uf-1). `LocalRing(A, prime)` inverts the complement of a certified proper prime; `local_ring(A, point)` constructs the local ring at an exact point. Its residue field is available separately. Fractions compare by cross multiplication in the appropriate localized quotient.

A vanishing displayed denominator can still represent a regular function, as z/x=y on xy=z. Regularity is solved using ideal quotients and membership witnesses. A genuine pole raises. Denominators that are zero divisors in the localized ring are rejected before cancellation: x/x is not an admissible total fraction at the crossing of xy=0. When that distinction matters, use `.fraction(numerator, denominator)` to preserve a denominator that SymPy would otherwise cancel before the call.

## Length and intersection multiplicity

`module_length(M)` returns the coefficient-space dimension of a finite-support module. At a specified field-valued point, it selects the primary local factor, giving the local module length; a point outside the support has length zero. Global coefficient-space dimension includes residue-field degrees. It rejects positive-dimensional support instead of treating a finite spanning approximation as a length.

`intersection_multiplicity(I, J, point)` supports proper finite intersections in a common affine space when one input is a complete intersection in its supplied equations. It computes the alternating lengths of the Koszul Tor homology, so it does not assume the answer is always the length of S/(I+J). Nonproper intersections, nonfinite intersection support, and inputs needing a general free resolution raise. This follows the [local intersection multiplicity convention](https://doc.sagemath.org/html/en/reference/schemes/sage/schemes/affine/affine_subscheme.html).

No completions, general local standard bases or arbitrary local-ring resolutions are implemented. Prime certification and primary-factor selection retain the [coefficient-field limits](capabilities.md) of the component interface. Principal localization witness searches and finite module linear algebra have explicit resource limits. Generic calculations over K(t) do not specialize exceptional parameter values.
