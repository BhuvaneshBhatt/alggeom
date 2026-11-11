# Finite fields and point searches

Finite fields are useful for exact point counts, modular examples and checking geometric constructions in different characteristics. The field size and characteristic are different: GF(9) has characteristic three. Construct it as `finite_field(3, 2)`, not as the composite-modulus ring `sympy.GF(9)`.

```python
import sympy as sp
from alggeom.fields import finite_field, field_order
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme, ProjectiveScheme

x, y, z = sp.symbols("x y z")
k = finite_field(3, 2)
assert field_order(k) == 9
assert k.generator**9 == k.generator
C = AffineScheme(PolynomialRing((x, y), k).ideal(y - x**2))
assert len(C.rational_points()) == 9
P = ProjectiveScheme(PolynomialRing((x, y, z), k))
assert len(P.rational_points()) == 91
```

An explicit irreducible modulus can be supplied, for example `finite_field(3, modulus=sp.Symbol("a")**2+1)`. Reducible moduli and composite characteristics raise. The extension generator is part of the coefficient field and must be distinct from geometric coordinates. Field embeddings from GF(p) to a declared extension are supported; embeddings between different extension presentations must be supplied by the caller and are not guessed.

Polynomial ideals, elimination, intersections, quotients and saturation, coordinate-ring/module arithmetic, Hilbert invariants and standard-graded polynomial resolutions accept these fields. Radicals use repeated Frobenius preimages in any dimension. Finite quotient algebras use Frobenius kernels and fixed subalgebras for nilradicals and primary factors, including products without a separating linear element. Trace and norm remain algebraic trace and norm; the trace pairing can vanish on a nonzero algebra in positive characteristic.

Positive-dimensional primary decomposition and prime certification in finite characteristic are not yet implemented. Operations requiring those certificates, such as a general positive-dimensional fraction field, can raise even when simpler polynomial arithmetic succeeds. Positive-characteristic rational-function coefficient fields support polynomial/finite-quotient arithmetic, but general radical and primary algorithms over imperfect fields are outside the interface. Functional decomposition, plane-curve genus/rationality and general affine normalization retain their characteristic-zero scope.

## Completeness and budgets

Over GF(q), `rational_points` or `point_search` enumerates all field-valued points in any dimension, checking every affine tuple or canonical projective tuple. Affine enumeration has q^n candidates; projective enumeration has (q^n-1)/(q-1). `max_candidates` bounds this ambient search before enumeration. Exhausting the budget raises; it does not return an incomplete tuple as a complete point set. This is enumeration, not a specialized point-counting algorithm.

```python
Q = ProjectiveScheme(PolynomialRing((x, y)))
assert {p.coordinates for p in Q.point_search(1)} == {(1, -1), (1, 0), (1, 1), (0, 1)}
parabola = AffineScheme(PolynomialRing((x, y)).ideal(y - x**2))
assert {p.coordinates for p in parabola.point_search(2)} == {(-1, 1), (0, 0), (1, 1)}
```

Over QQ (the rationals ℚ), affine height H bounds the numerator magnitude and positive denominator of each reduced coordinate by H. Projective height bounds the largest absolute entry of a primitive integral representative; the first nonzero entry is chosen positive. Projective points are returned in the usual normalized rational-coordinate representation. A bounded search may miss rational points and never certifies global absence. Number-field height searches are not provided. The existing characteristic-zero `rational_points(field=K)` remains a complete enumeration for zero-dimensional schemes over the requested supported field.

Characteristic-zero parameter-dependent algebraic fields are available as `alggeom.fields.AlgebraicFunctionField`. Their polynomial rings preserve the declared coefficients and support explicit native basis construction. See the [algroots guide](https://github.com/BhuvaneshBhatt/algroots/blob/main/docs/algebraic-function-fields.md).
