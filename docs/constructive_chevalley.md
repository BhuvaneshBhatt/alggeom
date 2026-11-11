# Constructive images of affine algebraic maps

## Why compute constructible images?

Eliminating the source variables from the graph of a polynomial map computes the **Zariski closure** of its image. This can include points with no preimage. For example, the map `(s, t) -> (s, s*t)` has dense image in the plane but misses every `(0, y)` with `y != 0`. A constructible image records this difference exactly.

## Exact polynomial projection

`alggeom.chevalley.project_constructible(source, variables, max_branches=4096)` projects a finite union of locally closed subsets of affine space by eliminating the specified coordinates. Every input piece has the form `V(I) ∩ D(g_1) ∩ ... ∩ D(g_r)`. The output is an explicitly represented finite union of locally closed sets, not only an existential predicate. The supported ground field is `QQ`, and the result describes points over its **algebraic closure**. It does not test membership over `QQ` or `RR` alone.

The procedure eliminates coordinates one at a time. For one variable `x`, the fiber equations specialize to univariate polynomials. Branching on leading-coefficient vanishing ensures that each branch has a fixed degree; pseudoremainder sequences then compute the specialized gcd without dividing by a coefficient that might vanish on the branch. If the equations vanish identically, a fiber point exists precisely when the product of the nonvanishing constraints is not identically zero. Otherwise the fiber exists when the gcd has a root outside the excluded divisor. To account for repeated roots, the algorithm compares the gcd with the **degree-th power of the exclusion polynomial**. These are exact calculations in characteristic zero.

Every split is recorded as polynomial equations and inequations in the remaining coordinates. Branches inconsistent over an algebraic closure are removed by ideal saturation. Repeated elimination provides a constructive Chevalley procedure for arbitrary affine polynomial equation and inequation presentations over `QQ`. `ProjectionLimitError` means the branch budget was exceeded; it is not evidence that the image is empty.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme

s, t, x, y = sp.symbols("s t x y")
source = AffineScheme(PolynomialRing((s, t)))
target = AffineScheme(PolynomialRing((x, y)))
image = source.morphism(target, (s**2, s * t)).constructible_image
strata = image.as_constructible(max_branches=4096)
assert strata.contains((0, 0))
assert strata.contains((4, 3))
assert not strata.contains((0, 3))
```

Rational maps with explicitly excluded denominator loci also use exact projection. For instance, the map `s -> 1/s`, defined for `s != 0`, has image `D(u)`.

## Computational limits

The procedure terminates mathematically because each univariate Euclidean step lowers degree and each zero-leading-coefficient branch lowers the specialized degree. It can nevertheless produce exponentially many locally closed pieces, with expensive Gröbner-basis saturation and polynomial coefficient growth. The `max_branches` limit is explicit, and exceeding it raises instead of returning an incomplete answer. Current exact computation is restricted to `QQ` coefficients; extensions and positive characteristic need separate arithmetic support. The output pieces are not necessarily disjoint or irredundant, and no minimal decomposition is promised.

For related interfaces, see [geometric constructions](geometric_constructions.md), [API reference](api.md), and [limitations](known_issues.md).
