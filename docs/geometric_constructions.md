# Exact images, geometric components, fibers, and normalized curves

An algebraic map can have a nonclosed image, a polynomial can be irreducible over its coefficient field but reducible over an extension, and the dimension or cardinality of a fiber can change with parameters. These operations address those questions over exact coefficient fields. Their hypotheses matter: the interfaces below do not claim general algorithms for every finite-type morphism or arbitrary algebraic variety.

## Local intersection multiplicity

`intersection_multiplicity(left, right, point)` computes proper isolated intersections locally, even if the global intersection has positive-dimensional components elsewhere. It selects primary components supported at the chosen point before evaluating the Koszul Tor lengths. One input must still have a complete-intersection presentation; general Tor from arbitrary free resolutions is not available.

For `left = (y*(x-1))` and `right = ((y-x)*(x-1))`, their common line `x=1` prevents a global zero-dimensional calculation, but their isolated transverse intersection at `(0,0)` has multiplicty `1`.

## Absolute components and primary factors

`is_geometrically_integral` and `geometric_components` distinguish base-field irreducibility from decomposition over an algebraic closure. The geometric tests handle principal linear and quadratic hypersurfaces in at most two variables and supported finite ideals.

`absolute_primary_decomposition(ideal, max_root_degree=4)` now preserves nilpotent primary multiplicities. It handles bounded univariate polynomials over `QQ`, low-degree principal hypersurfaces, and zero-dimensional polynomial ideals over `QQ` whose finite algebra splits through a bounded exact algebraic extension. The zero-dimensional algorithm searches for a separating linear form, computes its multiplication-matrix characteristic polynomial, and constructs primary factors over a field containing all represented eigenvalues. Every proposed factor is checked for a geometrically rational reduced residue field, and their ideal intersection must reconstruct the original ideal. This includes nonmonogenic local Artinian algebras such as `QQ[x,y]/(x**2,y**2)`.

```python
import sympy as sp
from alggeom.geometric import absolute_primary_decomposition
from alggeom.ideals import PolynomialRing

x, y = sp.symbols("x y")
ring = PolynomialRing((x, y))
components = absolute_primary_decomposition(ring.ideal(x**2 - 1, y**2 - 1))
assert len(components.components) == 4
```

The certified absolute-primary cases also include positive-dimensional ideals with algebraic coordinate constants whose extended associated primes become affine-linear, including certain embedded components; positive-dimensional monomial ideals, primary factors supported on affine-linear geometric primes, bounded homogeneous binary forms, and products of powers of hypersurfaces whose irreducibility follows from a projective smoothness check. The smoothness test checks every standard projective chart, and component intersections are verified. The root-degree and separating-form search bounds are explicit. General absolute primary decomposition for arbitrary positive-dimensional ideals, high-degree multivariate absolute factorization, and unavailable exact algebraic roots remain unsupported; these cases raise instead of treating field-relative factors as absolute components.

## Locally closed images and Chevalley stratification

`LocallyClosedSet(I, nonzero)` represents `V(I)` with the listed nonvanishing conditions. `ConstructibleSet` supports finite unions, intersections, differences, geometric point membership, and Zariski closure.

The `constructible_image` of an `AffineMorphism` or `RationalMap` always retains an exact existential fiber predicate and its image closure. Calling `image.as_constructible()` or `image.strata()` constructs explicit locally closed image pieces over `QQ`. It first attempts faster specialized methods:

- Polynomial maps to the affine line with arbitrary polynomial source equations. Their image is finite or cofinite. Block Gröbner leading coefficients give a finite exceptional set; each irreducible exceptional factor is checked by exact algebraic fiber membership, subject to an algebraic-root degree limit.
- Maps to affine spaces whose graph equations reduce to a parameter-dependent affine-linear system after eliminating source variables with unit-linear relations. Determinantal rank conditions describe the image as a finite union of locally closed sets.

For example, `(s,t) -> (s,s*t)` has image `D(x) ∪ V(x,y)` although its image closure is all of `A^2`:

```python
from alggeom.schemes import AffineScheme

s, t, x, y = sp.symbols("s t x y")
source = AffineScheme(PolynomialRing((s, t)))
target = AffineScheme(PolynomialRing((x, y)))
image = source.morphism(target, (s, s * t)).constructible_image
strata = image.as_constructible()
assert strata.contains((0, 0))
assert not strata.contains((0, 1))
```

A general constructive elimination method handles arbitrary polynomial graph equations and nonvanishing conditions over `QQ`, including nonlinear multivariate polynomial maps and rational maps. It computes a finite locally closed image representation using parameter-sensitive univariate pseudoremainder branches. `max_branches` bounds the computation, and exceeding the budget raises an explicit error. The result is not necessarily an irredundant partition. See [constructive images](constructive_chevalley.md) for the algorithm and its coefficient-field scope.

## Parameter-dependent fibers in several dimensions

For one equation with one fiber variable and one parameter, `fiber_stratification(f, x, a)` retains the discriminant-based exact generic open set and exceptional fibers.

For tuples of equations, fiber variables, and parameters, `fiber_stratification(equations, variables, parameters)` computes a **complete rank stratification** when every equation is affine-linear in the fiber variables and polynomial in the parameters. It supports any finite number of fiber and parameter coordinates (subject to a minor-enumeration budget). Rank conditions partition parameter space into nonempty affine fibers with dimension `n-rank` and empty fibers where the augmented matrix has larger rank. A zero-dimensional consistent affine-linear fiber contains exactly one geometric point; a positive-dimensional fiber contains infinitely many.

```python
from alggeom.families import fiber_stratification

a, b, x, y = sp.symbols("a b x y")
family = fiber_stratification((a * x + b * y - 1, a * y), (x, y), (a, b))
assert family.specialize((1, 0)).geometric_points == 1
assert family.specialize((0, 1)).dimension == 1
assert family.specialize((0, 0)).dimension == -1
```

A single nonlinear equation supports any number of polynomial parameters and fiber variables. In one fiber variable, Sylvester rank conditions give the exact number of distinct geometric roots even at repeated-root and degree-drop strata. In several variables, a nonzero nonconstant hypersurface has dimension one less than the ambient fiber space. Systems of independent nonlinear univariate equations also admit exact product strata. See [nonlinear fibers](nonlinear_families.md). Coupled nonlinear systems can be stratified by dimension, with optional exact distinct-point counts within explicit Bézout and projection budgets; the elimination-based method does not return comprehensive Gröbner bases. See [coupled nonlinear fibers](comprehensive_fibers.md).

## Degree, finiteness, dominance, and ramification

`is_dominant` compares reduced scheme-theoretic image closures for general affine morphisms. For exact affine morphisms, `is_finite` tests module finiteness over the target coordinate ring using pure source-coordinate powers in the lexicographic initial ideal of the graph. It supports nonreduced targets and positive characteristic. For source morphisms into affine space over a characteristic-zero exact field:

- `morphism_degree` computes the length of the exact generic fiber over the rational-function field of the affine-space target. It applies to dominant generically finite maps, including birational morphisms that are not finite. For nonreduced generic fibers, the degree counts multiplicity.
- `is_birational` requires integral source and target and tests whether the dominant function-field extension has degree one.

For finite affine morphisms, `ramification_locus` is the non-unramified subscheme defined by `Fitt_0(Ω_(B/A))`, where `A` and `B` are the target and source coordinate rings. In smooth affine spaces of equal dimension this recovers the Jacobian determinant criterion. `branch_locus` computes its closed image. On a singular source, this detects failure of unramifiedness; it does not characterize flatness or smoothness.

For `(s,t)->(s**2,t)`, the generic degree is `2`, ramification is `s=0`, and the branch locus is `x=0`. For `(s,t)->(s,s*t)`, the generic degree is `1` and the map is birational but not finite.

## Explicit projective normalization gluing

`projective_normalization_atlas` normalizes each nonempty standard affine chart of an integral projective curve. `atlas.transition(i,j)` now returns **explicit regular transition coordinates** on the principal overlap `D(atlas.overlap_guard(i,j))`. Coordinates of the normalized target chart are transported through rational expressions in the original curve's function field, then represented in the source chart's localization by exact ideal-membership certificates. The transition is checked against the target normalization equations and the original projective-coordinate change. `atlas.verify_cocycle()` tests identities and compositions on all ordered triple overlaps.

These charts and regular gluing maps give an abstract affine-atlas presentation of the normalized projective curve, including positive-genus cases for which affine normalization completes. They do not construct a single embedding into projective space or an explicit homogeneous coordinate model for the normalization. Nonintegral projective curves and cases exceeding integral-closure limits are rejected.

## Mathematical references

See [references](references.md) for Cox, Little & O'Shea on elimination, Greuel & Pfister on primary decomposition and local algebra, and Eisenbud on integral closure, multiplicity, and finite maps. The algorithms here are certified for the input classes described above; resource exhaustion and unsupported hypotheses are not mathematical counterexamples.

## Geometric constants versus generic fiber splitting

For an integral curve such as `x**2-y=0`, passing from `QQ(y)` to an algebraic closure of the **function field** splits the generic zero-dimensional fiber in the variable `x`. It does not split the curve into absolute components over `QQbar`, because `sqrt(y)` is not a constant algebraic extension of `QQ`. The absolute-primary interface therefore never treats splitting over an algebraic closure of a generic parameter field as a certificate of geometric reducibility. Its positive-dimensional constant-field extension algorithm requires algebraic coordinate relations and checks every extended associated prime for a geometrically prime affine-linear presentation.
