# Coupled nonlinear fiber families

## When to use this

A polynomial system can change dimension, become empty, or acquire multiple solutions as parameters vary. A generic Gröbner basis over the rational-function field of the parameters does not describe what happens when its leading coefficients vanish. The `comprehensive_fiber_stratification` interface constructs **exact constructible parameter conditions** for dimension and, when requested within a declared complexity bound, the number of distinct geometric points. It accepts coupled nonlinear polynomial systems, not only independent equations or hypersurfaces.

All statements concern geometric fibers over an algebraic closure of `QQ`. Input coefficients must be rational, though the parameters may specialize to algebraic numbers. A fiber is allowed to be empty, nonreduced, or positive-dimensional. Inequations describe Zariski-open exclusions in each fiber.

## Example: an exceptional line

```python
import sympy as sp
from alggeom.comprehensive_strata import comprehensive_fiber_stratification

x, y, a = sp.symbols("x y a")
family = comprehensive_fiber_stratification((x * y - a, y), (x, y), (a,))

assert family.specialize((0,)).dimension == 1
assert family.specialize((1,)).dimension == -1
```

The family has the entire line `y=0` when `a=0` and no solutions otherwise. `pieces` contains constructible parameter loci paired with `FiberInvariants(dimension, geometric_points)`; the dimension of an empty fiber is `-1`.

For one-dimensional zero-dimensional fibers, the geometric-point count can also be computed:

```python
x, a = sp.symbols("x a")
family = comprehensive_fiber_stratification(
    (x**2 - a,), (x,), (a,), count_points=True, max_points=2
)
assert family.specialize((0,)).geometric_points == 1
assert family.specialize((1,)).geometric_points == 2
```

The count is the number of **distinct geometric points**, not the scheme-theoretic length. For the fiber `x**2=0`, the point count is one even though the algebra has length two. Positive-dimensional fibers have infinitely many geometric points (`sp.oo`). With `count_points=False` (the default), zero-dimensional point counts are `None`, meaning *not computed*, not zero.

## Exact algorithm and mathematical guarantee

Let `X_a` be the constructible fiber in affine `n`-space. For each `d` from zero through `n`, consider all coordinate projections onto `d` coordinates. The elementary algebraic-matroid fact that a coordinate transcendence basis exists gives

`dim(X_a) >= d` exactly when **at least one** such coordinate projection is Zariski dense in affine `d`-space.

The package applies the constructive Chevalley algorithm to each projected family, obtaining finite unions of locally closed sets `V(I) ∩ D(g_1) ∩ ... ∩ D(g_r)` in parameter and projected coordinates. On a fixed parameter fiber, one piece is dense precisely when all its defining equations specialize to zero polynomials in the projected coordinates and each inequation specializes to a nonzero polynomial. Testing polynomial coefficients produces an exact constructible condition on the parameters. Subtracting successive conditions gives the dimension strata, including the empty-fiber locus. No numerical sampling or generic-fiber assumption enters this construction.

When point counts are requested, the condition `#X_a >= k` is expressed existentially with `k` solution tuples and a linear functional separating those tuples. A finite deterministic family of linear forms suffices over characteristic zero, because a nonzero polynomial in the form parameter has finitely many roots. Constructive Chevalley elimination then provides the exact `>=k` locus. The resulting loci partition the zero-dimensional part by the number of distinct solutions. The affine Bézout bound `D**n`, with `D` at least the maximum fiber-equation degree, supplies a certified upper bound on the size of a finite geometric fiber. The function refuses a complete point-count request if that bound exceeds `max_points`.

`max_branches` bounds each Chevalley projection. Exhaustion raises `ProjectionLimitError`; an excessive point-count degree bound raises `StratificationLimitError`. The routine never labels a partial search as a complete stratification. Cardinality computation may be dramatically more expensive than the dimension filtration.

## Relation to comprehensive Gröbner systems

**This is a fiber-invariant algorithm, not a Gröbner-cover algorithm.** A Gröbner cover supplies parameter regions together with polynomial data that specialize to Gröbner bases for every point of those regions. The current algorithm instead computes exact constructible loci for geometric dimensions and optional cardinalities using first-order elimination. It does **not** return parametric Gröbner bases, Hilbert functions, primary components, or arbitrary singularity strata. The separate [Gröbner-cover interface](groebner_covers.md) supplies certified bases for parameter branches it can represent exactly; nonlinear irreducible exceptional parameter varieties remain unsupported.

For small linear systems, use `linear_fiber_stratification`; for a single equation, use `hypersurface_fiber_stratification`. Both may be much faster. `families.fiber_stratification` chooses the existing specialized methods first and falls back to the coupled elimination algorithm for more general nonlinear systems.

## References

- Antonio Montes and Michael Wibmer (2010), *Gröbner Bases for Polynomial Systems with Parameters*, Journal of Symbolic Computation 45(12), 1391–1425. https://doi.org/10.1016/j.jsc.2010.06.017
- Deepak Kapur, Yao Sun, and Dingkang Wang (2013), *An Efficient Algorithm for Computing a Comprehensive Gröbner System of a Parametric Polynomial System*, Journal of Symbolic Computation 49, 27–44. https://doi.org/10.1016/j.jsc.2011.12.015

These references describe true comprehensive Gröbner systems and explain the additional polynomial-basis data needed for that stronger interface.
