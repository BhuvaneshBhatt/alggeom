# Nonlinear parameter-dependent fibers

## Why stratify parameter families?

The dimension and number of distinct geometric points of a polynomial fiber may change when leading coefficients vanish or roots collide. A generic Gröbner basis alone does not capture all exceptional parameter values. These algorithms describe parameter conditions on which the stated invariants are exact.

## Polynomial hypersurfaces with any number of parameters

`alggeom.nonlinear_strata.hypersurface_fiber_stratification(equation, variables, parameters)` handles one polynomial equation over `QQ`, with arbitrary polynomial dependence on any finite set of fiber and parameter variables.

For a **single fiber variable**, the algorithm partitions parameter space by the degree of the specialized polynomial and the rank of its Sylvester matrix with its derivative. In characteristic zero, if a nonzero specialized polynomial has degree `d`, its number of distinct geometric roots equals `rank(Sylvester(f, f')) - d + 1`. Rank conditions are expressed as determinantal equations and nonvanishing minors. The cases of a nonzero constant (empty fiber) and the zero polynomial (entire affine line) are included.

For **two or more fiber variables**, a single nonzero nonconstant polynomial cuts out a hypersurface of dimension `n-1` over the algebraic closure. The zero polynomial gives the whole `n`-dimensional affine space; a nonzero constant gives the empty set. Coefficient-vanishing conditions give complete parameter strata for these three possibilities. This computes dimension and whether the geometric point set is infinite, not finer invariants such as connectedness or singularity.

```python
import sympy as sp
from alggeom.families import fiber_stratification

x, a, b = sp.symbols("x a b")
family = fiber_stratification((a * x * x + b * x + 1,), (x,), (a, b))
assert family.specialize((1, 0)).geometric_points == 2
assert family.specialize((1, 2)).geometric_points == 1
assert family.specialize((0, 2)).geometric_points == 1
assert family.specialize((0, 0)).dimension == -1
```

## Products of independent nonlinear systems

`alggeom.nonlinear_strata.separable_fiber_stratification(equations, variables, parameters)` combines these exact univariate partitions for systems where each equation involves exactly one fiber coordinate and no coordinate appears in two equations. The resulting product strata preserve dimension and multiply finite geometric point counts. For example, `x**2-a = 0` and `y**2-b = 0` have four points when both parameters are nonzero, two when exactly one is zero, and one when both are zero.

The `fiber_stratification` dispatcher uses the specialized hypersurface, affine-linear, or separable-nonlinear algorithm when applicable, then falls back to exact constructible elimination for coupled nonlinear systems such as `(x**2+y**2-a, x*y-b)`. The fallback computes dimension strata; request exact finite cardinalities through the dedicated API with explicit budgets.

## Bounds and remaining cases

Enumerating minors grows combinatorially. The `max_minors` and `max_pieces` limits prevent unexpectedly large computations; a limit error is distinct from an empty fiber. Coupled nonlinear systems can be stratified by constructive elimination using [comprehensive fiber stratification](comprehensive_fibers.md). This gives exact dimension conditions and optionally exact distinct-point counts, subject to Bézout and projection budgets. A true comprehensive Gröbner cover, with specialized Gröbner bases on every parameter stratum, is not yet available.

See [constructive images](constructive_chevalley.md), [geometric constructions](geometric_constructions.md), and the [capability summary](capabilities.md).

## Coupled systems

For several nonlinear equations sharing fiber coordinates, [`comprehensive_fiber_stratification`](comprehensive_fibers.md) uses exact constructible projection to compute fiber dimensions. With explicit budgets it can also determine the number of distinct points of finite fibers. This is an elimination-based invariant stratification, not a comprehensive Gröbner basis or primary-decomposition stratification.
