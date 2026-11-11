# Input limitations and conventions

Use this guide to recognize supported inputs and explicit refusals. The [capability matrix](capabilities.md) gives the broader mathematical scope; the [coefficient-field matrix](coefficient_fields_support.md) describes arithmetic restrictions.

## Computational boundaries

| Area | Supported behavior | Boundary |
| --- | --- | --- |
| Local intersections | Components whose support misses the selected point are removed before the local Tor calculation. The intersection of `y*(x-1)=0` and `(y-x)*(x-1)=0` at `(0,0)` has multiplicity one. | The intersection must be proper and isolated locally, with one complete-intersection presentation. |
| Absolute primary decomposition | Finite Artinian quotients and documented positive-dimensional families have geometric-primary certificates. Supported algebraic coordinate constants can be split over certified extensions. | General coupled nonprincipal absolute decomposition is unavailable. Unsupported geometric primes raise instead of being labeled geometrically prime. |
| Constructible projection | Polynomial and rational maps over `QQ` have exact finite locally closed image pieces. | A `ProjectionLimitError` reports exhaustion of the branch budget; it does not establish an empty image. |
| Coupled nonlinear families | Constructible dimension strata, bounded geometric cardinality strata and polynomial bases on supported parameter branches are available. | General component and singularity strata are not returned. See [fiber stratification](comprehensive_fibers.md) and [Gröbner covers](groebner_covers.md). |

## Input handling

Exact `CRootOf` coefficients are algebraic scalars even when the bound defining-polynomial symbol has the same name as a polynomial coordinate. Inexact coefficients are rejected.

A rational-map composition whose entire source lies in a pole of the outer map raises an empty-domain error before evaluating the outer expression. Preserve exclusions before SymPy simplifies them away; an algebraic field identity alone does not preserve the domain of a real function.

For a finite affine morphism, `ramification_locus` uses the zeroth Fitting ideal of relative Kähler differentials. This describes the non-unramified locus, including supported singular sources and nonreduced targets. Smoothness and flatness require their own hypotheses.
