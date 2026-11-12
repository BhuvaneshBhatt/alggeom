# API reference

This reference is generated from public signatures and docstrings by `python tools/generate_api.py`. Use `python tools/generate_api.py --check` to detect stale documentation.

Start with the [geometry tutorial](geometry.md) for algebraic conventions, then [schemes and fibers](schemes.md), [components](components.md), [resolutions](resolutions.md), [curves](curves.md), [module arithmetic](modules.md), [local geometry](local_geometry.md), [finite fields](finite_fields.md) and [normalization](normalization.md), or the [decomposition guide](guide.md) for functional composition. See [result semantics](results.md) for decomposition certificates and the [coefficient-field tutorial](coefficient_fields.md) for implicit reconstruction. The decomposition `parameters` and `real` modules require the `real` extra. Generating the complete reference also requires that extra.

## `alggeom.absolute_factorization`

Absolute factorization of bivariate polynomials using logarithmic derivatives.

For a squarefree polynomial f in QQ[x, y], the space of closed rational one-forms with poles along f and bounded bidegrees has one dimension per geometric irreducible factor. Its multiplication algebra separates factors using exact rational linear algebra and algebraic eigenvalues. Every factor is recovered by polynomial gcd and checked by exact reconstruction.

### `AbsoluteFactors`

```text
class AbsoluteFactors(factors, field)
```

Geometrically irreducible factors and their multiplicities over a number field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/absolute_factorization.py#L19)

Stored fields: `factors`, `field`.

### `FactorizationLimitError`

```text
class FactorizationLimitError(NotImplementedError)
```

Exact factor reconstruction exceeds a declared algebraic complexity limit.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/absolute_factorization.py#L15)

### `absolute_factor_bivariate`

```text
absolute_factor_bivariate(expression, variables, *, max_unknowns=100, max_roots=12)
```

Factor a bivariate QQ polynomial into certified geometric factors.

Multiplicities are preserved. Factors are represented over a common exact number field. The method uses logarithmic derivatives and characteristic- zero residue-algebra splitting; it does not rely on factorization over a generic rational-function field, which could split a geometrically prime polynomial incorrectly.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/absolute_factorization.py#L158)

### `absolute_factor_hypersurface`

```text
absolute_factor_hypersurface(expression, variables, *, max_unknowns=90, max_roots=12)
```

Factor a rational-coefficient multivariate hypersurface geometrically.

For two variables the optimized bivariate backend is used. In higher dimensions, solve the closed-form equations for logarithmic differentials in all coordinates, then reconstruct geometric factors from residue eigenvalues and exact multivariate gcds. Only characteristic zero and rational coefficients are supported; expensive systems raise explicitly.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/absolute_factorization.py#L203)

## `alggeom.chevalley`

Constructive existential elimination over algebraically closed characteristic-zero fields.

The implementation uses comprehensive univariate polynomial gcds. Parameter coefficient vanishing is split into locally closed loci, and pseudoremainders are computed only on branches where the divisor's leading coefficient is nonzero. Repeating the one-variable projection eliminates arbitrary finite sequences of existential polynomial variables. All computations are exact over QQ and interpreted geometrically over its algebraic closure.

### `ProjectionLimitError`

```text
class ProjectionLimitError(RuntimeError)
```

The constructive projection exceeded its explicit branch budget.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/chevalley.py#L19)

### `project_constructible`

```text
project_constructible(constructible, variables, *, max_branches=4096)
```

Project a constructible set along specified coordinates exactly.

Input pieces are finite intersections of polynomial zero and nonzero conditions over QQ. The result is a finite union of locally closed sets, equal to the existential image over an algebraic closure of QQ. Every parameter-specialization branch is resolved using exact polynomial gcd arithmetic; hitting `max_branches` raises instead of returning an incomplete image. Complexity can be very high for nonlinear systems.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/chevalley.py#L106)

## `alggeom.components`

Exact radicals and primary components with coefficient-field-specific algorithms.

### `DecompositionLimitError`

```text
class DecompositionLimitError(RuntimeError)
```

An ideal decomposition exceeds a declared resource limit; no partial answer is returned.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L16)

### `radical`

```text
radical(ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64)
```

Compute the radical by characteristic-zero components or finite-field Frobenius roots.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L252)

### `minimal_primes`

```text
minimal_primes(ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64)
```

Return the prime ideals of the reduced irreducible components over the declared field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L244)

### `associated_primes`

```text
associated_primes(ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64)
```

Return all associated primes, including embedded primes, over the declared field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L232)

### `primary_decomposition`

```text
primary_decomposition(ideal, *, max_dimension=256, max_nodes=128, max_exponent=256, max_separating_attempts=64)
```

Return an irredundant primary decomposition, including embedded components.

Generic localization at independent coordinates reduces to finite algebras. Saturation contracts the generic primary factors; a certified exponent splits off the remaining closed locus. Components are over the declared field, not automatically over its algebraic closure. Resource exhaustion raises and never returns an apparently complete prefix. Decompositions of embedded components need not be unique.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L212)

### `is_prime`

```text
is_prime(ideal, **limits)
```

Decide whether a proper ideal is prime over its declared coefficient field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L291)

### `is_primary`

```text
is_primary(ideal, **limits)
```

Decide whether a proper ideal is primary, retaining nilpotents.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/components.py#L301)

## `alggeom.comprehensive_strata`

Exact geometric fiber stratification via constructive algebraic elimination.

Every input is a finite polynomial equation/inequation system over QQ, with arbitrary coupling among the fiber coordinates. The algorithm computes the fiber-dimension filtration by detecting dense coordinate projections. It can also partition the zero-dimensional locus by its number of geometric points, using an explicitly bounded sequence of existential distinct-point tests.

These are complete constructible conditions whenever the branch and cardinality budgets allow completion. This is a quantifier-elimination algorithm, not a comprehensive Groebner basis: it does not return parametric Groebner bases.

### `StratificationLimitError`

```text
class StratificationLimitError(RuntimeError)
```

The requested complete stratification exceeds a declared degree/branch limit.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/comprehensive_strata.py#L25)

### `GeometricFiberPiece`

```text
class GeometricFiberPiece(locus, invariants)
```

A constructible parameter locus with fixed dimension and point count.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/comprehensive_strata.py#L162)

Stored fields: `locus`, `invariants`.

### `ComprehensiveFiberStratification`

```text
class ComprehensiveFiberStratification(equations, variables, parameters, nonzero, pieces)
```

Exact constructible dimension and cardinality partition of affine fibers.

All geometric quantities refer to the algebraic closure of QQ. The pieces may be represented by overlapping locally closed terms internally, but distinct invariant values occur on disjoint parameter loci.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/comprehensive_strata.py#L170)

Stored fields: `equations`, `variables`, `parameters`, `nonzero`, `pieces`.

#### `ComprehensiveFiberStratification.specialize`

```text
specialize(self, values)
```

Return the certified invariants of an exact algebraic parameter tuple.

### `comprehensive_fiber_stratification`

```text
comprehensive_fiber_stratification(equations, variables, parameters, *, nonzero=(), max_branches=4096, count_points=False, max_points=4)
```

Partition arbitrary coupled polynomial families by exact fiber invariants.

Constructs the dimension filtration using Chevalley elimination for every coordinate projection. Over an algebraically closed field, a fiber has dimension at least d iff one projection onto d of its coordinates is dense. For zero-dimensional fibers, at least-k distinct-point loci are computed by existential elimination and refined into exact cardinality strata.

`max_points` bounds the Bézout upper bound, not a sampling or guessing cutoff. When a requested complete answer exceeds this bound or Chevalley exceeds `max_branches`, the function raises without returning partial data. Set `count_points=False` to compute dimensions alone for expensive inputs.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/comprehensive_strata.py#L205)

## `alggeom.constructible`

Locally closed algebraic sets and exact existential images.

ConstructibleSet stores a finite union of V(I) intersect D(g_1)...D(g_r). ConstructibleImage retains an existential presentation when a finite locally closed decomposition has not been computed. Point membership is over an algebraic closure of the declared exact coefficient field, not merely the rational points of that field.

### `LocallyClosedSet`

```text
class LocallyClosedSet(ideal, nonzero=())
```

The algebraic set V(ideal) with each supplied polynomial nonzero.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/constructible.py#L17)

Stored fields: `ideal`, `nonzero`.

#### `LocallyClosedSet.closure_ideal`

The Zariski closure, including removal of locally empty components.

#### `LocallyClosedSet.contains`

```text
contains(self, values, *, field=None)
```

Test an exact field-valued point, using a specified extension if needed.

### `ConstructibleSet`

```text
class ConstructibleSet(pieces)
```

A finite union of locally closed subsets of one affine ambient space.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/constructible.py#L53)

Stored fields: `pieces`.

#### `ConstructibleSet.ring`

The common polynomial ring, if this set has a nonempty representation.

#### `ConstructibleSet.closure_ideal`

Compute the Zariski closure by saturation and ideal intersection.

#### `ConstructibleSet.contains`

```text
contains(self, values, *, field=None)
```

Check membership in at least one of the locally closed pieces.

### `ConstructibleImage`

```text
class ConstructibleImage(mapping)
```

An exact existential image of a polynomial or rational affine map.

Membership is decided by solving the exact fiber equations over an algebraic closure. The representation is not yet a finite decomposition into locally closed subsets; `closure_ideal` is merely its closure.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/constructible.py#L128)

Stored fields: `mapping`.

#### `ConstructibleImage.closure_ideal`

Return the Zariski closure ideal, not a substitute for the exact image.

#### `ConstructibleImage.as_constructible`

```text
as_constructible(self, *, max_minors=4096, max_exceptional_degree=4, max_branches=4096)
```

Return explicit exact locally closed image pieces over QQ.

Uses structural linearization and rank conditions when possible, followed by constructive Chevalley elimination for nonlinear polynomial or rational maps. Branch limits raise explicitly.

#### `ConstructibleImage.strata`

```text
strata(self, **options)
```

Return a finite locally closed image decomposition when computable.

#### `ConstructibleImage.contains`

```text
contains(self, values, *, field=None)
```

Decide geometric-fiber nonemptiness for an exact target point.

## `alggeom.cover_proofs`

Exact parameter-partition proofs and replayable polynomial reduction identities.

The verifier uses ideal arithmetic on constructible loci, exact polynomial identities, and Buchberger's criterion. It does not re-run the cover dispatcher.

### `CoverProofError`

```text
class CoverProofError(ValueError)
```

A cover condition or certificate has no supported exact proof.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/cover_proofs.py#L15)

### `ReductionIdentity`

```text
class ReductionIdentity(dividend, divisors, quotients, remainder)
```

A polynomial equality dividend = sum(quotient*divisor) + remainder.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/cover_proofs.py#L166)

Stored fields: `dividend`, `divisors`, `quotients`, `remainder`.

#### `ReductionIdentity.verify`

```text
verify(self)
```

Check the stored polynomial identity exactly over its coefficient field.

### `StratumProof`

```text
class StratumProof(condition, variables, parameters, basis, source, basis_membership, source_reductions, critical_pairs, nonzero_guards, parameter_bindings=())
```

An exact Gröbner-basis proof on a constructible parameter stratum.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/cover_proofs.py#L193)

Stored fields: `condition`, `variables`, `parameters`, `basis`, `source`, `basis_membership`, `source_reductions`, `critical_pairs`, `nonzero_guards`, `parameter_bindings`.

#### `StratumProof.verify`

```text
verify(self)
```

Check all identities, leading terms and guards without Buchberger search.

### `GroebnerCoverProof`

```text
class GroebnerCoverProof(pieces, parameters, equations, variables)
```

A parameter partition with exact polynomial proofs for every basis.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/cover_proofs.py#L458)

Stored fields: `pieces`, `parameters`, `equations`, `variables`.

#### `GroebnerCoverProof.verify`

```text
verify(self)
```

Replay source identities, reductions and the exact parameter partition.

### `check_parameter_partition`

```text
check_parameter_partition(pieces, parameters)
```

Prove coverage and pairwise disjointness over an algebraic closure.

The proof is set-theoretic and uses saturation, ideal intersections and exact constructible-set difference; it is not based on test points.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/cover_proofs.py#L140)

### `prove_stratum`

```text
prove_stratum(piece, equations, variables, parameters)
```

Construct identities using a traced Buchberger computation on one branch.

The trace expresses every basis element in the original generators. Subsequent source and S-polynomial reductions certify ideal equality and Buchberger's criterion. All coefficient denominators are recorded.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/cover_proofs.py#L336)

## `alggeom.curves`

Plane-curve genus, normalization analysis and genus-zero parametrization.

Models are reduced irreducible plane equations over QQ or number fields. Normalization analysis checks geometric integrality and resolves singularities to compute the smooth projective genus and adjoints. It does not construct a general affine integral closure. Rational curves additionally have explicit birational parametrizations and projective normalization morphisms.

### `PlaneCurve`

```text
class PlaneCurve(equation, variables, domain=QQ)
```

A reduced irreducible affine plane model over QQ or a number field.

The affine equation and its projective closure share a function field. Genus refers to the smooth projective normalization, so singular plane degree alone does not determine it. Geometrically reducible models are rejected by normalization analysis. Rational-function constant fields, arbitrary space curves and higher-dimensional normalization are outside this interface's scope.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/curves.py#L21)

Stored fields: `equation`, `variables`, `domain`.

#### `PlaneCurve.degree`

The total degree of the plane model, not the normalized curve's genus.

#### `PlaneCurve.affine_scheme`

The affine closed scheme defined by this equation in A^2.

#### `PlaneCurve.projective_scheme`

The projective closure, with a fresh final homogenizing coordinate.

#### `PlaneCurve.genus`

The genus of the smooth projective normalization over the declared constants.

#### `PlaneCurve.normalization`

```text
normalization(self, *, max_steps=256, max_dimension=256)
```

Return algroots' exact CurveNormalization analysis certificate.

It records genus, arithmetic_genus, delta_sum, singularities including infinitely near points, and canonical_adjoints, with verify() replay. This is resolution and adjoint data, not an integral-closure coordinate ring. Computation limits raise without a partial genus claim.

#### `PlaneCurve.rationality`

```text
rationality(self, parameter=None, *, point=None, extend_constants=False, max_steps=256, max_dimension=256)
```

Return an exact CurveRationalityResult with maps or a rationality obstruction.

Status is rational, rational_over_extension, positive_genus, or no_ground_field_point. A genus-zero obstruction over QQ is decided by a conic model. Other number fields need a supplied point witness or extend_constants=True; failure to find a point is not called an obstruction. Any constants extension and the inverse birational map remain available in this specialist certificate.

#### `PlaneCurve.parametrize`

```text
parametrize(self, parameter=None, **options)
```

Return a RationalMap from an open subset of A^1 onto this plane curve.

Denominator exclusions describe this affine parameter chart. Use rationality() for inverse coordinates and exceptional points, or normalization_map() to include the projective parameter at infinity. A positive genus or missing ground-field point raises ValueError.

#### `PlaneCurve.affine_normalization`

```text
affine_normalization(self, **limits)
```

Return the affine integral-closure scheme, including positive-genus plane models.

#### `PlaneCurve.projective_normalization_atlas`

```text
projective_normalization_atlas(self, **limits)
```

Return exact affine normalization maps covering the projective model.

The returned atlas records normalized affine patches and overlap opens. It does not provide explicit transition coordinates for newly adjoined integral generators or a single projective source embedding.

#### `PlaneCurve.normalization_map`

```text
normalization_map(self, *, parameters=None, **options)
```

Return the projective normalization morphism P^1 -> C for a rational curve.

The target is base changed if extend_constants=True supplies a parametrization over a larger field. Positive-genus curves currently supply normalization analysis only, without a normalized scheme model.

### `genus`

```text
genus(curve)
```

Return the smooth projective genus of a PlaneCurve.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/curves.py#L149)

### `normalization`

```text
normalization(curve, **limits)
```

Return exact plane-curve normalization analysis, including singularities and adjoints.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/curves.py#L154)

### `parametrize`

```text
parametrize(curve, parameter=None, **options)
```

Return a rational affine parametrization of a genus-zero plane curve.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/curves.py#L159)

## `alggeom.decomposition`

Exact compositions of rational maps and implicit algebraic functions.

### `common_right_component`

```text
common_right_component(functions, variables, *, constants=QQ, parameter=None)
```

Return a Lambda generating the joint field of the supplied rational functions.

Return None for constant-only inputs or a joint field of transcendence degree greater than one. A degree-one Lambda is a valid generator. Use common.common_components for outer maps and evidence.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_api.py#L11)

### `decompose`

```text
decompose(expression, variables, *, constants=QQ, **limits)
```

Return an outer-to-inner tuple of SymPy Lambdas for one resolved chain.

`variables` is a Symbol or an ordered sequence of input coordinates. Only the last map may have several arguments. Constant and indecomposible inputs give a one-element tuple. Search limits are those of `rational_chain`. An unfinished chain raises `IncompleteDecompositionError` with its certificate in `error.result`.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_api.py#L21)

### `decompositions`

```text
decompositions(expression, variables, *, parameter=None, constants=QQ, **limits)
```

Return every nontrivial two-map pair over the declared constant field.

Each pair is `(outer, inner)` with SymPy Lambda maps. Answers are distinct up to a fractional-linear change of the intermediate coordinate. An empty tuple proves that no such pair exists in this scope. Incomplete searches raise `IncompleteDecompositionError`; their partial answers remain in `error.result`. Limits are those of `rational.search_rational`.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_api.py#L38)

### `decompose_relation`

```text
decompose_relation(relation, variables, *, constants=QQ, **limits)
```

Return an ordered tuple of implicit component polynomials.

`variables=(output,input)` fixes the orientation. Adjacent polynomials share an intermediate coordinate. Their iterated resultant is a nonzero constant times a positive power of the original relation. Completeness concerns the selected chain under coefficient-field and rational-graph search. An unresolved leaf raises `IncompleteDecompositionError`; the exact partial chain remains available in `error.result`.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_api.py#L58)

### `decomposition_condition`

```text
decomposition_condition(polynomial, variable, *, degrees, parameters=())
```

Return the exact real parameter condition for fixed polynomial degrees.

`degrees=(outer_degree,inner_degree)` fixes the requested composition. The condition excludes coefficient poles and degree drops. For explicit components or all degree strata use `parameters.polynomial_pair` or `parameters.polynomial_family`. Requires the real-domain dependency.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_api.py#L74)

### `real_equivalent`

```text
real_equivalent(chain, **options)
```

Decide equality of a chain's selected real graph on a declared domain.

`options` are those of `real.verify_chain`. Equality of graphs is independent of totality and uniqueness; use that specialist interface for their individual results. Unsupported exact decisions raise their solver error. Requires the real-domain dependency.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_api.py#L87)

### `DecompositionError`

```text
class DecompositionError(ValueError)
```

Invalid decomposition input or a failed exact reconstruction.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_validation.py#L12)

### `DegreeLimitError`

```text
class DegreeLimitError(DecompositionError)
```

The input exceeds a pre-computation degree limit.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_validation.py#L16)

### `IncompleteDecompositionError`

```text
class IncompleteDecompositionError(result)
```

A requested complete mathematical result exceeds the search scope.

The partial search or chain certificate is available as `result`.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/_validation.py#L78)

## `alggeom.decomposition.algebraic`

Implicit composition certificates and supplied-factor reconstruction.

### `AlgebraicDecomposition`

```text
class AlgebraicDecomposition(relation, outer, inner, exponent, unit)
```

An irreducible resultant decomposition of an implicit relation.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/algebraic.py#L28)

Stored fields: `relation`, `outer`, `inner`, `exponent`, `unit`.

#### `AlgebraicDecomposition.output`

Output coordinate of the implicit relation.

#### `AlgebraicDecomposition.input`

Input coordinate of the implicit relation.

#### `AlgebraicDecomposition.intermediate`

Coordinate shared by the two components.

#### `AlgebraicDecomposition.constants`

Declared constant field of the implicit relation.

#### `AlgebraicDecomposition.verify`

```text
verify(self)
```

Check all polynomials, nontriviality, degrees, power and scalar unit.

### `AlgebraicReconstruction`

```text
class AlgebraicReconstruction(relation, ambient, candidate, parameter, coefficient_field, rationality, parametrization, decomposition, reason)
```

Reconstruction from one supplied factor; never a complete discovery claim.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/algebraic.py#L136)

Stored fields: `relation`, `ambient`, `candidate`, `parameter`, `coefficient_field`, `rationality`, `parametrization`, `decomposition`, `reason`.

#### `AlgebraicReconstruction.complete`

False: one supplied factor does not establish discovery completeness.

#### `AlgebraicReconstruction.status`

Whether the supplied factor yields a nontrivial composition.

#### `AlgebraicReconstruction.decompositions`

The reconstructed composition, or an empty tuple for a rejected candidate.

#### `AlgebraicReconstruction.verify`

```text
verify(self)
```

Replay the supplied-factor embedding, rational maps and final identity.

### `reconstruct_algebraic`

```text
reconstruct_algebraic(relation, variables, *, candidate, ambient, parameter=None, constants=QQ, point=None, max_steps=256, max_dimension=256)
```

Reconstruct an implicit composition from a factor in a declared finite field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/algebraic.py#L264)

## `alggeom.decomposition.chains`

Recursive rational and implicit chains with replayable split evidence.

### `RationalComponent`

```text
class RationalComponent(expression, variables)
```

A leaf map; only the last component can have several input coordinates.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/chains.py#L89)

Stored fields: `expression`, `variables`.

### `ChainNode`

```text
class ChainNode(function, variables, constants, search=None, rational_search=None, supplied=None, split=None, outer=None, inner=None, stop_reason=None, detail=None)
```

A certified split, or a leaf with a proof or an explicit stopping reason.

`search` retains the original two-factor result. For implicit rational graphs, `rational_search` records the additional graph-map search. Unsupported diagnostics certify no negative mathematical claim.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/chains.py#L97)

Stored fields: `function`, `variables`, `constants`, `search`, `rational_search`, `supplied`, `split`, `outer`, `inner`, `stop_reason`, `detail`.

#### `ChainNode.leaves`

Leaves in outer-to-inner composition order.

#### `ChainNode.resolved`

Whether this leaf has an exhaustive negative proof in the search scope.

### `RationalChain`

```text
class RationalChain(root, max_depth, max_nodes, max_degree)
```

One ordered outer-to-inner chain over the declared constant field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/chains.py#L335)

Stored fields: `root`, `max_depth`, `max_nodes`, `max_degree`.

#### `RationalChain.original_expression`

Input expression before rational cancellation.

#### `RationalChain.variables`

Ordered input coordinates, or (output,input) for an implicit chain.

#### `RationalChain.constants`

Declared constant field.

#### `RationalChain.leaves`

Leaves in outer-to-inner composition order.

#### `RationalChain.components`

Component maps or polynomials in outer-to-inner order.

#### `RationalChain.complete`

Whether every leaf is resolved in the declared search scope.

#### `RationalChain.status`

Complete or incomplete according to leaf evidence.

#### `RationalChain.input_guard`

Nonvanishing conditions for visible original denominators.

#### `RationalChain.compose`

```text
compose(self)
```

Return the function-field composition, without pointwise pole claims.

#### `RationalChain.verify`

```text
verify(self)
```

Replay every split/search and the final rational composition.

### `AlgebraicChain`

```text
class AlgebraicChain(root, max_depth, max_nodes, max_degree, max_rational_degree)
```

One implicit chain with exact resultant powers and unresolved leaves.

Completeness inherits the minimal coefficient-field scope of automatic algebraic search, augmented by rational graph-map decomposition. It does not certify arbitrary algebraic-function indecomposability or all chains.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/chains.py#L418)

Stored fields: `root`, `max_depth`, `max_nodes`, `max_degree`, `max_rational_degree`.

#### `AlgebraicChain.relation`

The normalized input relation.

#### `AlgebraicChain.constants`

Declared constant field.

#### `AlgebraicChain.variables`

Ordered input coordinates, or (output,input) for an implicit chain.

#### `AlgebraicChain.leaves`

Leaves in outer-to-inner composition order.

#### `AlgebraicChain.components`

Component maps or polynomials in outer-to-inner order.

#### `AlgebraicChain.complete`

Whether every leaf is resolved in the declared search scope.

#### `AlgebraicChain.status`

Complete or incomplete according to leaf evidence.

#### `AlgebraicChain.resultant`

```text
resultant(self)
```

Eliminate the ordered chain's internal coordinates, inner first.

#### `AlgebraicChain.exponent`

Power of the input relation in the iterated resultant.

#### `AlgebraicChain.unit`

Constant scalar multiplying that relation power.

#### `AlgebraicChain.verify`

```text
verify(self)
```

Replay the tree and its global resultant = unit * relation**exponent.

### `rational_chain`

```text
rational_chain(expression, variables, *, constants=QQ, initial=None, max_depth=8, max_nodes=63, max_candidates=256, max_degree=32)
```

Recursively choose the first certified split, for one or several inputs.

`max_nodes` bounds the entire stored binary tree, `max_depth` bounds recursive edges. Two-factor budgets apply separately at each visited node; neither factorization nor certificate replay has a wall-clock bound. A valid incomplete chain remains usable as an exact composition certificate. `initial` optionally supplies a verified two-component certificate at the root; recursive searches then refine both of its components.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/chains.py#L652)

### `algebraic_chain`

```text
algebraic_chain(relation, variables, *, constants=QQ, initial=None, max_depth=8, max_nodes=63, max_candidates=64, max_degree=6, max_rational_degree=32, max_factor_dimension=32, max_separating_attempts=64, max_steps=256, max_dimension=256)
```

Build one recursive implicit chain, preserving resultant powers and units.

Linear-coordinate relations first use rational decomposition in their graph orientation. Other splits use automatic algebraic discovery. Repeated ancestor relations are open leaves, never indecomposability proofs. Completeness has the two-factor engine's declared coefficient-field scope. `initial` accepts an already verified algebraic decomposition, including results reconstructed from a supplied factor and extension field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/chains.py#L695)

## `alggeom.decomposition.common`

A shared rational intermediate coordinate for a tuple of functions.

### `CommonDecomposition`

```text
class CommonDecomposition(functions, field)
```

Constructive equality K(f_1,...,f_m)=K(h), including negative field evidence.

Individual outer maps may be linear. The inner map generates exactly the joint field, so a degree-one inner map is a valid answer. Constant inputs have no distinguished nonconstant generator; higher transcendence degree rules out a common scalar rational intermediate.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/common.py#L13)

Stored fields: `functions`, `field`.

#### `CommonDecomposition.status`

Rational, constant, or higher_transcendence_degree.

#### `CommonDecomposition.inner`

A SymPy Lambda generating the joint field, or None when none is selected.

#### `CommonDecomposition.outers`

The univariate Lambda maps f_i=g_i(h), or an empty tuple without a generator.

#### `CommonDecomposition.input_guard`

Denominator exclusions visible in the original input functions.

#### `CommonDecomposition.verify`

```text
verify(self)
```

Replay field generation, inverse coordinates and each original function.

### `common_components`

```text
common_components(functions, variables, *, parameter=None, constants=QQ)
```

Compute all outer maps through one generator of the joint rational function field.

Equality is in the function field; input_guard records original affine poles. This is a generator of the joint field, not an enumeration of all intermediate fields or all possible shared substitutions.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/common.py#L66)

## `alggeom.decomposition.discovery`

Automatic algebraic decomposition and candidate coverage accounting.

### `CandidateAttempt`

```text
class CandidateAttempt(index, result, error=None)
```

A reconstruction result or an explicitly unresolved candidate diagnostic.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/discovery.py#L23)

Stored fields: `index`, `result`, `error`.

### `AlgebraicSearch`

```text
class AlgebraicSearch(discovery, parameter, attempts)
```

Automatic minimal coefficient-field search, with exact coverage accounting.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/discovery.py#L32)

Stored fields: `discovery`, `parameter`, `attempts`.

#### `AlgebraicSearch.relation`

The normalized input relation.

#### `AlgebraicSearch.candidate_count`

Number of factor orbits discovered for reconstruction.

#### `AlgebraicSearch.examined`

Number of candidate reconstruction attempts retained.

#### `AlgebraicSearch.decompositions`

Every successful composition among the retained attempts.

#### `AlgebraicSearch.complete`

Whether every relevant factor degree and orbit has a resolved reconstruction.

#### `AlgebraicSearch.status`

Decomposability within the declared minimal coefficient-field scope.

#### `AlgebraicSearch.verify`

```text
verify(self)
```

Replay orbit coverage and every claimed answer/rejection.

Unresolved diagnostic text is not a proof of unsupportedness; it always prevents completeness and an indecomposability claim.

### `search_algebraic`

```text
search_algebraic(relation, variables, *, parameter=None, constants=QQ, factor_degrees=None, max_degree=6, max_factor_dimension=32, max_separating_attempts=64, max_candidates=64, exhaustive=False, max_steps=256, max_dimension=256)
```

Search factor orbits for minimal coefficient-field compositions, retaining coverage evidence.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/discovery.py#L121)

## `alggeom.decomposition.factors`

Finite coefficient algebras and exact Galois orbits of polynomial factors.

### `FactorOrbit`

```text
class FactorOrbit(relation, degree, ambient, coefficients, cofactor)
```

One monic-factor orbit, in its explicitly constructed coefficient field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/factors.py#L71)

Stored fields: `relation`, `degree`, `ambient`, `coefficients`, `cofactor`.

#### `FactorOrbit.factor`

The monic factor represented by this orbit.

#### `FactorOrbit.extension_degree`

Degree of the factor coefficient field over constants(input).

#### `FactorOrbit.verify`

```text
verify(self)
```

Replay field embeddings, division identities and claimed orbit coverage.

### `FactorDegree`

```text
class FactorDegree(relation, degree, quotient, parameter, weights, defining_polynomial, coordinate_polynomials, component_polynomials, orbits)
```

An exhaustive, cyclic presentation of all factors of one fixed degree.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/factors.py#L122)

Stored fields: `relation`, `degree`, `quotient`, `parameter`, `weights`, `defining_polynomial`, `coordinate_polynomials`, `component_polynomials`, `orbits`.

#### `FactorDegree.geometric_factor_count`

Number of factors of this degree over an algebraic closure.

#### `FactorDegree.verify`

```text
verify(self)
```

Replay universal equations, inverse maps and irreducible orbit coverage.

### `FactorIssue`

```text
class FactorIssue(degree, reason)
```

An uncompleted factor degree, never a negative decomposition proof.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/factors.py#L243)

Stored fields: `degree`, `reason`.

### `FactorDiscovery`

```text
class FactorDiscovery(relation, requested_degrees, levels, issues, max_factor_dimension, max_separating_attempts)
```

Factor extensions for requested degrees, grouped by Galois orbit.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/factors.py#L251)

Stored fields: `relation`, `requested_degrees`, `levels`, `issues`, `max_factor_dimension`, `max_separating_attempts`.

#### `FactorDiscovery.orbits`

All discovered factor orbits in requested degree order.

#### `FactorDiscovery.complete`

Coverage of requested factor degrees, including every conjugate orbit.

#### `FactorDiscovery.verify`

```text
verify(self)
```

Replay field embeddings, division identities and claimed orbit coverage.

### `discover_factors`

```text
discover_factors(relation, variables, *, degrees=None, constants=QQ, max_degree=6, max_factor_dimension=32, max_separating_attempts=64)
```

Construct coefficient extensions and all monic factor orbits of requested degrees.

Completeness is over constants(input), up to its Galois action, not numerical root selection. Limits bound input degree, factor-algebra dimension and the search for a separating form; they are not wall-clock limits.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/factors.py#L400)

## `alggeom.decomposition.parameters`

Polynomial decomposition on exact real parameter strata.

### `PolynomialSystem`

```text
class PolynomialSystem(polynomial, variable, parameters, degrees, component_variable, outer, inner, unknowns, equations, guard)
```

Normalized coefficient equations for fixed outer and inner degrees.

`inner` is monic with constant zero. `outer` has the original leading and constant coefficients. This fixes affine changes of the middle variable. `formula` is existential in `unknowns` and free in `parameters`.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L42)

Stored fields: `polynomial`, `variable`, `parameters`, `degrees`, `component_variable`, `outer`, `inner`, `unknowns`, `equations`, `guard`.

#### `PolynomialSystem.formula`

Exact normalized coefficient relation, including degree and pole guards.

### `PolynomialPair`

```text
class PolynomialPair(system, condition, outer, inner, coefficient_values, method='normalized_coefficient_reconstruction')
```

Necessary and sufficient condition and explicit normalized components.

Components are valid only on `condition`. The coefficient construction is deterministic in characteristic zero; no parameter sampling is used.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L68)

Stored fields: `system`, `condition`, `outer`, `inner`, `coefficient_values`, `method`.

#### `PolynomialPair.composition`

Candidate composition, to be interpreted under `condition`.

### `PolynomialCase`

```text
class PolynomialCase(degree, polynomial, components, status)
```

All normalized fixed-degree components on one guarded degree stratum.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L285)

Stored fields: `degree`, `polynomial`, `components`, `status`.

#### `PolynomialCase.subs`

```text
subs(self, assignments)
```

Specialize symbolic components after selecting their parameter guard.

### `PolynomialFamily`

```text
class PolynomialFamily(polynomial, variable, parameters, component_variable, parameter_domain, stratification)
```

A complete guarded partition, including exceptional degrees and poles.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L325)

Stored fields: `polynomial`, `variable`, `parameters`, `component_variable`, `parameter_domain`, `stratification`.

#### `PolynomialFamily.branches`

Disjoint guarded cases covering the declared parameter domain.

#### `PolynomialFamily.select`

```text
select(self, assignments)
```

Specialize the unique case containing an exact parameter assignment.

#### `PolynomialFamily.verify`

```text
verify(self)
```

Replay the parameter partition and all guarded component identities.

### `polynomial_system`

```text
polynomial_system(polynomial, variable, *, degrees, parameters=(), component_variable=None)
```

Build exact coefficient equations for `polynomial = outer(inner(x))`.

`degrees=(outer_degree, inner_degree)` must contain integers at least two. Coefficients may be rational functions of the declared real parameters with exact real algebraic constants. Leading-coefficient and pole guards exclude specializations with smaller degree or undefined input coefficients.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L89)

### `polynomial_pair`

```text
polynomial_pair(polynomial, variable, *, degrees, parameters=(), component_variable=None)
```

Return a real parameter condition and explicit fixed-degree components.

Normalize the inner polynomial to be monic with constant zero. Its remaining coefficients are uniquely determined, top down, by the highest coefficients of the input. Then recover the outer coefficients top down. Vanishing of the residual coefficients is necessary and sufficient on the degree/pole guard. This avoids an unnormalized nonlinear QE search and works in characteristic zero. It does not search rational or general implicit algebraic components.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L165)

### `polynomial_family`

```text
polynomial_family(polynomial, variable, *, parameters=(), parameter_domain=True, component_variable=None)
```

Reconstruct normalized components across every exact degree-drop stratum.

Branches include decomposable, indecomposable, degree-below-two, and undefined coefficient-pole cases. Every ordered nontrivial factorization of the actual degree is considered. Guards are exact formulas, never sample substitutions. Empty guards may remain symbolically; they assert no result at any real point.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L355)

### `verify_family`

```text
verify_family(result)
```

Replay branch coverage, actual degrees, component guards and identities.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/parameters.py#L454)

## `alggeom.decomposition.rational`

Exact rational compositions over one or several input coordinates.

### `RationalDecomposition`

```text
class RationalDecomposition(original_expression, variables, parameter, constants, outer, inner)
```

An exact f(X)=outer(inner(X)), with a univariate outer function.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/rational.py#L134)

Stored fields: `original_expression`, `variables`, `parameter`, `constants`, `outer`, `inner`.

#### `RationalDecomposition.function`

The reduced rational function over the declared constants.

#### `RationalDecomposition.degrees`

Univariate outer degree and rational total degree of the inner.

#### `RationalDecomposition.input_guard`

Nonvanishing conditions for denominators visible in the original expression.

#### `RationalDecomposition.verify`

```text
verify(self)
```

Check exact composition, declared coordinates and nontrivial component degrees.

### `RationalSearch`

```text
class RationalSearch(original_expression, variables, parameter, constants, comparison_variables, factors, examined, candidate_count, decompositions)
```

A searched prefix of all two-factor classes over the declared constants.

Completeness is modulo Möbius changes of intermediate coordinate, with inner total degree >1 and outer degree >1. It does not assert recursive chains, larger constant fields, or other multivariate decomposition notions.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/rational.py#L210)

Stored fields: `original_expression`, `variables`, `parameter`, `constants`, `comparison_variables`, `factors`, `examined`, `candidate_count`, `decompositions`.

#### `RationalSearch.complete`

Whether every candidate in the declared two-factor scope was examined.

#### `RationalSearch.status`

Decomposability or the precise reason no conclusion is available.

#### `RationalSearch.verify`

```text
verify(self)
```

Replay exact component identities or the searched candidate prefix.

### `search_rational`

```text
search_rational(expression, variables, *, parameter=None, constants=QQ, max_candidates=256, max_degree=32)
```

Find GRS02 uni-multivariate two-factor decompositions up to a budget.

Degrees are rational total degrees in the input variables. Candidate products have equal proper-divisor degrees in both variable blocks and vanish on the diagonal. Counting is polynomial in the degree bound and number of factors; enumeration is lazy. Factorization and wall-clock time are not bounded by `max_candidates`. Unsupported exact operations raise.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/rational.py#L305)

## `alggeom.decomposition.real`

Exact real graph comparison, branch selection, totality and uniqueness.

### `verify_chain`

```text
verify_chain(chain, *, output=None, intermediates=None, domain=True, target_selector=True, component_selectors=None, check_functions=True)
```

Compare a decomposition chain with its selected affine real composition.

Components and selectors are ordered outermost to innermost. Rational chains may name their intermediate coordinates explicitly; implicit chains already fix those coordinates in `chain.components`. Each selector may use only its component's input and output coordinates. Chain search need not be complete, but its exact composition certificate must verify.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/real.py#L47)

### `verify_decomposition`

```text
verify_decomposition(decomposition, *, output=None, domain=True, target_selector=True, outer_selector=True, inner_selector=True, check_functions=True)
```

Compare a two-component certificate on a selected real input domain.

Rational certificates require a separate output coordinate when selectors mention it. Their intermediate coordinate is `decomposition.parameter`. Implicit certificates fix all three coordinates in their polynomials. Selectors use only the input and output of the corresponding component.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/decomposition/real.py#L154)

## `alggeom.families`

Parameter-sensitive invariants of a univariate polynomial fiber family.

### `FiberInvariants`

```text
class FiberInvariants(dimension, geometric_points)
```

Geometric dimension and distinct-point count of one fiber over an algebraic closure.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/families.py#L8)

Stored fields: `dimension`, `geometric_points`.

### `FiberStratum`

```text
class FiberStratum(condition, invariants)
```

A parameter condition with exact geometric fiber invariants.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/families.py#L16)

Stored fields: `condition`, `invariants`.

### `FiberStratification`

```text
class FiberStratification(equation, variable, parameter, exceptional_polynomial, generic)
```

Generic invariants and a complete finite exceptional locus in A^1.

The complement of `exceptional_polynomial=0` is a certified open stratum. `strata()` isolates its distinct algebraic zeros and computes their exact fiber invariants, covering all geometric parameter points.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/families.py#L24)

Stored fields: `equation`, `variable`, `parameter`, `exceptional_polynomial`, `generic`.

#### `FiberStratification.generic_condition`

A sufficient Zariski-open condition for the generic invariants.

#### `FiberStratification.strata`

```text
strata(self, *, max_exceptional_roots=32, max_root_degree=2)
```

Enumerate a full stratification of the one-dimensional parameter line.

Over an algebraic closure the exceptional polynomial has finitely many roots. Each distinct root forms one closed stratum, with its fiber invariants computed by exact specialization. The complement is one open stratum. A declared root-count limit prevents accidental large algebraic isolations. The default degree limit excludes higher-degree RootOf specializations that can be prohibitively expensive.

#### `FiberStratification.specialize`

```text
specialize(self, value)
```

Compute the geometric dimension and number of distinct points exactly.

### `fiber_stratification`

```text
fiber_stratification(equation, variable, parameter)
```

Stratify univariate polynomial or multidimensional affine-linear fibers.

A sequence of equations and sequences of fiber/parameter symbols select the determinantal algorithm, which supports any finite number of variables and polynomial parameters but requires linear fiber equations. Scalar arguments select the discriminant-based univariate algorithm.

The discriminant of the generic squarefree polynomial, its coefficient denominators and leading coefficient, and the content of the original family determine a certified open set with constant geometric root count. This avoids confusing repeated generic roots with specialization failures.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/families.py#L87)

## `alggeom.fields`

Prime and finite-extension coefficient fields for exact algebraic geometry.

### `AlgebraicFunctionField`

```text
class AlgebraicFunctionField(modulus)
```

A certified finite field extension of QQ(parameters) or K(parameters).

Construct from an irreducible univariate `Poly` over QQ or a rational function field over QQ/a number field. `extend` accepts a modulus over this field, checks irreducibility through a squarefree norm, and compresses the result. Generator names remain valid coefficient expressions. A symbol denotes a formal algebraic root; real or complex embeddings are chosen explicitly when specializing.

[Source](https://github.com/BhuvaneshBhatt/algroots/blob/main/src/algroots/algebraic_function_fields.py#L103)

### `finite_field`

```text
finite_field(characteristic, degree=1, *, name='a', modulus=None)
```

Return GF(p) or GF(p^degree) with a validated irreducible defining polynomial.

A supplied modulus determines the extension; otherwise python-flint chooses one. The returned object is a SymPy coefficient domain. GF(p^degree) is never interpreted as the composite residue ring Z/(p^degree).

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/fields.py#L8)

### `field_order`

```text
field_order(field)
```

Return the cardinality of a supported finite coefficient field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/fields.py#L18)

## `alggeom.function_field_constants`

Certified algebraic constants of finitely generated function fields.

For a prime ideal P over QQ, a rational function p/q in Frac(QQ[x]/P) may satisfy a polynomial over QQ. Its minimal polynomial is obtained by eliminating the coordinates from the saturated graph ideal (P, t*q-p):q^infinity. A proposed constant field is accepted as complete only when every conjugate graph component is certified geometrically prime and their intersection equals the extended original ideal.

For a prime hypersurface quadratic in one coordinate, the square class of its discriminant gives a complete, terminating constant-field algorithm: nonconstant squarefree factors prove geometric integrality, while a constant nonsquare square class gives a quadratic constant field. A standalone proof replays the discriminant identity and, in the split case, the saturated graph ideals. Other presentations still have a bounded candidate search; failure to discover a suitable primitive constant never means the constant field is QQ.

### `AlgebraicConstants`

```text
class AlgebraicConstants(field, splitting_field, minimal_polynomial, numerator, denominator, geometric_primes, certificate=None)
```

A certified full constant field and its geometric generic components.

The fraction `numerator/denominator` represents a primitive algebraic constant in the original function field. Its `minimal_polynomial` is irreducible over QQ. Every `geometric_primes` entry is proven absolutely prime over `field`; exact intersection reconstructs the base-changed input ideal. When the constants are QQ, numerator is zero and denominator is one.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/function_field_constants.py#L27)

Stored fields: `field`, `splitting_field`, `minimal_polynomial`, `numerator`, `denominator`, `geometric_primes`, `certificate`.

#### `AlgebraicConstants.verify`

```text
verify(self, prime)
```

Replay a complete discriminant certificate, when available.

Discriminant-derived results carry an independently replayable certificate. Search results without that evidence return False.

### `QuadraticConstantCertificate`

```text
class QuadraticConstantCertificate(coordinate, discriminant, square_class, square_factor)
```

Recheckable completeness proof for a quadratic function field over QQ.

If `square_factor` is None, the discriminant has a nonconstant factor with odd multiplicity, proving irreducibility in QQbar(t)[y]. Otherwise `discriminant = square_class * square_factor**2`. Since the defining prime has degree two over QQ(t), its full constant field is QQ(sqrt(c)). All divisions and component graph constructions are independently reconstructed by the verifier; the producer's Gröbner basis is not trusted.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/function_field_constants.py#L56)

Stored fields: `coordinate`, `discriminant`, `square_class`, `square_factor`.

#### `QuadraticConstantCertificate.verify`

```text
verify(self, prime, result)
```

Independently recheck discriminant, full constants, and geometric graphs.

### `algebraic_constant_field`

```text
algebraic_constant_field(prime, *, max_degree=4, numerator_degree=1, max_candidates=60, max_points=200)
```

Certify the full algebraic constant field of a supported QQ function field.

Searches low-degree rational functions for a primitive algebraic constant. A candidate is accepted only when all conjugate graph ideals are certified geometrically prime and their intersection recovers the extended prime. Hence the result describes the *full* algebraic closure of QQ in the function field, not merely a known algebraic subfield.

A geometrically prime field has constants QQ. Unsupported primes or a search that does not discover a full primitive constant raise `NotImplementedError`. The bounded search is not a general decision procedure for all finitely generated function fields.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/function_field_constants.py#L424)

### `certify_geometric_prime`

```text
certify_geometric_prime(ideal, *, max_points=200)
```

Prove a prime ideal stays prime over the algebraic closure.

Supports arbitrary-dimensional nonlinear ideals when they have a smooth rational point and are provably prime over the coefficient field. Also recognizes linear and polynomial-graph primes. Raises NotImplementedError when the available criteria cannot certify geometric primality; it does not treat an inconclusive search as evidence of reducibility.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/function_field_constants.py#L325)

## `alggeom.geometric`

Geometric integrality and exact primary factors over algebraic extensions.

Absolute algorithms handle hypersurfaces, bounded finite algebras, and complete intersections in disjoint variable blocks over QQ. Unrestricted positive-dimensional nonprincipal splitting is not supplied.

### `GeometricComponents`

```text
class GeometricComponents(components, field)
```

Absolute components explicitly represented over a splitting number field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/geometric.py#L16)

Stored fields: `components`, `field`.

### `AbsolutePrimaryDecomposition`

```text
class AbsolutePrimaryDecomposition(components, field)
```

Geometric primary components after an exact finite algebraic field extension.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/geometric.py#L209)

Stored fields: `components`, `field`.

#### `AbsolutePrimaryDecomposition.verify`

```text
verify(self, original, **limits)
```

Independently verify primary status, absolute supports and reconstruction.

Every associated prime of each component is recomputed over the declared splitting field and must equal that component's radical. Each radical must have its own geometric-primality certificate. Finally, the intersection must exactly recover the extended input, and no component may be redundant. Embedded primes are retained. No computed certificate, candidate Groebner backend, or component classification flag is trusted here.

### `is_geometrically_integral`

```text
is_geometrically_integral(ideal)
```

Decide geometric integrality for principal linear/quadratic hypersurfaces.

For an affine plane conic, homogenizing the equation gives a ternary quadratic form. Over an algebraic closure its zero set is irreducible exactly when that form has rank three. A univariate polynomial describes one geometrically integral point precisely when its degree is one.

Supported finite zero-dimensional ideals are also tested by splitting their geometric primary components and checking the reduced quotient.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/geometric.py#L32)

### `geometric_components`

```text
geometric_components(ideal, *, extension=None)
```

Return certified geometric components of supported finite ideals and hypersurfaces.

Quadratic splitting fields over QQ are inferred from a discriminant. An explicit extension may also be supplied. Every returned factor is checked to have degree one, preventing an incomplete decomposition.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/geometric.py#L97)

### `absolute_primary_decomposition`

```text
absolute_primary_decomposition(ideal, *, max_root_degree=4, max_components=128)
```

Compute exact absolute primary factors when geometric supports can be certified.

A univariate polynomial over QQ splits into powers of linear maximal ideals over its constructed splitting field; all multiplicities are preserved. Principal ideals use exact absolute hypersurface factorization. Certain nonprincipal ideals are handled through disjoint-variable complete intersections, monomial decomposition, and certified geometric supports. The max_components bound limits the number of separated intersections. General nonprincipal absolute primary decomposition remains unsupported.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/geometric.py#L768)

## `alggeom.graded`

Standard-graded Hilbert invariants, dimensions and projective closure.

### `is_homogeneous`

```text
is_homogeneous(ideal)
```

Decide homogeneity for the standard grading, treating coefficient-field elements as scalars.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L33)

### `hilbert_series`

```text
hilbert_series(ideal, *, parameter=None)
```

Return the exact rational Hilbert series of S/I, not of I as an S-module.

The standard grading assigns degree one to each coordinate. The numerator is computed on the leading monomial ideal by lcm inclusion-exclusion. Arbitrary weights and multigradings are not inferred.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L45)

### `hilbert_function`

```text
hilbert_function(ideal, degree)
```

Return dim_K(S/I)_degree for a nonnegative integer degree.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L59)

### `hilbert_polynomial`

```text
hilbert_polynomial(ideal, *, parameter=None)
```

Return the eventual polynomial of the graded-piece dimensions as a QQ Poly.

This is not the cumulative Hilbert-Samuel polynomial. An Artinian graded quotient has a finite Hilbert series and eventual Hilbert polynomial zero.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L76)

### `krull_dimension`

```text
krull_dimension(ideal)
```

Return dim(S/I), also for inhomogeneous ideals; -1 denotes the zero ring.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L97)

### `projective_dimension`

```text
projective_dimension(ideal)
```

Return dim Proj(S/I), with -1 for an empty projective scheme.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L102)

### `degree`

```text
degree(ideal)
```

Return standard-graded Hilbert multiplicity; an Artinian quotient gives its length.

For a nonempty projective scheme this is its degree. The zero ring gives zero. No smoothness, radicality or equidimensionality is assumed.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L109)

### `projective_closure`

```text
projective_closure(ideal, *, homogenizer=None)
```

Return the homogenized ideal saturated at a new projective coordinate.

Saturation removes components at infinity introduced by homogenizing only the supplied affine generators. The resulting ring appends the homogenizer to the original coordinates. Scheme structure is retained.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/graded.py#L130)

## `alggeom.groebner_certificates`

Exact certificates for polynomial bases supplied by algroots.

Basis construction and membership-witness production belong to algroots. This module checks the resulting ideal identities, reducedness and S-pair remainders before a certificate is used in geometric calculations.

### `GroebnerBasisProof`

```text
class GroebnerBasisProof(basis, memberships, compact_proof=None)
```

A reduced basis with exact generator-combination evidence.

Each memberships row expresses one basis element using the original generators. Coefficients may be polynomial expressions or sparse rational maps. A compact_proof stores a rational arithmetic DAG instead of rows. Supply one representation and call verify to check its claims.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/groebner_certificates.py#L104)

Stored fields: `basis`, `memberships`, `compact_proof`.

#### `GroebnerBasisProof.verify`

```text
verify(self, generators, ring, *, normal_form_backend='sympy')
```

Replay the proof against the supplied original generators and ring.

### `construct_membership_witnesses`

```text
construct_membership_witnesses(generators, candidate, ring, *, max_basis=128, max_pairs=500)
```

Construct exact expression witnesses for a candidate's ideal membership.

First attempt division by input generators, which is often inexpensive. If it fails, run a targeted, tracked Buchberger process with leading-pair priority and Buchberger product criterion. Stop as soon as all requested membership identities are available, without completing an unnecessary second Gröbner basis. Proof production can still be expensive. `max_basis` and `max_pairs` are hard safety budgets, not mathematical completeness bounds; a budget failure raises without claiming a proof. The resulting witnesses must still be independently checked by `verify`.

[Source](https://github.com/BhuvaneshBhatt/algroots/blob/main/src/algroots/groebner_proofs.py#L121)

### `verify_reduced_groebner_basis`

```text
verify_reduced_groebner_basis(generators, candidate, ring, *, membership_witnesses=None, compact_proof=None, normal_form_backend='sympy')
```

Verify a reduced monic basis independently of its construction backend.

Membership witnesses have one row per candidate and one column per input generator. Exact replay establishes candidate inclusion in the input ideal; without witnesses, an independent reference basis checks this inclusion. Input reduction establishes the reverse inclusion. Reducedness and every required S-pair remainder are checked in both paths.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/groebner_certificates.py#L132)

## `alggeom.groebner_cover`

Certified parameter-specialized Groebner systems over exact coefficient fields.

The computation branches on leading coefficients encountered by Buchberger's algorithm. An open branch retains a basis whose reductions remain valid after specialization. Exceptional fibers are recomputed, not evaluated from a potentially invalid generic basis. One-parameter families admit a complete finite decomposition when the exceptional algebraic roots can be represented. For multiple parameters, exceptional divisors with polynomial algebraic parameter-graph charts can be recursively split. More general implicit parameter varieties raise an explicit coverage error.

### `GroebnerCover`

```text
class GroebnerCover(equations, variables, parameters, pieces)
```

Certified basis selection for parameter-specialized polynomial ideals.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/groebner_cover.py#L34)

Stored fields: `equations`, `variables`, `parameters`, `pieces`.

#### `GroebnerCover.verify_partition`

```text
verify_partition(self)
```

Prove exact coverage and pairwise disjointness of parameter strata.

#### `GroebnerCover.certificate`

```text
certificate(self)
```

Produce replayable ideal-membership and Buchberger certificates.

Each proof stores polynomial reduction quotients, S-pair reductions, and coefficient-denominator guards checked on its entire stratum.

#### `GroebnerCover.specialize`

```text
specialize(self, values)
```

Return a specialized basis at an exact algebraic parameter tuple.

#### `GroebnerCover.verify_at`

```text
verify_at(self, values)
```

Independently compare the selected basis with a specialized Groebner basis.

### `GroebnerCoverPiece`

```text
class GroebnerCoverPiece(condition, basis, order, parameter_bindings=(), free_parameters=(), algebraic_constants=())
```

A constructible parameter condition and a valid Groebner basis on it.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/groebner_cover.py#L22)

Stored fields: `condition`, `basis`, `order`, `parameter_bindings`, `free_parameters`, `algebraic_constants`.

### `GroebnerCoverLimitError`

```text
class GroebnerCoverLimitError(NotImplementedError)
```

The parameter branching requires an unsupported algebraic specialization.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/groebner_cover.py#L18)

### `comprehensive_groebner_cover`

```text
comprehensive_groebner_cover(equations, variables, parameters, *, max_pairs=150, max_pieces=64, max_exceptional_degree=6)
```

Build specialization-valid polynomial bases for tractable parameter strata.

A complete answer is returned for univariate parameter spaces when all exceptional roots admit exact representations, and for multivariate parameter spaces when all exceptional branches admit polynomial parameter charts over an algebraic constant extension, or a supported generic finite extension. Unsupported implicit parameter varieties raise instead of returning an incomplete cover.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/groebner_cover.py#L218)

## `alggeom.ideals`

Exact polynomial ideals and coordinate rings over supported exact coefficient fields.

### `PolynomialRing`

```text
class PolynomialRing(variables, domain=QQ, order='grevlex')
```

An ordered polynomial ring over a supported exact coefficient field.

Variables may be empty, representing the coefficient field itself. The term order controls normal forms; it is not part of mathematical ideal equality. Coefficient parameters and polynomial variables must be disjoint.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/ideals.py#L15)

Stored fields: `variables`, `domain`, `order`.

#### `PolynomialRing.convert`

```text
convert(self, expression)
```

Return the canonical exact expression, rejecting poles in the ring variables.

#### `PolynomialRing.ideal`

```text
ideal(self, *generators, groebner_backend='native', groebner_limits=None)
```

Generate an ideal with an explicit verified Gröbner construction policy.

Accept expressions or one list/tuple. Backend selection is stored on this ideal and controls its cached basis and normal forms.

#### `PolynomialRing.quotient`

```text
quotient(self, *generators, groebner_backend='native', groebner_limits=None)
```

Return a coordinate ring with the selected verified basis backend.

#### `PolynomialRing.with_order`

```text
with_order(self, order)
```

Use another term order with the same coefficient field and coordinates.

### `Ideal`

```text
class Ideal(ring, generators=(), groebner_backend='native', groebner_limits=None)
```

A finitely generated ideal, with mathematical equality and exact normal forms.

Zero generators are omitted; the empty generator tuple is the zero ideal. Membership witnesses refer to the stored generators, not only a computed Groebner basis. Equality requires identical coordinates and coefficients, but permits different term orders.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/ideals.py#L88)

Stored fields: `ring`, `generators`, `groebner_backend`, `groebner_limits`.

#### `Ideal.groebner_basis_from`

```text
groebner_basis_from(self, backend='native', *, limits=None)
```

Return a reduced basis constructed and verified through algroots.

Native requests over a prime field run directly in that characteristic. Rational requests use modular lifting. Every accelerated result passes an independent exact check before it enters geometric calculations.

#### `Ideal.groebner_basis_proof_from`

```text
groebner_basis_proof_from(self, backend='native', *, limits=None, max_witness_basis=128, max_witness_pairs=500)
```

Return a basis with independently checked generator-combination witnesses.

Algroots constructs the candidate and its membership witnesses. Exact replay in alggeom establishes both ideal inclusions before the proof is returned. Use `backend="sympy"` for a portable candidate calculation.

#### `Ideal.groebner_basis`

The reduced monic Groebner basis as a tuple of expressions.

#### `Ideal.is_zero`

Whether the ideal is zero.

#### `Ideal.is_unit`

Whether the ideal is the entire polynomial ring.

#### `Ideal.krull_dimension`

The dimension of the quotient S/I; -1 denotes the zero ring.

#### `Ideal.divide`

```text
divide(self, expression)
```

Return quotients against groebner_basis and the unique reduced remainder.

#### `Ideal.normal_form`

```text
normal_form(self, expression)
```

Return the reduced polynomial representative modulo this ideal.

#### `Ideal.contains`

```text
contains(self, expression)
```

Decide exact ideal membership; invalid coefficient/ring inputs raise.

#### `Ideal.is_subset`

```text
is_subset(self, other)
```

Decide containment in another ideal of the same ambient ring.

#### `Ideal.radical_contains`

```text
radical_contains(self, expression)
```

Decide f in radical(I) by the exact saturation criterion I:f^infinity=S.

#### `Ideal.radical`

```text
radical(self, **limits)
```

Return the radical over the declared coefficient field.

#### `Ideal.minimal_primes`

```text
minimal_primes(self, **limits)
```

Return the prime ideals of the reduced irreducible components.

#### `Ideal.associated_primes`

```text
associated_primes(self, **limits)
```

Return all associated primes, including embedded primes.

#### `Ideal.primary_decomposition`

```text
primary_decomposition(self, **limits)
```

Return an irredundant primary decomposition retaining embedded components.

#### `Ideal.membership`

```text
membership(self, expression)
```

Return coefficients expressing membership in the stored generators, or None.

#### `Ideal.change_ring`

```text
change_ring(self, ring)
```

Embed the generators in another explicitly supplied polynomial ring.

#### `Ideal.eliminate`

```text
eliminate(self, *variables)
```

Contract to the ring of retained coordinates by lexicographic elimination.

Eliminating every coordinate returns an ideal in the coefficient field. Elimination describes a Zariski closure, not an exact image of points with inequalities or pole exclusions.

#### `Ideal.intersection`

```text
intersection(self, other)
```

Return the exact intersection using an auxiliary-coordinate elimination.

#### `Ideal.quotient`

```text
quotient(self, other)
```

Return I:J; a zero divisor ideal J=0 gives the unit ideal.

#### `Ideal.saturate`

```text
saturate(self, other)
```

Return I:J^infinity, supporting a polynomial or a finitely generated ideal.

Principal saturation eliminates t from I+(1-t*f). For several generators intersect the principal saturations; saturating by their product would compute a different operation.

### `CoordinateRing`

```text
class CoordinateRing(ideal)
```

The possibly positive-dimensional algebra S/I.

krull_dimension measures chains of prime ideals. vector_space_dimension is the length over the coefficient field and is infinite for a nonempty positive-dimensional algebra. The zero ring has Krull dimension -1 here.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/ideals.py#L399)

Stored fields: `ideal`.

#### `CoordinateRing.ambient`

The ambient polynomial ring S.

#### `CoordinateRing.variables`

The ordered ambient coordinate symbols.

#### `CoordinateRing.domain`

The coefficient field, independent of the quotient's Krull dimension.

#### `CoordinateRing.normal_form`

```text
normal_form(self, expression)
```

Return the canonical representative of a coordinate-ring element.

#### `CoordinateRing.element`

```text
element(self, expression)
```

Construct an element with this parent and a reduced representative.

#### `CoordinateRing.gens`

The residue classes of the ordered ambient coordinates.

#### `CoordinateRing.zero`

The additive identity, including in the zero ring.

#### `CoordinateRing.one`

The multiplicative identity; it equals zero in the zero ring.

#### `CoordinateRing.hom`

```text
hom(self, images, target)
```

Define a coefficient-preserving homomorphism by the coordinate images.

#### `CoordinateRing.quotient`

```text
quotient(self, *generators, groebner_backend='native', groebner_limits=None)
```

Form a further quotient, retaining the ambient polynomial coordinates.

#### `CoordinateRing.base_change`

```text
base_change(self, domain)
```

Extend the coefficient field through its specified exact embedding.

This is field extension, not specialization of transcendental parameters.

#### `CoordinateRing.krull_dimension`

The Krull dimension of S/I, with -1 for the zero ring.

#### `CoordinateRing.vector_space_dimension`

Finite algebra length, or infinity for positive-dimensional algebras.

#### `CoordinateRing.is_reduced`

The field-relative is reduced property of this coordinate ring.

#### `CoordinateRing.is_irreducible`

The field-relative is irreducible property of this coordinate ring.

#### `CoordinateRing.is_integral`

The field-relative is integral property of this coordinate ring.

#### `CoordinateRing.is_equidimensional`

The field-relative is equidimensional property of this coordinate ring.

#### `CoordinateRing.is_unmixed`

The field-relative is unmixed property of this coordinate ring.

#### `CoordinateRing.codimension`

The field-relative codimension property of this coordinate ring.

#### `CoordinateRing.fraction_field`

```text
fraction_field(self)
```

Return the rational-function field of this integral coordinate ring.

#### `CoordinateRing.localize`

```text
localize(self, *denominators)
```

Invert the supplied functions with a parent-aware localization.

#### `CoordinateRing.normalize`

```text
normalize(self, **limits)
```

Return the integral-closure coordinate ring of a reduced affine algebra.

#### `CoordinateRing.is_normal`

Whether this coordinate ring is normal.

#### `CoordinateRing.finite_algebra`

```text
finite_algebra(self, *, max_dimension=None)
```

Convert a finite coordinate algebra to algroots, rejecting positive dimension.

A scalar coefficient field uses a dummy coordinate constrained to zero because the finite-algebra backend requires a nonempty coordinate list.

## `alggeom.linear_strata`

Exact determinantal partitions for parameter-dependent affine-linear systems.

### `LinearFiberPiece`

```text
class LinearFiberPiece(locus, dimension, geometric_points, rank)
```

Constructible parameter locus with fixed affine-fiber dimension.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/linear_strata.py#L28)

Stored fields: `locus`, `dimension`, `geometric_points`, `rank`.

### `LinearFiberStratification`

```text
class LinearFiberStratification(parameters, variables, coefficient_matrix, constants, pieces)
```

Exact rank partition of A(a)x=b(a) over an algebraically closed field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/linear_strata.py#L38)

Stored fields: `parameters`, `variables`, `coefficient_matrix`, `constants`, `pieces`.

#### `LinearFiberStratification.specialize`

```text
specialize(self, values)
```

Return fiber dimension and distinct point count at an exact parameter tuple.

### `linear_fiber_stratification`

```text
linear_fiber_stratification(equations, variables, parameters, *, max_minors=4096)
```

Partition a polynomial-parameter family of affine-linear systems exactly.

The equations must be affine linear in the fiber coordinates, but may be polynomial of arbitrary degree in the parameters. Each returned piece is a finite union of locally closed determinantal conditions. Inconsistent systems have empty fibers; positive-dimensional affine fibers have infinitely many points over the algebraic closure.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/linear_strata.py#L65)

## `alggeom.local`

Tangent cones, local dimensions, multiplicities and isolated intersection lengths.

### `tangent_cone`

```text
tangent_cone(value, point)
```

Return the homogeneous tangent-cone ideal at an exact point.

Dilating coordinates about the point and saturating the deformation parameter computes initial forms of the entire ideal, including forms created by cancellation. Taking only initial terms of input generators would miss such equations. Projective points use their first nonzero chart.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/local.py#L34)

### `local_dimension`

```text
local_dimension(value, point)
```

Return local dimension after extending to the point's coefficient field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/local.py#L53)

### `embedding_dimension`

```text
embedding_dimension(value, point)
```

Return the dimension of the Zariski tangent space at the point.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/local.py#L58)

### `multiplicity`

```text
multiplicity(value, point)
```

Return Hilbert-Samuel multiplicity from the graded tangent cone.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/local.py#L69)

### `module_length`

```text
module_length(module, point=None, *, max_dimension=256, max_vector_dimension=4096)
```

Return the length of a finite-support module, globally or at an exact point.

Global length is vector-space dimension over the coefficient field. At a field-valued point the module is base changed first and its primary local factor is selected. Nonzero-dimensional support is rejected.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/local.py#L74)

### `intersection_multiplicity`

```text
intersection_multiplicity(left, right, point, *, max_dimension=256)
```

Compute a proper isolated local intersection using a Koszul Tor complex.

One input must be a complete intersection in its supplied equations. Higher Tor lengths are included, so the answer is not assumed to equal the length of the sum ideal. General intersections needing a non-Koszul free resolution are outside this interface.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/local.py#L128)

## `alggeom.localized_linear_cover`

Complete, certified parameter covers for linear univariate fiber systems.

Each fiber equation has the form f_i(a) * x + g_i(a), over an arbitrary QQ[a]/I and any explicitly localized basic open D(h).  The recursion splits on coefficient ideals and residual determinants without requiring I prime, reduced, equidimensional, or free from embedded components.  It terminates in at most two coefficient branches per equation plus a disjoint partition of the determinant exceptional ideal, subject to a declared piece budget.

This is a *pointwise* comprehensive cover on geometric parameter points. On nonreduced parameter schemes, equality of geometric fibers does not imply flat or scheme-theoretic constancy across nilpotent infinitesimal directions. Higher-degree fiber polynomials are not covered by this module.

### `LinearCoverLimitError`

```text
class LinearCoverLimitError(RuntimeError)
```

The exact branch or verification budget was exceeded.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/localized_linear_cover.py#L23)

### `LinearCoverPiece`

```text
class LinearCoverPiece(localization, kind, pivot=None, second=None)
```

A disjoint geometric parameter locus and one certified fiber basis.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/localized_linear_cover.py#L38)

Stored fields: `localization`, `kind`, `pivot`, `second`.

#### `LinearCoverPiece.condition`

Exact polynomial equalities and nonvanishing guards for this piece.

### `LocalizedLinearCover`

```text
class LocalizedLinearCover(equations, variable, base, pieces)
```

Recursive complete cover for finitely many degree <= 1 equations in x.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/localized_linear_cover.py#L53)

Stored fields: `equations`, `variable`, `base`, `pieces`.

#### `LocalizedLinearCover.basis`

```text
basis(self, piece)
```

Provide (1), (), or a monic linear polynomial over the localization.

#### `LocalizedLinearCover.verify`

```text
verify(self, *, max_pairs=400)
```

Independently replay every fiber identity and cover partition.

No candidate Groebner backend is invoked. Nonvanishing witnesses are explicit inverted parameters; every vanishing identity is reduced modulo the exact saturation I:S^infinity. Coverage and disjointness use symbolic constructible-set differences and ideal saturation.

#### `LocalizedLinearCover.specialize`

```text
specialize(self, values)
```

Select the exact fiber basis for an algebraic parameter point.

### `comprehensive_localized_linear_cover`

```text
comprehensive_localized_linear_cover(equations, variable, base, *, max_pieces=128)
```

Complete recursion for linear fiber equations over arbitrary QQ[a]/I.

Branches first on the coefficient f_i, then on its constant g_i. If f_i is invertible, the remaining equations reduce to determinant conditions f_i*g_j-f_j*g_i. Splitting their *entire arbitrary ideal* yields disjoint unit branches and one pivot branch. Every leaf is independently verified.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/localized_linear_cover.py#L190)

## `alggeom.maps`

Graphs, image closures and affine geometric diagnostics for exact polynomial ideals.

### `RationalMap`

```text
class RationalMap(source, coordinates, targets, excluded=())
```

Rational expressions on an explicit principal open subset of an affine scheme.

Original denominator exclusions are retained even after cancellation. Graph and image ideals describe Zariski closures. They do not encode an exact constructible image, a real locus, or selected real branches.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L35)

Stored fields: `source`, `coordinates`, `targets`, `excluded`.

#### `RationalMap.functions`

Reduced coordinate expressions over the ambient polynomial fraction field.

#### `RationalMap.denominator_factors`

Polynomial nonvanishing factors defining the declared affine input domain.

#### `RationalMap.domain_guard`

Nonvanishing conditions, independent of the source's defining equations.

#### `RationalMap.graph_ideal`

The closure of the graph over its declared denominator-free input domain.

#### `RationalMap.image_closure`

The defining ideal of the scheme-theoretic image closure in target space.

#### `RationalMap.constructible_image`

Return the exact geometric image, retaining excluded denominator loci.

#### `RationalMap.evaluate`

```text
evaluate(self, point)
```

Evaluate an exact coefficient-field point, checking equations and input poles.

#### `RationalMap.compose`

```text
compose(self, inner, *, targets=None)
```

Compose maps while retaining poles from both stages.

The inner image must lie in the outer source. New output symbols are chosen when the outer targets would capture an inner source coordinate.

### `implicitize`

```text
implicitize(expressions, parameters, *, targets=None, source=None, domain=None)
```

Return the image-closure ideal, inferring the field from source or defaulting to QQ.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L200)

### `graph_ideal`

```text
graph_ideal(expressions, parameters, *, targets=None, source=None, domain=None)
```

Return the saturated graph closure, with source coordinates preceding targets.

The field is inferred from source, or defaults to QQ when no source is given.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L205)

### `jacobian`

```text
jacobian(value)
```

Return the defining-equation Jacobian of an ideal, ring or embedded scheme.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L213)

### `tangent_space`

```text
tangent_space(value, point)
```

Return a column basis of tangent directions after extending the point field.

Projective inputs use the affine chart of the first nonzero point coordinate. Nonreduced and mixed-dimensional schemes retain their defining equations.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L225)

### `tangent_ideal`

```text
tangent_ideal(value, point)
```

Return the embedded affine tangent space, translated through the chosen point.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L242)

### `singular_locus`

```text
singular_locus(value, *, reduced=False, assume_equidimensional=False)
```

Compute a singular subscheme, inferring equidimensionality when necessary.

Equidimensional inputs retain the Jacobian scheme structure unless reduced is requested. Mixed-dimensional inputs require reduced=True: component dimension strata are saturated away from higher-dimensional components, then their singular sets are joined. Scheme inputs return schemes; projective results are saturated at the irrelevant ideal.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/maps.py#L265)

## `alggeom.modules`

Finite module presentations, polynomial syzygies, and Kahler differentials.

### `Submodule`

```text
class Submodule(ring, generators)
```

The column span of a matrix in a finite free module over S or S/I.

Columns are retained in their supplied order, including zero columns, so membership coefficients and syzygy coordinates refer to the input matrix.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L85)

Stored fields: `ring`, `generators`.

#### `Submodule.ambient_rank`

The number of rows of the generator matrix, not the module's rank.

#### `Submodule.normal_form`

```text
normal_form(self, vector)
```

Reduce a column modulo the submodule using its module Groebner basis.

#### `Submodule.contains`

```text
contains(self, vector)
```

Decide whether a column belongs to this submodule.

#### `Submodule.coordinates`

```text
coordinates(self, vector)
```

Return coefficients in the input columns, or None for a nonmember.

Over S/I the resulting identity holds modulo I in each row. Auxiliary coefficients for the lifted ideal relations are omitted.

#### `Submodule.syzygies`

```text
syzygies(self)
```

Return a matrix whose columns generate the relations among the input columns.

### `PresentedModule`

```text
class PresentedModule(ring, presentation)
```

The cokernel of a matrix A^s -> A^r over a polynomial or coordinate ring.

Rows index the r module generators and columns the s defining relations. Presentations are finite; they are not assumed minimal or free resolutions.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L192)

Stored fields: `ring`, `presentation`.

#### `PresentedModule.element`

```text
element(self, vector)
```

Construct a residue class in this presented module.

#### `PresentedModule.gens`

The classes of the ordered presentation generators.

#### `PresentedModule.zero`

The additive identity in this module.

#### `PresentedModule.hom`

```text
hom(self, matrix, target)
```

Define a linear map by columns of generator images.

#### `PresentedModule.tensor_product`

```text
tensor_product(self, other)
```

Return the tensor product over the current base ring.

#### `PresentedModule.hom_module`

```text
hom_module(self, target)
```

Return the module of linear maps to target.

#### `PresentedModule.annihilator`

```text
annihilator(self)
```

Return the annihilator in the ambient polynomial ring.

#### `PresentedModule.generator_count`

The number of generators in this presentation.

#### `PresentedModule.relation_count`

The number of relation columns in this presentation.

#### `PresentedModule.relations`

The submodule of the free module generated by the presentation columns.

#### `PresentedModule.normal_form`

```text
normal_form(self, vector)
```

Return a representative of the class of a column in this cokernel.

#### `PresentedModule.is_zero`

Decide whether the presented module is zero.

#### `PresentedModule.syzygies`

```text
syzygies(self)
```

Return generators for the relations among the presentation's relation columns.

#### `PresentedModule.fitting_ideal`

```text
fitting_ideal(self, index)
```

Return the indexed Fitting ideal, as an ideal in the ambient polynomial ring.

#### `PresentedModule.restrict_scalars`

```text
restrict_scalars(self)
```

Present this module over its ambient polynomial ring S.

For a module over S/I, append I times every free generator to its relation matrix. This changes the base of a subsequent resolution: an S-resolution is finite, whereas an S/I-resolution may be infinite.

#### `PresentedModule.free_resolution`

```text
free_resolution(self, *, generator_degrees=None, max_length=None)
```

Compute a minimal standard-graded resolution over a polynomial ring.

### `ModuleElement`

```text
class ModuleElement(module, vector)
```

A parent-aware residue class in a finitely presented module.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L369)

Stored fields: `module`, `vector`.

### `ModuleHomomorphism`

```text
class ModuleHomomorphism(source, target, matrix)
```

A linear map of presented modules, with columns the images of source generators.

Relations are checked in the target cokernel. Kernel and image objects retain presentations; kernel_map and image_map supply their natural inclusions. Cokernel_map is the natural target quotient.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L424)

Stored fields: `source`, `target`, `matrix`.

#### `ModuleHomomorphism.compose`

```text
compose(self, inner)
```

Return self after a map into its source module.

#### `ModuleHomomorphism.kernel`

The kernel as a finitely presented module over the same base.

#### `ModuleHomomorphism.kernel_map`

The natural inclusion of kernel into source.

#### `ModuleHomomorphism.image`

The image, presented by relations among its supplied generator images.

#### `ModuleHomomorphism.image_map`

The natural inclusion of the presented image into target.

#### `ModuleHomomorphism.cokernel`

The target divided by this map's image.

#### `ModuleHomomorphism.cokernel_map`

The natural quotient from target to cokernel.

#### `ModuleHomomorphism.is_injective`

Whether the kernel is zero.

#### `ModuleHomomorphism.is_surjective`

Whether the cokernel is zero.

#### `ModuleHomomorphism.lift`

```text
lift(self, element)
```

Return a source preimage of a target element, raising for a nonmember.

#### `ModuleHomomorphism.inverse`

```text
inverse(self)
```

Return the inverse of an isomorphism; nonisomorphisms raise.

### `HomModule`

```text
class HomModule(source, target)
```

The module Hom(M,N), with conversions between elements and actual linear maps.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L567)

Stored fields: `ring`, `presentation`, `source`, `target`.

#### `HomModule.element`

```text
element(self, vector)
```

Construct a residue class in this presented module.

#### `HomModule.gens`

The classes of the ordered presentation generators.

#### `HomModule.zero`

The additive identity in this module.

#### `HomModule.hom`

```text
hom(self, matrix, target)
```

Define a linear map by columns of generator images.

#### `HomModule.tensor_product`

```text
tensor_product(self, other)
```

Return the tensor product over the current base ring.

#### `HomModule.hom_module`

```text
hom_module(self, target)
```

Return the module of linear maps to target.

#### `HomModule.annihilator`

```text
annihilator(self)
```

Return the annihilator in the ambient polynomial ring.

#### `HomModule.generator_count`

The number of generators in this presentation.

#### `HomModule.relation_count`

The number of relation columns in this presentation.

#### `HomModule.relations`

The submodule of the free module generated by the presentation columns.

#### `HomModule.normal_form`

```text
normal_form(self, vector)
```

Return a representative of the class of a column in this cokernel.

#### `HomModule.is_zero`

Decide whether the presented module is zero.

#### `HomModule.syzygies`

```text
syzygies(self)
```

Return generators for the relations among the presentation's relation columns.

#### `HomModule.fitting_ideal`

```text
fitting_ideal(self, index)
```

Return the indexed Fitting ideal, as an ideal in the ambient polynomial ring.

#### `HomModule.restrict_scalars`

```text
restrict_scalars(self)
```

Present this module over its ambient polynomial ring S.

For a module over S/I, append I times every free generator to its relation matrix. This changes the base of a subsequent resolution: an S-resolution is finite, whereas an S/I-resolution may be infinite.

#### `HomModule.free_resolution`

```text
free_resolution(self, *, generator_degrees=None, max_length=None)
```

Compute a minimal standard-graded resolution over a polynomial ring.

#### `HomModule.homomorphism`

```text
homomorphism(self, element)
```

Interpret a Hom-module element as a map M -> N.

#### `HomModule.from_homomorphism`

```text
from_homomorphism(self, mapping)
```

Express a map in this Hom presentation, modulo target relations.

### `syzygies`

```text
syzygies(matrix, ring)
```

Return column generators of the exact kernel of a polynomial matrix.

The matrix acts on column vectors. For a coordinate ring S/I, relations are computed modulo I by lifting to a free S-module with the I-relations included. The output need not be a minimal generating set or a basis.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L160)

### `fitting_ideal`

```text
fitting_ideal(module, index)
```

Return Fitt_index(M), generated by minors of size generator_count-index.

For A=S/I the result is the inverse-image ideal in S and contains I. Minors of size <=0 generate the unit ideal; a size larger than either matrix dimension gives the zero ideal of A, whose inverse image is I.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L305)

### `differentials`

```text
differentials(ring)
```

Present Omega_(A/K) by the transposed Jacobian of A=S/I.

The generators are the ordered ambient coordinate differentials. Both defining equations and their differential relations are respected over A.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L339)

### `differential`

```text
differential(expression, ring)
```

Return the coordinate column representing df in Omega_(A/K).

Coefficient-field parameters are constants for this relative differential. Equal representatives in A give equal classes in the differential module.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L354)

### `direct_sum`

```text
direct_sum(*modules)
```

Return the direct sum of one or more modules over the same ring.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L526)

### `tensor_product`

```text
tensor_product(left, right)
```

Present a tensor product with basis order e_i tensor f_j, j varying fastest.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L544)

### `annihilator`

```text
annihilator(module)
```

Return the annihilator, lifted to the ambient polynomial ring for a quotient base.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/modules.py#L553)

## `alggeom.morphism_invariants`

Degree, finiteness and critical loci of exact affine polynomial maps.

### `is_dominant`

```text
is_dominant(mapping)
```

Test density of an affine morphism's image in its reduced target.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/morphism_invariants.py#L9)

### `is_finite`

```text
is_finite(mapping)
```

Test module finiteness of an affine morphism by its graph initial ideal.

With source variables preceding target variables in lexicographic order, the graph coordinate ring is finite over the target coordinate ring exactly when its initial ideal contains a pure power of every source variable. The criterion also applies to nonreduced targets and positive characteristic.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/morphism_invariants.py#L35)

### `morphism_degree`

```text
morphism_degree(mapping, *, max_dimension=4096)
```

Return the generic finite-fiber length of a dominant affine morphism.

The target must be affine space over a characteristic-zero exact field. Generic length agrees with geometric fiber cardinality when the generic fiber is reduced, but counts multiplicities for nonreduced fibers. The morphism need not be finite: generically finite open maps have degree.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/morphism_invariants.py#L88)

### `is_birational`

```text
is_birational(mapping)
```

Test whether an integral dominant morphism induces degree-one fields.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/morphism_invariants.py#L99)

### `ramification_locus`

```text
ramification_locus(mapping)
```

Return the non-unramified locus of a finite affine morphism.

The closed subscheme is defined by Fitt_0(Omega_{B/A}), where B is the source coordinate algebra and A is the target coordinate algebra. Its presentation has one row per source coordinate and columns given by the Jacobian of the source equations and the pullbacks of target coordinates. This also applies when the source or target is singular or nonreduced. It detects failure of unramifiedness, not the more general failure of smoothness or flatness.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/morphism_invariants.py#L108)

### `branch_locus`

```text
branch_locus(mapping)
```

Return the closed image of the non-unramified locus of a finite affine map.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/morphism_invariants.py#L136)

## `alggeom.noether`

Constructive Noether normalization with finite-module and integrality data.

### `NoetherNormalization`

```text
class NoetherNormalization(coordinate_map, inverse_map, base_embedding, integral_relations, module_generators)
```

A finite polynomial subring and an explicit isomorphic coordinate presentation.

coordinate_map goes from the input ring to the transformed ring. The polynomial base embeds through base_embedding. Integral relations are triangular and monic in the nonbase coordinates. module_generators span over the polynomial base; they are not asserted to be a free basis.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/noether.py#L75)

Stored fields: `coordinate_map`, `inverse_map`, `base_embedding`, `integral_relations`, `module_generators`.

#### `NoetherNormalization.ring`

The transformed coordinate ring, isomorphic to the input.

#### `NoetherNormalization.base`

The embedded polynomial coordinate ring.

#### `NoetherNormalization.projection`

The finite affine-scheme morphism to the polynomial base.

#### `NoetherNormalization.generic_rank`

The finite vector-space length after extending the polynomial base to its field.

#### `NoetherNormalization.verify`

```text
verify(self)
```

Check inverse coordinates, base injection and the monic finite-module certificate.

### `noether_normalization`

```text
noether_normalization(value, *, names=None, max_linear_attempts=64, max_transform_degree=256, max_generators=1024)
```

Construct a certified finite projection to affine space over the coefficient field.

Successive coordinate automorphisms make an elimination equation monic. The fallback uses triangular polynomial changes when the field is too small for a suitable linear direction. Resource exhaustion raises without claiming that an uncertified generically finite projection is finite.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/noether.py#L165)

## `alggeom.nonlinear_strata`

Exact parameter strata for nonlinear polynomial hypersurface fibers.

A univariate fiber is stratified by its specialized polynomial degree and Sylvester-matrix rank, which determines the number of distinct geometric roots. Multivariate hypersurfaces are stratified by coefficient vanishing: a nonzero nonconstant hypersurface has dimension n-1 over an algebraic closure.

### `NonlinearFiberPiece`

```text
class NonlinearFiberPiece(locus, invariants)
```

A constructible parameter locus with fixed geometric fiber invariants.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/nonlinear_strata.py#L20)

Stored fields: `locus`, `invariants`.

#### `NonlinearFiberPiece.dimension`

Geometric dimension of a fiber on this parameter locus.

#### `NonlinearFiberPiece.geometric_points`

Number of distinct geometric points, or infinity for positive dimension.

### `NonlinearFiberStratification`

```text
class NonlinearFiberStratification(equation, variables, parameters, pieces)
```

A complete coefficient-rank partition for one polynomial equation over QQ.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/nonlinear_strata.py#L38)

Stored fields: `equation`, `variables`, `parameters`, `pieces`.

#### `NonlinearFiberStratification.specialize`

```text
specialize(self, values)
```

Evaluate exact geometric dimension and distinct-point cardinality.

### `hypersurface_fiber_stratification`

```text
hypersurface_fiber_stratification(equation, variables, parameters, *, max_minors=4096)
```

Stratify arbitrary polynomial hypersurface fibers over QQ exactly.

Fiber equations may be nonlinear in any number of variables and contain any number of polynomial parameters. In one fiber variable, Sylvester ranks yield exact counts of distinct geometric points, including multiple roots and degree drops. With at least two fiber coordinates, the hypersurface dimension is exact; positive-dimensional fibers have infinitely many geometric points. The determinant budget is explicit.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/nonlinear_strata.py#L86)

### `SeparableFiberStratification`

```text
class SeparableFiberStratification(equations, variables, parameters, factors, pieces)
```

Exact product stratification of uncoupled nonlinear fiber equations.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/nonlinear_strata.py#L183)

Stored fields: `equations`, `variables`, `parameters`, `factors`, `pieces`.

#### `SeparableFiberStratification.specialize`

```text
specialize(self, values)
```

Combine exact invariants of the independent univariate fibers.

### `separable_fiber_stratification`

```text
separable_fiber_stratification(equations, variables, parameters, *, max_pieces=4096)
```

Stratify products of independent nonlinear univariate equations.

Each fiber coordinate may appear in at most one defining equation, and each equation may depend on at most one fiber coordinate. Coefficients can be arbitrary polynomials in all declared parameters. This includes nonreduced and singular zero-dimensional fibers. Coupled systems require a comprehensive multivariate Gröbner algorithm and are not accepted.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/nonlinear_strata.py#L203)

## `alggeom.normalization`

Affine integral closure by radical test ideals and finite endomorphism rings.

### `NormalizationStep`

```text
class NormalizationStep(source, target, test_ideal, denominator, numerators)
```

One finite birational extension End_A(J)=(gJ:J)/g, with replay data.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/normalization.py#L62)

Stored fields: `source`, `target`, `test_ideal`, `denominator`, `numerators`.

#### `NormalizationStep.verify`

```text
verify(self)
```

Replay the test ideal, colon module and saturated graph presentation.

### `IntegralClosure`

```text
class IntegralClosure(mapping, steps=(), inverse_functions=(), module_generators=(), components=())
```

A normalization ring and its finite birational coordinate-ring map.

For a domain, inverse_functions express all normalized coordinates in the original function field, and module_generators span the integral closure over the input ring. For a reduced ring with several components, components contains their domain normalizations and ring is their explicit product.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/normalization.py#L106)

Stored fields: `mapping`, `steps`, `inverse_functions`, `module_generators`, `components`.

#### `IntegralClosure.source`

The original reduced coordinate ring.

#### `IntegralClosure.ring`

The integral-closure coordinate ring.

#### `IntegralClosure.conductor`

The conductor ideal in the original ambient ring, for a domain input.

#### `IntegralClosure.verify`

```text
verify(self)
```

Check finite-extension steps, birational coordinates and final normality.

### `is_normal`

```text
is_normal(value, **limits)
```

Decide normality over characteristic-zero fields, including reduced disconnected rings.

Distinct irreducible components of a normal Noetherian ring are disjoint. Each component is then tested by equality A=End_A(J) for a radical Jacobian test ideal J. The zero ring is normal by the empty-scheme convention.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/normalization.py#L36)

### `integral_closure`

```text
integral_closure(value, *, max_steps=16, max_new_variables=64, max_module_generators=1024, max_components=32, **decomposition_limits)
```

Compute the affine normalization with a finite birational map and exact evidence.

Inputs must be reduced algebras over supported characteristic-zero fields. Domains use the Grauert-Remmert endomorphism-ring algorithm. Several components are normalized separately and assembled as a product. Limits raise without returning an intermediate ring as a completed normalization.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/normalization.py#L279)

### `normalize`

```text
normalize(value, **limits)
```

Return the integral-closure coordinate ring of a reduced affine algebra.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/normalization.py#L326)

### `normalization_map`

```text
normalization_map(value, **limits)
```

Return the normalization morphism for an affine scheme, or its ring homomorphism.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/normalization.py#L331)

## `alggeom.parameter_cover`

Recursive Gröbner covers over prime parameter ideals.

A branch works in Frac(QQ[a]/P) for a certified prime ideal P. Exceptional leading coefficients enlarge P; their minimal primes are recursively visited. The covers describe algebraic-closure parameter points, but the coefficient fraction field and associated primes are computed over QQ.

### `ParameterCoverLimitError`

```text
class ParameterCoverLimitError(RuntimeError)
```

A recursive algebraic parameter cover exceeds the declared budget.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_cover.py#L21)

### `PrimeCoverPiece`

```text
class PrimeCoverPiece(condition, certificate)
```

Constructible locus and a generic-basis proof valid everywhere on it.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_cover.py#L25)

Stored fields: `condition`, `certificate`.

#### `PrimeCoverPiece.basis`

Return the certified generic lex Gröbner basis on this locus.

#### `PrimeCoverPiece.verify`

```text
verify(self)
```

Check exact polynomial identities and the declared open guards.

### `PrimeParameterCover`

```text
class PrimeParameterCover(pieces, parameters, variables, equations)
```

A partition into parameter loci with replayable generic bases.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_cover.py#L49)

Stored fields: `pieces`, `parameters`, `variables`, `equations`.

#### `PrimeParameterCover.verify`

```text
verify(self)
```

Replay all proofs and verify symbolic coverage/disjointness.

#### `PrimeParameterCover.specialize`

```text
specialize(self, values)
```

Evaluate the unique stratum basis at exact algebraic parameters.

### `comprehensive_prime_cover`

```text
comprehensive_prime_cover(equations, variables, parameters, *, max_nodes=40, max_pieces=80, max_pairs=150)
```

Recursively certify a Gröbner cover using prime quotient function fields.

For each branch V(P) minus specified closed loci, compute a generic basis in Fr(QQ[a]/P). Every new inverted coefficient defines an exceptional variety. Recursion proceeds through its field-relative minimal primes; overlaps are partitioned by exact polynomial equations/inequations. Failure to certify an exceptional ideal raises instead of returning a partial cover. The result is bounded because recursive Gröbner covers can have exponential complexity.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_cover.py#L104)

## `alggeom.parameter_fields`

Exact fraction fields of prime polynomial parameter quotients.

A prime ideal P in QQ[a] defines the integral domain A=QQ[a]/P and its fraction field Frac(A). Field elements are fractions of polynomial normal forms, with equality tested by cross multiplication modulo P. This makes rational arithmetic possible without solving P for a particular parameter.

### `ParameterFieldError`

```text
class ParameterFieldError(ValueError)
```

The supplied parameter ideal has no certified fraction field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_fields.py#L16)

### `PrimeParameterField`

```text
class PrimeParameterField(prime_ideal)
```

Fraction field of `QQ[parameters]/prime_ideal`.

The constructor proves that the ideal is prime before allowing division; a primary but nonprime ideal is rejected. Computations use exact Gröbner normal forms and work on implicit nonlinear parameter varieties without selecting a radical or rational parametrization.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_fields.py#L49)

#### `PrimeParameterField.reduce`

```text
reduce(self, expression)
```

Return the canonical polynomial normal form modulo the prime ideal.

#### `PrimeParameterField.element`

```text
element(self, numerator, denominator=1)
```

Construct a quotient-field element, rejecting zero divisors.

#### `PrimeParameterField.coerce`

```text
coerce(self, value)
```

Convert a rational expression or compatible field element.

#### `PrimeParameterField.records`

```text
records(self)
```

List polynomial conditions required by prior field inversions.

### `ParameterFraction`

```text
class ParameterFraction(field, numerator, denominator)
```

An exact fraction in a certified prime parameter quotient field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_fields.py#L98)

Stored fields: `field`, `numerator`, `denominator`.

#### `ParameterFraction.inverse`

```text
inverse(self)
```

Invert a nonzero residue class and record its specialization guard.

#### `ParameterFraction.as_expr`

```text
as_expr(self)
```

Return a representative rational expression in the parameters.

### `QuotientGroebnerResult`

```text
class QuotientGroebnerResult(basis, denominator_guards, parameter_ideal, variables, equations, membership_witnesses)
```

A generic lex basis with independently replayable membership witnesses.

Identities hold in the fraction field of the declared prime quotient. Every coefficient inverse used by the construction is recorded as a nonvanishing guard. Guard conditions must be imposed on any parameter stratum before specializing the generic basis.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_fields.py#L254)

Stored fields: `basis`, `denominator_guards`, `parameter_ideal`, `variables`, `equations`, `membership_witnesses`.

#### `QuotientGroebnerResult.verify`

```text
verify(self)
```

Verify source membership, Buchberger's criterion and all guards.

The verifier never calls the basis constructor. Polynomial identities are checked modulo the parameter prime using exact fraction-field arithmetic; the quotient field must itself be certified prime.

### `groebner_over_prime_parameters`

```text
groebner_over_prime_parameters(equations, variables, parameter_ideal, *, max_pairs=200)
```

Compute a lex Gröbner basis over `Frac(QQ[a]/P)`.

`P` must be a certified prime ideal. Every inverted residue numerator and coefficient denominator is returned as a nonvanishing guard; the result is a generic fiber basis on V(P) away from those guards. It makes no claim about fibers on the vanishing locus of a guard.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_fields.py#L335)

## `alggeom.parameter_localization`

Exact arithmetic on principal opens of arbitrary parameter quotient rings.

For A = QQ[a]/I and h = product(inverted), the localization A[h^-1] has kernel (I:h^infinity)/I. Consequently equality is decidable by reducing cross products modulo the *saturated* ideal. No primality, reducedness, or non-zero-divisor hypothesis is imposed: components on which h vanishes identically disappear upon localization, as algebraic geometry requires.

This is a general, finitely-generated localization, NOT the total fraction ring of a non-domain, and it does not perform a comprehensive cover search.

### `LocalizedParameterElement`

```text
class LocalizedParameterElement(localization, numerator, powers)
```

A fraction represented by a numerator and explicit inversion powers.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_localization.py#L147)

Stored fields: `localization`, `numerator`, `powers`.

#### `LocalizedParameterElement.is_zero`

```text
is_zero(self)
```

Decide whether the fraction is zero on this localized quotient.

#### `LocalizedParameterElement.as_expr`

```text
as_expr(self)
```

Return a rational-expression representative; equality uses saturation.

### `ParameterLocalization`

```text
class ParameterLocalization(ideal, inverted=())
```

QQ[a]/I localized at finitely many explicitly inverted polynomials.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/parameter_localization.py#L21)

Stored fields: `ideal`, `inverted`.

#### `ParameterLocalization.kernel_ideal`

Preimage of zero under A -> S^-1 A, computed by exact saturation.

#### `ParameterLocalization.reduce`

```text
reduce(self, polynomial)
```

Compute the canonical class modulo the localization kernel.

#### `ParameterLocalization.element`

```text
element(self, numerator, powers=None)
```

Represent numerator / product(inverted[i]**powers[i]).

#### `ParameterLocalization.coerce`

```text
coerce(self, value)
```

Coerce a polynomial or an element of this exact localization.

#### `ParameterLocalization.inverse_generator`

```text
inverse_generator(self, index)
```

Return the certified inverse of one declared localization generator.

#### `ParameterLocalization.localize_at`

```text
localize_at(self, *polynomials)
```

Refine to the intersection with additional principal opens.

#### `ParameterLocalization.condition`

Constructible condition of this quotient's principal open.

#### `ParameterLocalization.split_at`

```text
split_at(self, polynomial)
```

Partition this locus into D(polynomial) and V(polynomial).

Each branch is a genuine localization or quotient, even when the parameter ideal has embedded components or zero divisors. An empty branch is represented by None; coverage/disjointness follows from f != 0 or f == 0 on every geometric parameter point.

#### `ParameterLocalization.split_ideal`

```text
split_ideal(self, exceptional)
```

Disjointly split D(J) and V(J) for any finitely generated ideal J.

D(J) = D(f_1) union [V(f_1) intersect D(f_2)] union ...; the pieces are disjoint by construction. The last branch is V(J), without assuming J is radical, prime or equidimensional. This provides the exact branching infrastructure, not a Gröbner cover on zero-divisor quotient rings by itself.

#### `ParameterLocalization.restrict`

```text
restrict(self, element)
```

Transport a fraction from an earlier prefix localization to this one.

## `alggeom.points`

Complete finite-field enumeration and bounded-height rational point searches.

### `point_search`

```text
point_search(scheme, height=None, *, max_candidates=1000000)
```

Return all finite-field points, or all QQ points in a specified height box.

Affine rational height is max(|numerator|,denominator) in each reduced coordinate. Projective height is the maximum absolute entry of a primitive integer representative. A bounded empty QQ search is not a proof that no rational point exists. Exhausting the candidate budget raises without returning a prefix as a complete bounded enumeration.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/points.py#L34)

## `alggeom.projective_normalization`

Explicit localization transitions between normalized projective-curve charts.

### `NormalizationPatch`

```text
class NormalizationPatch(chart, mapping, closure)
```

A normalized affine chart and its finite birational morphism.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/projective_normalization.py#L10)

Stored fields: `chart`, `mapping`, `closure`.

#### `NormalizationPatch.normalized_scheme`

The normal affine source of the finite chart map.

### `NormalizationTransition`

```text
class NormalizationTransition(source_index, target_index, guard, coordinates)
```

A regular transition on D(guard), with coordinates in its localization.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/projective_normalization.py#L24)

Stored fields: `source_index`, `target_index`, `guard`, `coordinates`.

### `ProjectiveNormalizationAtlas`

```text
class ProjectiveNormalizationAtlas(original, patches)
```

Finite normalized charts, with explicit regular transition coordinates.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/projective_normalization.py#L34)

Stored fields: `original`, `patches`.

#### `ProjectiveNormalizationAtlas.overlap_guard`

```text
overlap_guard(self, source_index, target_index)
```

Return the pullback of the target's distinguished homogeneous coordinate.

#### `ProjectiveNormalizationAtlas.transition`

```text
transition(self, source_index, target_index, *, max_localization_power=32)
```

Construct normalized chart transitions explicitly in A[guard^-1].

The normalization's inverse rational functions express generators of the target patch's integral closure in its original curve field. Transporting these fractions through projective chart coordinates produces functions regular on the overlap. Ideal membership gives an explicit numerator divided by a power of the overlap guard; failure to establish regularity raises instead of returning a rational map.

#### `ProjectiveNormalizationAtlas.verify_cocycle`

```text
verify_cocycle(self, *, max_localization_power=32)
```

Verify identity and composition of transition maps in chart function fields.

Algebraic equality on each nonempty triple overlap is checked in the normalized coordinate ring's exact fraction field. Regularity of the individual transitions was checked on their principal-open overlaps.

### `projective_normalization_atlas`

```text
projective_normalization_atlas(curve, **limits)
```

Normalize standard projective curve charts and retain integral coordinates.

Projective normalization is represented as an affine atlas with finite birational maps; each transition is computed and verified as a regular function on the localized overlap. No global projective embedding is constructed by this representation.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/projective_normalization.py#L161)

## `alggeom.properties`

Geometric predicates with explicit field-relative and empty-scheme conventions.

### `is_reduced`

```text
is_reduced(value, **limits)
```

Decide whether the scheme or coordinate ring has no nonzero nilpotents.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L23)

### `is_irreducible`

```text
is_irreducible(value, **limits)
```

Decide irreducibility over the declared field; the empty scheme is not irreducible.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L29)

### `is_integral`

```text
is_integral(value, **limits)
```

Decide whether the coordinate ring is a nonzero domain over its declared field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L34)

### `is_equidimensional`

```text
is_equidimensional(value, **limits)
```

Decide equality of minimal-component dimensions; the empty scheme satisfies it.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L41)

### `is_unmixed`

```text
is_unmixed(value, **limits)
```

Decide equality of all associated-prime dimensions, including embedded primes.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L46)

### `codimension`

```text
codimension(value)
```

Return embedded codimension, with infinity for the empty scheme.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L51)

### `is_smooth`

```text
is_smooth(value, point=None)
```

Decide smoothness globally, or at a supplied exact field-valued point.

Pointwise computations use the local tangent cone and embedding dimension. Global mixed-dimensional calculations use the reduced singular locus.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/properties.py#L57)

## `alggeom.rational_functions`

Parent-aware rational functions, principal localizations and local rings.

### `FractionField`

```text
class FractionField(base)
```

The field of rational functions of a nonzero integral coordinate ring.

Fractions are compared by cross multiplication modulo the defining ideal. Polynomial gcd cancellation is only a representation simplification; a unique numerator/denominator form is not assumed in a nonfactorial ring.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rational_functions.py#L123)

Stored fields: `base`.

#### `FractionField.fraction`

```text
fraction(self, numerator, denominator=1)
```

Construct a fraction, checking its base ring and allowed denominators.

Coordinate elements must belong to the base coordinate ring. Use as_expr() explicitly to reinterpret a representative in another ring.

#### `FractionField.gens`

The images of the base coordinate generators.

#### `FractionField.zero`

The additive identity.

#### `FractionField.one`

The multiplicative identity.

#### `FractionField.domain`

The coefficient field.

#### `FractionField.kernel`

The zero ideal of the base, lifted to its polynomial ambient ring.

### `Localization`

```text
class Localization(base, denominators, max_power=256)
```

The ring obtained by inverting finitely many functions in a coordinate ring.

An explicit affine model adjoins u with u*product(denominators)=1. Kernel saturation preserves correct equality when localization removes components.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rational_functions.py#L148)

Stored fields: `base`, `denominators`, `max_power`.

#### `Localization.fraction`

```text
fraction(self, numerator, denominator=1)
```

Construct a fraction, checking its base ring and allowed denominators.

Coordinate elements must belong to the base coordinate ring. Use as_expr() explicitly to reinterpret a representative in another ring.

#### `Localization.gens`

The images of the base coordinate generators.

#### `Localization.zero`

The additive identity.

#### `Localization.one`

The multiplicative identity.

#### `Localization.domain`

The coefficient field.

#### `Localization.inverted`

The product of the specified invertible functions.

#### `Localization.kernel`

The contraction of zero from the localization: I:f^infinity.

#### `Localization.affine_model`

An affine coordinate-ring presentation of the localization.

#### `Localization.natural_map`

The canonical map from the base ring to the affine localization model.

### `LocalRing`

```text
class LocalRing(base, prime)
```

The localization of a coordinate ring at a specified prime ideal.

The prime is represented in the ambient polynomial ring and must contain the defining ideal. Component contraction computes equality with nilpotents and zero divisors. General local standard bases and completion are not provided by this arithmetic interface.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rational_functions.py#L213)

Stored fields: `base`, `prime`.

#### `LocalRing.fraction`

```text
fraction(self, numerator, denominator=1)
```

Construct a fraction, checking its base ring and allowed denominators.

Coordinate elements must belong to the base coordinate ring. Use as_expr() explicitly to reinterpret a representative in another ring.

#### `LocalRing.gens`

The images of the base coordinate generators.

#### `LocalRing.zero`

The additive identity.

#### `LocalRing.one`

The multiplicative identity.

#### `LocalRing.domain`

The coefficient field.

#### `LocalRing.kernel`

The contraction of zero after localizing at the prime.

#### `LocalRing.maximal_ideal`

The ambient inverse image of the maximal ideal of this local ring.

#### `LocalRing.residue_field`

The fraction field of the integral residue coordinate ring.

### `RationalFunction`

```text
class RationalFunction(parent, numerator, denominator=1)
```

A fraction in an explicit fraction field, localization or local ring.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rational_functions.py#L263)

Stored fields: `parent`, `numerator`, `denominator`.

#### `RationalFunction.as_expr`

```text
as_expr(self)
```

Return an ambient rational expression, removing the quotient-ring parent.

#### `RationalFunction.inverse`

```text
inverse(self)
```

Invert a unit; a nonunit or zero raises ZeroDivisionError.

#### `RationalFunction.is_unit`

Whether the inverse belongs to the same parent.

#### `RationalFunction.to_coordinate`

```text
to_coordinate(self)
```

Return a base-ring element if this rational function is globally regular.

### `fraction_field`

```text
fraction_field(value)
```

Return the rational-function field of an integral coordinate ring.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rational_functions.py#L368)

### `local_ring`

```text
local_ring(value, point_or_prime)
```

Return a local ring at a prime or at an exact affine/projective point.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rational_functions.py#L373)

## `alggeom.resolutions`

Minimal standard-graded polynomial-module resolutions and Betti invariants.

### `GradedResolution`

```text
class GradedResolution(module, generator_degrees, differentials, shifts, augmentation, complete)
```

A minimal graded complex F_i -> F_(i-1) augmenting to a presented module.

`shifts[i]` lists j for the summands S(-j) in F_i. `differentials[i-1]` maps F_i to F_(i-1). The augmentation maps F_0 to the original presentation generators. A truncated complex has `complete=False` and does not supply projective dimension, regularity or a module Hilbert series.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/resolutions.py#L89)

Stored fields: `module`, `generator_degrees`, `differentials`, `shifts`, `augmentation`, `complete`.

#### `GradedResolution.betti_numbers`

The dictionary beta_(i,j), indexed by homological and internal degrees.

#### `GradedResolution.betti_table`

Rows j-i mapped to tuples of beta_(i,j), with homological degree in columns.

A truncated complex supplies the known Betti columns through its last term; higher homological degrees have not been computed.

#### `GradedResolution.projective_dimension`

The module's projective dimension; the zero module has value -1.

#### `GradedResolution.regularity`

Castelnuovo-Mumford regularity; the zero module has value -oo.

#### `GradedResolution.hilbert_series`

```text
hilbert_series(self, *, parameter=None)
```

Return the alternating graded free-module Hilbert series of the resolved module.

#### `GradedResolution.verify`

```text
verify(self)
```

Check grading, minimality, augmentation and exactness by independent kernel tests.

For a truncated complex, exactness is checked at every term except its leftmost term. Completeness additionally requires a zero terminal kernel. Verification recomputes syzygies and can cost as much as construction.

### `free_resolution`

```text
free_resolution(module, *, generator_degrees=None, max_length=None)
```

Compute a minimal standard-graded free resolution over S=K[x_1,...,x_n].

An Ideal means the quotient S/I. A PresentedModule must be over S; use `restrict_scalars()` explicitly for a module over S/I. Integer generator shifts default to zero. Relations must be homogeneous for these shifts. Successive exact syzygies are minimized by graded membership. The default length bound is n, the Hilbert syzygy bound. A smaller explicit bound returns a truncated complex marked incomplete when another term is needed.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/resolutions.py#L211)

## `alggeom.rings`

Elements and coefficient-preserving homomorphisms of affine coordinate rings.

### `CoordinateElement`

```text
class CoordinateElement(ring, expression)
```

An element of S/I with exact arithmetic and equality modulo I.

Arithmetic requires the same parent coordinate ring, or an expression coercible into it. as_expr() removes the parent and returns a polynomial representative. Elements are unhashable because equal quotient rings may use different term orders and hence different representatives.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rings.py#L43)

Stored fields: `ring`, `expression`.

#### `CoordinateElement.as_expr`

```text
as_expr(self)
```

Return the reduced expression without its coordinate-ring parent.

#### `CoordinateElement.is_unit`

Whether an inverse exists in this coordinate ring, without assuming finiteness.

#### `CoordinateElement.inverse`

```text
inverse(self)
```

Return the inverse in S/I, raising ZeroDivisionError for a nonunit.

### `RingHomomorphism`

```text
class RingHomomorphism(source, target, images)
```

A homomorphism S/I -> T/J defined by images of the S coordinates.

Source and target have the same declared coefficient field. Construction checks all source relations. Ideal outputs are inverse-image ideals in the relevant ambient polynomial ring, including that ring's defining ideal.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/rings.py#L158)

Stored fields: `source`, `target`, `images`.

#### `RingHomomorphism.graph_ideal`

The graph in target-first lex order, with fresh source coordinate symbols.

#### `RingHomomorphism.kernel`

The kernel as an ambient ideal containing the source defining ideal.

#### `RingHomomorphism.preimage`

```text
preimage(self, ideal)
```

Contract a target ideal along this homomorphism.

#### `RingHomomorphism.extend`

```text
extend(self, ideal)
```

Generate the target ideal from images of a source ideal's generators.

#### `RingHomomorphism.is_injective`

Whether the map is injective, using exact kernel contraction.

#### `RingHomomorphism.lift`

```text
lift(self, element)
```

Return one polynomial preimage of a target element, or raise ValueError.

#### `RingHomomorphism.is_surjective`

Whether all target coordinate classes have polynomial preimages.

#### `RingHomomorphism.inverse`

```text
inverse(self)
```

Return the inverse of a bijective homomorphism.

#### `RingHomomorphism.compose`

```text
compose(self, inner)
```

Return self after inner, requiring the intermediate coordinate rings to agree.

#### `RingHomomorphism.base_change`

```text
base_change(self, domain)
```

Extend source, target and coordinate images to the same exact coefficient field.

## `alggeom.schemes`

Affine and projective schemes, exact points, charts and polynomial morphisms.

### `AffineScheme`

```text
class AffineScheme(ideal)
```

The affine scheme Spec(S/I), retaining nilpotents and the field of definition.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L54)

Stored fields: `ideal`.

#### `AffineScheme.ambient`

The ambient polynomial ring.

#### `AffineScheme.variables`

The ordered ambient coordinate symbols.

#### `AffineScheme.domain`

The exact coefficient field.

#### `AffineScheme.coordinate_ring`

The coordinate ring S/I, including its Krull and vector-space dimensions.

#### `AffineScheme.dimension`

The Krull dimension; the empty scheme has dimension -1.

#### `AffineScheme.is_empty`

Scheme emptiness, which is different from having no rational or real points.

#### `AffineScheme.is_reduced`

The field-relative is reduced property of this scheme.

#### `AffineScheme.is_irreducible`

The field-relative is irreducible property of this scheme.

#### `AffineScheme.is_integral`

The field-relative is integral property of this scheme.

#### `AffineScheme.is_equidimensional`

The field-relative is equidimensional property of this scheme.

#### `AffineScheme.is_unmixed`

The field-relative is unmixed property of this scheme.

#### `AffineScheme.codimension`

The field-relative codimension property of this scheme.

#### `AffineScheme.tangent_cone`

```text
tangent_cone(self, point)
```

Return tangent cone at an exact point of this scheme.

#### `AffineScheme.local_dimension`

```text
local_dimension(self, point)
```

Return local dimension at an exact point of this scheme.

#### `AffineScheme.embedding_dimension`

```text
embedding_dimension(self, point)
```

Return embedding dimension at an exact point of this scheme.

#### `AffineScheme.multiplicity`

```text
multiplicity(self, point)
```

Return multiplicity at an exact point of this scheme.

#### `AffineScheme.tangent_space`

```text
tangent_space(self, point)
```

Return a basis of tangent directions in an affine chart.

#### `AffineScheme.is_normal`

Decide normality; projective schemes are checked on their affine charts.

#### `AffineScheme.is_smooth`

```text
is_smooth(self, point=None)
```

Decide smoothness globally or at a specified point.

#### `AffineScheme.singular_locus`

```text
singular_locus(self, *, reduced=False)
```

Return the singular subscheme, with automatic dimension checks.

#### `AffineScheme.point_search`

```text
point_search(self, height=None, *, max_candidates=1000000)
```

Enumerate finite-field points or QQ points within an explicit height bound.

#### `AffineScheme.point`

```text
point(self, values, *, field=None)
```

Validate a point over the base field or a specified exact extension.

#### `AffineScheme.subscheme`

```text
subscheme(self, *equations)
```

Impose additional equations with their scheme multiplicities.

#### `AffineScheme.intersection`

```text
intersection(self, other)
```

Return the scheme-theoretic intersection in the same affine space.

#### `AffineScheme.union`

```text
union(self, other)
```

Return the scheme-theoretic union in the same affine space.

#### `AffineScheme.base_change`

```text
base_change(self, field)
```

Extend the exact coefficient field through its declared embedding.

#### `AffineScheme.morphism`

```text
morphism(self, target, images)
```

Construct a polynomial morphism to an affine or projective target.

#### `AffineScheme.normalize`

```text
normalize(self, **limits)
```

Return the normalized affine scheme with its reduced integral-closure ring.

#### `AffineScheme.normalization_map`

```text
normalization_map(self, **limits)
```

Return the finite birational morphism from the normalization to this scheme.

#### `AffineScheme.local_ring`

```text
local_ring(self, point)
```

Return the local ring after extending coefficients to an exact point.

#### `AffineScheme.noether_normalization`

```text
noether_normalization(self, **options)
```

Construct a finite polynomial subring and the corresponding affine projection.

#### `AffineScheme.projective_closure`

```text
projective_closure(self, *, homogenizer=None)
```

Return the saturated projective closure with one appended coordinate.

#### `AffineScheme.product`

```text
product(self, other, *, variables=None)
```

Form the affine product, with disjoint coordinates in factor order.

#### `AffineScheme.principal_open`

```text
principal_open(self, polynomial, *, inverse=None)
```

Realize D(f) as Spec(S[u]/(I,u*f-1)), retaining an explicit inverse coordinate.

#### `AffineScheme.rational_points`

```text
rational_points(self, *, field=None, max_dimension=256, max_candidates=1000000)
```

Enumerate all points in a specified field for a zero-dimensional scheme.

Linear factors of an exact rational-univariate representation select precisely the points in that field. Multiplicity is retained in the scheme, not by repeating points. Positive-dimensional point searches and general Diophantine solving are outside this method.

#### `AffineScheme.reduced`

```text
reduced(self, **limits)
```

Return the reduced closed subscheme defined by the radical ideal.

#### `AffineScheme.irreducible_components`

```text
irreducible_components(self, **limits)
```

Return reduced irreducible components over the declared coefficient field.

### `ProjectiveScheme`

```text
class ProjectiveScheme(ideal)
```

A closed subscheme of ordinary projective space, with saturated homogeneous ideal.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L364)

Stored fields: `ideal`.

#### `ProjectiveScheme.ambient`

The homogeneous ambient polynomial ring.

#### `ProjectiveScheme.variables`

The ordered homogeneous coordinates.

#### `ProjectiveScheme.domain`

The coefficient field.

#### `ProjectiveScheme.homogeneous_coordinate_ring`

The graded ring S/I, whose dimension differs from dimension of Proj.

#### `ProjectiveScheme.dimension`

The geometric dimension of Proj; -1 denotes the empty scheme.

#### `ProjectiveScheme.degree`

The projective degree from the standard-graded Hilbert series.

#### `ProjectiveScheme.is_empty`

Whether the projective scheme is empty after irrelevant saturation.

#### `ProjectiveScheme.is_reduced`

The field-relative is reduced property of this scheme.

#### `ProjectiveScheme.is_irreducible`

The field-relative is irreducible property of this scheme.

#### `ProjectiveScheme.is_integral`

The field-relative is integral property of this scheme.

#### `ProjectiveScheme.is_equidimensional`

The field-relative is equidimensional property of this scheme.

#### `ProjectiveScheme.is_unmixed`

The field-relative is unmixed property of this scheme.

#### `ProjectiveScheme.codimension`

The field-relative codimension property of this scheme.

#### `ProjectiveScheme.tangent_cone`

```text
tangent_cone(self, point)
```

Return tangent cone at an exact point of this scheme.

#### `ProjectiveScheme.local_dimension`

```text
local_dimension(self, point)
```

Return local dimension at an exact point of this scheme.

#### `ProjectiveScheme.embedding_dimension`

```text
embedding_dimension(self, point)
```

Return embedding dimension at an exact point of this scheme.

#### `ProjectiveScheme.multiplicity`

```text
multiplicity(self, point)
```

Return multiplicity at an exact point of this scheme.

#### `ProjectiveScheme.tangent_space`

```text
tangent_space(self, point)
```

Return a basis of tangent directions in an affine chart.

#### `ProjectiveScheme.is_normal`

Decide normality; projective schemes are checked on their affine charts.

#### `ProjectiveScheme.is_smooth`

```text
is_smooth(self, point=None)
```

Decide smoothness globally or at a specified point.

#### `ProjectiveScheme.singular_locus`

```text
singular_locus(self, *, reduced=False)
```

Return the singular subscheme, with automatic dimension checks.

#### `ProjectiveScheme.point_search`

```text
point_search(self, height=None, *, max_candidates=1000000)
```

Enumerate finite-field points or QQ points within an explicit height bound.

#### `ProjectiveScheme.point`

```text
point(self, values, *, field=None)
```

Validate homogeneous coordinates and normalize their first nonzero value.

#### `ProjectiveScheme.affine_chart`

```text
affine_chart(self, index, *, variables=None)
```

Dehomogenize on the chart where the indexed coordinate is nonzero.

#### `ProjectiveScheme.base_change`

```text
base_change(self, field)
```

Extend the coefficient field and saturate the resulting homogeneous ideal.

#### `ProjectiveScheme.subscheme`

```text
subscheme(self, *equations)
```

Impose additional homogeneous equations and saturate at the irrelevant ideal.

#### `ProjectiveScheme.intersection`

```text
intersection(self, other)
```

Return the scheme-theoretic projective intersection.

#### `ProjectiveScheme.union`

```text
union(self, other)
```

Return the scheme-theoretic projective union.

#### `ProjectiveScheme.morphism`

```text
morphism(self, target, forms)
```

Construct a morphism given by homogeneous forms of equal degree.

#### `ProjectiveScheme.rational_points`

```text
rational_points(self, *, field=None, max_dimension=256, max_candidates=1000000)
```

Enumerate field-valued points of a zero-dimensional projective scheme by charts.

#### `ProjectiveScheme.reduced`

```text
reduced(self, **limits)
```

Return the reduced projective closed subscheme.

#### `ProjectiveScheme.irreducible_components`

```text
irreducible_components(self, **limits)
```

Return the reduced projective components over the declared coefficient field.

### `AffinePoint`

```text
class AffinePoint(scheme, coordinates, field=None)
```

An exact field-valued point of an affine scheme, with a specified embedding.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L325)

Stored fields: `scheme`, `coordinates`, `field`.

#### `AffinePoint.evaluate`

```text
evaluate(self, expression)
```

Evaluate a regular coordinate-ring element at this point.

#### `AffinePoint.change_field`

```text
change_field(self, field)
```

Embed this point's values into a specified larger exact field.

### `ProjectivePoint`

```text
class ProjectivePoint(scheme, coordinates, field=None)
```

A projective point over an exact field, with coordinates canonical up to scaling.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L579)

Stored fields: `scheme`, `coordinates`, `field`.

### `AffineChart`

```text
class AffineChart(projective_scheme, index, variables=None)
```

A projective affine patch with its embedding and transition coordinates.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L612)

Stored fields: `projective_scheme`, `index`, `variables`.

#### `AffineChart.embedding_coordinates`

Homogeneous coordinate expressions for the affine-chart embedding.

#### `AffineChart.scheme`

The affine scheme obtained by setting the selected homogeneous coordinate to one.

#### `AffineChart.embedding`

The polynomial morphism from this chart to its projective scheme.

#### `AffineChart.to_projective`

```text
to_projective(self, point)
```

Embed a point of this affine chart into projective space.

#### `AffineChart.from_projective`

```text
from_projective(self, point)
```

Convert a point in this open chart to its affine coordinates.

#### `AffineChart.transition_coordinates`

```text
transition_coordinates(self, other)
```

Return rational coordinates into another chart on their overlap.

#### `AffineChart.overlap_guard`

```text
overlap_guard(self, other)
```

The nonvanishing condition defining overlap with another chart.

### `AffineMorphism`

```text
class AffineMorphism(source, target, coordinates)
```

A polynomial morphism of affine schemes over one coefficient field.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L699)

Stored fields: `source`, `target`, `coordinates`.

#### `AffineMorphism.pullback`

The contravariant map from target to source coordinate rings.

#### `AffineMorphism.compose`

```text
compose(self, inner)
```

Return self after inner as a polynomial morphism.

#### `AffineMorphism.graph_ideal`

The graph ideal, with fresh target symbols following the source coordinates.

#### `AffineMorphism.image_closure`

The scheme-theoretic image closure, computed from the pullback kernel.

#### `AffineMorphism.is_dominant`

Decide whether the image is topologically dense in the target.

#### `AffineMorphism.is_finite`

Decide module finiteness over the affine target coordinate ring.

#### `AffineMorphism.generic_degree`

Return generic fiber length for a dominant generically finite affine morphism.

#### `AffineMorphism.ramification_locus`

```text
ramification_locus(self)
```

Return the non-unramified subscheme of a finite affine morphism.

#### `AffineMorphism.branch_locus`

```text
branch_locus(self)
```

Return the image closure of the relative non-unramified locus.

#### `AffineMorphism.constructible_image`

Return the exact geometric image with an existential fiber predicate.

#### `AffineMorphism.preimage`

```text
preimage(self, subscheme)
```

Return the scheme-theoretic inverse image of a target closed subscheme.

#### `AffineMorphism.fiber`

```text
fiber(self, point)
```

Return the exact fiber over a field-valued target point, including nilpotents.

#### `AffineMorphism.fiber_product`

```text
fiber_product(self, other, *, variables=None)
```

Form the affine fiber product of two morphisms to the same target.

#### `AffineMorphism.base_change`

```text
base_change(self, field)
```

Extend source, target and coordinate expressions to a larger exact field.

#### `AffineMorphism.inverse`

```text
inverse(self)
```

Return the inverse when the coordinate-ring pullback is an isomorphism.

### `ProjectiveMorphism`

```text
class ProjectiveMorphism(source, target, coordinates)
```

A morphism to ordinary projective space represented by global polynomial coordinates.

A projective source requires homogeneous forms of equal degree. Coordinates must have no common zero on the source. Affine sources allow arbitrary polynomial coordinates with no common affine zero. This presentation does not represent every line bundle or every abstract scheme morphism.

[Source](https://github.com/BhuvaneshBhatt/alggeom/blob/main/src/alggeom/schemes.py#L839)

Stored fields: `source`, `target`, `coordinates`.

#### `ProjectiveMorphism.compose`

```text
compose(self, inner)
```

Return self after another morphism into its projective source.

#### `ProjectiveMorphism.fiber`

```text
fiber(self, point)
```

Return the fiber over a projective point, as an affine or projective scheme.

#### `ProjectiveMorphism.preimage`

```text
preimage(self, subscheme)
```

Pull back the homogeneous equations of a projective closed subscheme.

#### `ProjectiveMorphism.graph_ideal`

The graph ideal in the product coordinates, saturated in both projective blocks.

#### `ProjectiveMorphism.image_closure`

The projective scheme-theoretic image closure by saturated graph elimination.

#### `ProjectiveMorphism.base_change`

```text
base_change(self, field)
```

Extend the coefficient field of this polynomial presentation.
