# Capabilities and limits

| Capability | Available interface | Boundary |
| --- | --- | --- |
| Prime parameter quotients and Gröbner covers | `parameter_fields.PrimeParameterField`, `parameter_fields.groebner_over_prime_parameters`, `parameter_cover.comprehensive_prime_cover` | Exact rational-coefficient parameter ideals with certified field-relative primality; implicit nonlinear exceptional branches, guard and membership certificates, symbolic partition verification; bounded by ideal decomposition and branch resources; not an absolute-primality decision procedure |
| Polynomial ideals | `ideals.PolynomialRing`, `ideals.Ideal` | Exact characteristic-zero and supported finite fields; global term orders; normal forms, witnesses, elimination, intersection, quotient, saturation |
| Coordinate rings and maps | `ideals.CoordinateRing`, `rings.CoordinateElement`, `rings.RingHomomorphism` | Exact arithmetic, units, kernels, contractions, image lifting; coefficient fields fixed by homomorphisms |
| Schemes and points | `schemes.AffineScheme`, `ProjectiveScheme`, `AffinePoint`, `ProjectivePoint` | Closed schemes over fields; complete finite-field enumeration in any dimension; characteristic-zero enumeration in dimension zero; bounded QQ searches; nilpotents retained |
| Charts, morphisms and fibers | `schemes.AffineChart`, `AffineMorphism`, `ProjectiveMorphism` | Projective charts and overlaps; polynomial and homogeneous morphisms, image closures, fibers and affine fiber products |
| Constructible geometry | `constructible.ConstructibleSet`, `ConstructibleImage` | Exact existential membership and finite locally closed Chevalley images for arbitrary polynomial and rational maps over QQ, subject to an explicit branch budget |
| Geometric integrality and factorization | `absolute_factorization.absolute_factor_bivariate`, `geometric.absolute_primary_decomposition`, `geometric.is_geometrically_integral` | Exact logarithmic-derivative absolute factors for rational-coefficient hypersurfaces in multiple variables, including multiplicities; selected nonprincipal positive-dimensional ideals, including disjoint-variable complete intersections, polynomial-graph primes, and smooth-rational-point geometric-primality certificates; certified algebraic constants for supported positive-dimensional prime function fields; arbitrary coupled nonprincipal decomposition remains unsupported |
| Parameter fibers | `families.fiber_stratification`, `comprehensive_strata.comprehensive_fiber_stratification` | Exact constructible fiber-dimension strata and optional bounded point counts over QQ |
| Specialization-aware Gröbner bases | `groebner_cover.comprehensive_groebner_cover` | Certified bases over one parameter and multivariate branches with rational linear specialization and some nonlinear exceptional hypersurfaces; some nested nonlinear exceptional strata with algebraic polynomial-graph parametrizations; general implicit exceptional varieties unsupported |
| Morphism invariants | `morphism_invariants.morphism_degree`, `ramification_locus`, `branch_locus` | Generic degree for dominant maps to affine space; graph-based finiteness of affine morphisms and Fitting-ideal non-unramified loci for finite affine maps |
| Radicals and components | `components.radical`, `minimal_primes`, `associated_primes`, `primary_decomposition` | Generic localization in characteristic zero; components over the declared field; embedded primary factors included |
| Ring and scheme predicates | `properties.is_reduced`, `is_irreducible`, `is_integral`, `is_equidimensional`, `is_unmixed`, `is_smooth` | Field-relative components; pointwise smoothness from local dimension and tangent rank; finite-characteristic component limits apply |
| Rational functions and local rings | `rational_functions.FractionField`, `Localization`, `LocalRing` | Fraction fields of domains, principal localizations, local rings at certified primes, exact arithmetic and regularity; no completions |
| Local invariants | `local.tangent_cone`, `local_dimension`, `embedding_dimension`, `multiplicity`, `module_length`, `intersection_multiplicity` | Exact field-valued points; finite-support module lengths; proper isolated local intersections with one complete-intersection presentation; primary components away from the point are removed |
| Noether normalization | `noether.noether_normalization` | Certified finite polynomial subring, explicit coordinate automorphisms, monic integral relations, spanning module generators and generic rank; polynomial changes over small finite fields |
| Affine and projective normalization charts | `normalization.integral_closure`, `normalize`, `normalization_map`, `projective_normalization.projective_normalization_atlas` | Reduced affine algebras in characteristic zero; projective curve normalization atlas has explicit regular transition coordinates and verified cocycles, but no single projective embedding |
| Finite coefficient fields | `fields.finite_field`, `field_order` | Validated GF(p) and GF(p^d), exact polynomial/module/Hilbert arithmetic and resolutions; Frobenius radicals; dimension-zero primary decomposition |
| Point searches | `points.point_search`, scheme `rational_points` | Complete finite-field enumeration with budget; bounded-height QQ searches do not prove absence of rational points |
| Plane curves | `curves.PlaneCurve`, `curves.genus`, `curves.parametrize` | Reduced irreducible plane models over QQ/number fields; normalization analysis and rational P^1 normalization maps |
| Graphs and images | `maps.RationalMap`, `maps.implicitize`, `maps.graph_ideal` | Saturated Zariski closures; input exclusions retained separately |
| Tangents and singularities | `maps.tangent_space`, `maps.tangent_ideal`, `maps.singular_locus` | General exact-point tangents; equidimensional Jacobian singular subschemes; reduced singular sets for mixed-dimensional inputs |
| Hilbert invariants | `graded.hilbert_series`, `graded.hilbert_function`, `graded.hilbert_polynomial`, `graded.degree` | Standard grading of S/I; no weighted or multigraded series |
| Projective closure | `graded.projective_closure`, `graded.projective_dimension` | Homogenization plus saturation; projective dimension differs from ring dimension |
| Modules | `modules.Submodule`, `modules.PresentedModule`, `modules.syzygies` | Elements, maps, kernels/images/cokernels, direct sums, tensor products, Hom, annihilators and exact finite presentations over S or S/I |
| Graded resolutions | `resolutions.free_resolution`, `resolutions.GradedResolution` | Minimal standard-graded resolutions over polynomial rings; Betti tables, regularity and module projective dimension; marked prefixes under a length limit |
| Differentials and Fitting ideals | `modules.differentials`, `modules.differential`, `modules.fitting_ideal` | Relative to the coefficient field; Fitting ideals over S/I returned as inverse images in S |
| Common right components | `decomposition.common_right_component`, `common.common_components` | Generator of a joint transcendence-degree-one rational field, with input exclusions and outer maps |
| Rational two-map decomposition | `rational.search_rational`, `decomposition.decompositions` | QQ and exact declared number fields; both degrees greater than one |
| Multivariate input | Same rational interface | One scalar intermediate and a univariate outer map |
| Recursive decomposition | `chains.rational_chain`, `chains.algebraic_chain` | One selected chain, with evidence on each leaf |
| Automatic algebraic factor extensions | `factors.discover_factors` | Requested factor degrees; finite algebra dimensions grow as binomial coefficients |
| Implicit coefficient-field reconstruction | `discovery.search_algebraic`, `algebraic.reconstruct_algebraic` | Rational coefficient curves supported by the backend over the declared constants |
| Parametric polynomial components | `parameters.polynomial_pair`, `parameters.polynomial_family` | Real parameters and polynomial components in characteristic zero |
| Branch, pole and domain checks | `real.verify_decomposition`, `real.verify_chain`; generic graphs via `semialg.function_verification` | Exact supported semialgebraic formulas in affine real coordinates |
| Number-field coefficients | `constants=QQ.algebraic_field(alpha)` or a ring's `domain` | Native coefficient embeddings over K and K(z); tested real and imaginary quadratic, cubic and composite fields |

## Mathematical scope

### Function decomposition

The rational search implements the uni-multivariate notion of Gutiérrez, Rubio and Sevilla: `f(X)=g(h(X))` with both rational degrees greater than one. In particular, `(x+y)**2` has no nontrivial split under this convention because `x+y` has total degree one. The search is complete over the declared constant field and returns classes modulo fractional-linear changes of the intermediate coordinate. It does not search vector-valued intermediate maps or characteristic-positive fields.

Algebraic search examines minimal fields generated by coefficients of factors. A supported negative result means that this construction has no accepted component for every covered orbit. It does not rule out every algebraic correspondence, every constant extension, or every alternative chain. Recursive rational graph splits extend this scope for relations linear in either coordinate.

### Curves and real domains

A genus-zero curve need not have a rational parameter over the declared ground field: a ground-field point is also needed. The [public curve interface](curves.md) uses algroots normalization analysis and parametrization with explicit model and resource limits. The curve analysis certificate supplies singularity/adjoint evidence. Actual affine integral-closure rings and finite birational maps are available separately through [normalization](normalization.md). An unsupported operation or an exhausted reconstruction remains unresolved. The package does not turn an unresolved diagnostic into an indecomposability proof.

Resultant identities describe algebraic composition and can include multiplicity. They do not by themselves select real branches or guarantee equality of affine domains. Real verification retains visible input denominators and internal component poles. It cannot reconstruct exclusions that were simplified away before the expression was passed in.

### Geometric structures and coefficient fields

The [geometry tutorial](geometry.md) specifies Krull dimension, vector-space length, graded-piece Hilbert polynomials and image-closure semantics. General scheme gluing and sheaves, a single global projective embedding of positive-genus normalizations, weighted/multigraded resolutions, quotient-ring resolutions, general Ext and Tor are not implemented. The intersection interface computes its required Koszul Tor groups internally. See [schemes](schemes.md), [components](components.md), [resolutions](resolutions.md) and [curves](curves.md) for their precise supported cases.

Singular-locus calculations infer equidimensionality from minimal primes when needed. Mixed-dimensional inputs require `reduced=True`; this returns their reduced singular set. Finite-field positive-dimensional primary decomposition and prime certification remain unsupported, so predicates or constructions needing them can raise. Rational-function fields in positive characteristic support polynomial arithmetic but not general radical, primary or normalization algorithms.

### Parameter fibers and exceptional loci

Generic calculations over K(t) do not retain exceptional parameter strata automatically, while `comprehensive_fiber_stratification` provides exact constructible dimension and optionally bounded geometric cardinality strata for coupled nonlinear equations over QQ. It does not compute component strata or singularity strata. The separate Gröbner-cover interface supplies polynomial bases on supported parameter branches. `families.fiber_stratification` also computes exact nonlinear hypersurface and independent-product strata. Coupled nonlinear systems have exact geometric dimension strata, and optionally bounded point-count strata, through [constructive elimination](comprehensive_fibers.md); prime-quotient Gröbner covers also handle supported implicit nonlinear exceptional parameter varieties with replayable identities and a symbolic partition proof; resource-limited cases still raise and no unconditional general termination bound is claimed. Constructive images of arbitrary polynomial and rational maps over QQ are available through [Chevalley projection](constructive_chevalley.md), subject to a branch limit. The [number-field guide](number_fields.md) covers tested embeddings and native conversions.

## Resource limits

| Limit | Applies to | Effect |
| --- | --- | --- |
| `max_degree` | Rational total degree or implicit output degree | Oversized two-map searches raise; chains retain an open leaf |
| `max_candidates` | Candidate products or discovered factor orbits | Stores the examined prefix and leaves coverage incomplete |
| `max_factor_dimension` | Universal monic factor algebra | Skips and records oversized factor degrees before constructing them |
| `max_separating_attempts` | Linear forms separating the finite factor algebra | Records unresolved degrees on exhaustion |
| `max_steps`, `max_dimension` | Curve and reconstruction backend | Backend limits remain visible as unresolved work or backend errors |
| `max_depth`, `max_nodes` | Recursive chain | Retains exact components with explicit open leaves |
| `max_dimension`, `max_nodes`, `max_exponent`, `max_separating_attempts` | Ideal primary/radical decomposition | Raises on exhaustion; no partial tuple is returned |
| `max_length` | Graded resolution | Returns a marked prefix when a higher syzygy is needed; complete invariants raise |
| `max_candidates` | Finite-field enumeration and bounded-height point search | Raises before exceeding the ambient candidate budget |
| `max_linear_attempts`, `max_transform_degree`, `max_generators` | Noether normalization | Bounds linear search, polynomial coordinate degree and finite module spanning set; raises on exhaustion |
| `max_steps`, `max_new_variables`, `max_module_generators`, `max_components` | Affine integral closure | Raises before returning an intermediate extension as a completed normalization |
| `max_vector_dimension` | Finite-support module length | Bounds the matrix dimension of the finite coefficient-space calculation |
| `max_rational_degree` | Rational graphs inside an implicit chain | Bounds the separate graph decomposition engine |

Defaults target small exact examples. Candidate budgets do not impose a wall-clock or memory limit on factorization, Groebner bases, resultants, quantifier elimination, or certificate replay. These operations can be expensive before any candidate is examined. Arbitrary parametric rational or implicit decomposition is outside the parameter-family API.

## Gröbner construction and verification

Ideal and quotient construction use the native algroots backend by default. `Ideal.groebner_basis_from("native", limits=...)` requests a separately computed basis and verifies it in the declared field. A portable algroots installation requires explicit `groebner_backend="sympy"` selection on the ideal or quotient factory. Algebraic and parameter fields retain their exact reconstruction and verification contracts. See [backend selection](backends.md) for installation, availability errors and the scope of resource limits.
