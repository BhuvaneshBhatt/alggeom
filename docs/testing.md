# Testing and validation

Install matching backends and the `real` and `test` extras, then run these commands from the package directory:

```sh
python tools/install_backends.py --published --profile real
python -m pip install '.[real,test]'
python -m pytest -p no:cacheprovider
ruff check --no-cache .
ruff format --no-cache --check .
python tools/generate_api.py --check
SYMPY_GROUND_TYPES=python python -m pytest -p no:cacheprovider
python -m build
```

## What the tests establish

| Test area | Independent mathematical expectation | Main failure modes |
| --- | --- | --- |
| Ideals | Explicit monomial intersections, quotients and multi-generator saturations; original-generator witnesses | Confusing saturation by an ideal with saturation by a product; incorrect elimination contraction |
| Coordinate rings | Known Krull dimensions, Artinian lengths, coefficient-field and zero-ring cases | Treating positive-dimensional rings as finite quotient algebras |
| Scheme geometry | Double fibers, extension-field points, conic chart transitions, Veronese images and base-locus rejection | Losing nilpotents, confusing rational-point absence with emptiness, using unsaturated projective graphs |
| Components | Known primary intersections, embedded and mixed-dimensional primes, nonmonogenic nilpotents, field-relative factorization | Dropping embedded structure, incorrect generic contraction, treating resource exhaustion as an answer |
| Absolute geometric components | Independent reconstruction of primary intersections, disjoint-variable factors, polynomial graphs and nonreduced supports | Treating field-relative primeness as geometric primeness, or a reconstructed intersection as proof of primariness |
| Gröbner covers | Specialized exact Gröbner bases at algebraic parameter tuples on generic and exceptional branches; bounded piece counts | Missing exceptional branches, incorrect leading monomials after specialization, overlapping pieces and exceeding budgets |
| Graded resolutions | Koszul and complete-intersection shifts, twisted-cubic Betti numbers, independent Hilbert series and exactness replay | Confusing a complex with a resolution, wrong grading, nonminimal scalar entries, false completeness |
| Local geometry and fractions | Tangent-cone cancellations, mixed dimensions, local multiplicity, local Tor lengths, zero divisors and regular fraction cancellation | Confusing tangent dimension with local dimension, or cancelling a zero divisor |
| Noether normalization and integral closure | Finite projections over QQ/number/finite fields, nonlinear finite-field coordinates, curve and surface integral closures, product normalizations, certificate mutations and bounds | Confusing generic finiteness with finiteness, spanning sets with free bases, or partial extensions with a normalization |
| Module arithmetic | Quotient module elements, relation validation, kernel/image/cokernel, tensor, Hom and empty modules | Losing quotient relations or assuming generator matrices are unique |
| Finite fields and point searches | Frobenius nilradicals with zero trace, nonmonogenic products, extension-field arithmetic, point sets and candidate limits | Using a trace kernel in bad characteristic or calling a bounded QQ search complete |
| Plane curves | Cusp delta/genus, projective normalization maps, conic obstructions and quadratic/cubic coefficient fields | Calling a missing point proof of nonrationality, losing infinity or algebraic constants |
| Maps and diagnostics | Cusp and hyperbola equations, canceled poles, explicit Jacobian kernels | Spurious graph components, lost exclusions, incorrect source equations |
| Graded invariants | Twisted-cubic Hilbert polynomial and brute monomial counts | Affine/projective dimension confusion, unsaturated homogenization, incorrect Artinian convention |
| Modules | Polynomial kernel witnesses, nilpotent quotient kernels, basis-invariant Fitting ideals and Leibniz rule | Ignoring base-ring equations, row/column convention errors, incorrect differential relations |
| Rational decompositions | Exact substitution and small candidate-count oracle | Missing classes, incorrect degree filters, false negatives |
| Polynomial chains | SymPy's separate polynomial algorithm plus direct composition | Incorrect complete chains or degree factors |
| Generated functions | Known compositions and explicit bidirectional field coordinates | Dependence on normal forms, failure under fractional-linear changes |
| Implicit components | Exact resultants including powers and units | Extra components, wrong multiplicity, altered coefficient embeddings |
| Search evidence | Complete orbit accounting and altered certificates | Incomplete work mistaken for indecomposability |
| Real graphs | Domain, pole, branch and coordinate-scope obligations | Canceled holes, unattainable branches, variable capture |
| Parameter families | Exact strata and independent specializations | Degree drops, undefined coefficients, missed exceptional equations |
| Applications | Known roots, implicit curve inverse, response equation, coefficient criteria | Correct low-level identities used incorrectly in a workflow |
| Documentation | Execute every Python fence and script; resolve links into the project | Stale snippets, broken source links, mismatched examples |
| Optional dependencies | Fresh interpreter blocks semialg and backend decomposition modules | Core imports accidentally depending on the real extra or copied entry points |
| Number fields | Prescribed quadratic intermediate fields, explicit cubic/composite embeddings and independently checked resultants | Lost algebraic embeddings, incorrect constant-field completeness, unsupported quotient conversions |
| Certificate mutations | Systematically alter identities, coverage, coordinate maps and proof trees | False acceptance of omitted orbits, changed scalar/power, substituted fields, fabricated stopping reasons |
| Development tools | Generated reference equality, dependency requirement consistency, executable benchmark JSON | Stale signatures/docstrings, mismatched backend inputs, failed workloads presented as measurements |

The polynomial comparison uses [SymPy's polynomial decomposition API](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.decompose) over QQ. Both implementations must compose to the input and produce the same multiset of degrees for a complete polynomial chain in characteristic zero. It does not compare printed components, demand a particular component order, or treat SymPy's one polynomial chain as an oracle for all rational two-map classes. Rational search is independently exercised on known intermediate fields and fractional-linear changes of the outer map.

The generated tests use bounded exact inputs and reproducible settings. Numerical sampling is not used to prove an identity or an exhaustive negative result. Examples are executed once per test session and their namespaces are shared with the application checks. Python snippets in each Markdown document share a namespace in document order, matching the way a reader follows a tutorial.

The number-field corpus includes imaginary quadratic rational searches and implicit reconstruction, real cubic and composite fields, automatic factor coefficient extensions, and real embedding checks. Quotients over K(z) check multiplication, inversion, traces, norms and minimal polynomials against prescribed coefficients. The imaginary-field conversion regression is a passing test; see the [number-field guide](number_fields.md).

## Algebraic-geometry certificates and independent checks

The parameter-cover proof layer verifies exact polynomial identities, zero-remainder basis-membership witnesses, standard-representation leading-monomial bounds for every critical pair, nonvanishing of coefficient denominators on the entire stratum, and correctness of algebraic parameter bindings throughout their claimed loci. Constructible-set differences and intersections establish coverage and pairwise disjointness over an algebraic closure; no finite set of parameter samples can substitute for this check. `tests/test_groebner_proofs.py` changes stored remainders, binding conditions and denominator guards to ensure invalid proofs are rejected. It also tests the intersection of the nonlinear exceptional graphs `a = ±ib` at the origin.

The optional Macaulay2 tests call independently implemented primary decomposition, associated primes, primariness and prime-radical checks. Examples include nonprincipal ideals with embedded and nonreduced components, and components appearing only after extension to `QQ(i)`. Each oracle fixture also has certifiably affine-linear radical supports: for these fixtures, the linear quotient presentation establishes geometric primality. Prime over a number field alone would **not** establish geometric primality for arbitrary nonlinear ideals. The distinction between prime over `QQ` and remaining prime after base extension is tested explicitly. Embedded primary components need not be uniquely determined, so equality of particular components is a stronger property of these fixtures, not a general uniqueness assumption.

Metamorphic tests transport ideals and primary components through invertible linear coordinate changes and exact field extensions. They check Krull dimension, Hilbert series and functions, graded degree, local multiplicity, and initial monomial ideals when variables *and the term order* are permuted together. Initial ideals generally change under arbitrary linear coordinate changes, so that is not asserted as an invariant.

These checks establish the mathematical obligations tested by the supported algorithms. They do not prove that unsupported parameter-ideal branches or arbitrary associated-prime decompositions are implemented. External oracle tests are optional, and a missing external executable yields a skip, not an oracle success.

### Additional adversarial tests

- Expand field extensions beyond quadratics, including fields with nontrivial automorphism groups and embedded primes of different dimensions.
- Exercise certificates containing chained rational parameter bindings and several interacting nonlinear exceptional ideals; check that a modified branch condition is rejected.
- Test constructible complements involving unequal-dimensional components and specializations where both leading coefficients and denominator factors vanish.
- Extend independent oracles to difficult singular and nonradical positive-dimensional ideals, comparing support primes separately from nonunique embedded primary components.
- Run the complete core/real backend matrix and timed performance comparisons before release; selected passing suites are not evidence that unrun tests pass.

## Backend ownership and adapter tests

`tests/test_algroots_adapter.py` checks delegation to algroots: coefficient characteristic, field embeddings, backend and budget propagation, and independently computed reference bases. `tests/test_basis_verification.py` and the witness/replay tests check alggeom's geometric certificate boundary. They do not implement or benchmark a second native Gröbner kernel. F4 scheduling, monomial storage and matrix elimination belong to algroots and its test suite.

## Independent Singular Gröbner-cover comparisons

`tests/test_singular_cover_oracle.py` compares whole specialized ideals with independently computed standard bases in Singular at rational parameter values, including intersections of exceptional conditions. It checks inclusion in both directions using Singular's polynomial remainder operation. This is an independent oracle for those samples, **not** a proof that the entire constructible stratum is correct; that obligation belongs to the exact certificate and partition verifier.

Set `SINGULAR_BIN` to the executable path and run `python -m pytest -m oracle tests/test_singular_cover_oracle.py`. If the executable is unavailable, these tests skip. When using a standalone executable packaged with shared libraries, configure `LD_LIBRARY_PATH` for that installation.

The symbolic partition verifier first tries a polynomial-factor Boolean proof, valid because a product vanishes over a field exactly when one of its factors vanishes. A successful propositional partition proof implies geometric coverage and disjointness. If factor relations are insufficient to prove the partition, the verifier retains the original ideal-theoretic constructible-set calculation. Mutation tests must exercise both the shortcut and the fallback.

## Continuous integration

The [test workflow](../.github/workflows/tests.yml) configures 24 FLINT jobs: Python 3.11, 3.12, 3.13 and 3.14 on Linux, macOS and Windows, each with core and real dependency profiles. Two Linux/Python 3.11 jobs also constrain the arithmetic dependencies to the minimum supported versions, one for each profile. Four further real-profile jobs select Python coefficients through `SYMPY_GROUND_TYPES=python`: Linux on Python 3.11 and 3.14, plus macOS and Windows on Python 3.14. The quality job checks Ruff, generated documentation, benchmark execution, wheel/sdist construction, installed examples and a wheel rebuilt from the source distribution. A configured job is not evidence that it has run successfully; inspect its actual GitHub Actions result for that claim.

Each job installs the [published backend requirements](../ci/backends.toml), without local source overrides. This exercises the dependency versions users can install from the index. The [backend guide](backends.md#install-matching-dependencies) describes the source-root option and dependency profiles. `tools/check_backends.py` verifies the active coefficient backend, exact function-field arithmetic, and requested real-adapter availability before tests run. Installed-wheel checks execute outside the source checkout.

Core jobs install algroots without semialg and run `python -m pytest --without-real`. That option excludes `tests/real` before collection and skips marked real examples and documentation. This ensures optional modules are never imported just to discover tests. Real jobs require semialg, run the full suite, and fail if it is missing. Core import tests additionally block the old decomposition modules in the backend so accidental delegation cannot pass unnoticed.

For a core local environment:

```sh
python tools/install_backends.py --published --profile core
python -m pip install '.[test]'
python tools/check_backends.py --profile core --ground-types flint
python -m pytest --without-real -p no:cacheprovider
```

## Certificate mutation coverage

`tests/test_mutations.py` starts with verified certificates and changes one mathematical or structural obligation at a time. The table covers rational identities and candidate counts; universal factor coordinates, separating weights and orbit coverage; coefficient generators and rationality evidence; resultant powers and units; reconstruction attempt ordering; and nested rational/implicit trees. Every changed record must return `False` from verification. The tests exercise altered degree, node and depth budgets, missing children, branch coordinates and field declarations as well as arithmetic data.

Open leaves for depth, node, degree and candidate limits are tested separately from repeated components and unsupported diagnostics. An unsupported diagnostic can be a valid unresolved record; the tests require that it remain incomplete and cannot be relabeled as indecomposable. Resolution mutations also alter differentials, augmentation and completeness; verification must reject incorrect exactness, grading and surjectivity. These are semantic mutations of stored certificates, not a source-code mutation score. They do not establish that every possible malformed object has been enumerated.

## Generated documentation

`python tools/generate_api.py` imports each public module and reads its `__all__`, signatures, docstrings, dataclass fields, properties, cached properties and inherited public methods. It writes [api.md](api.md), including links to defining source lines. The complete reference requires the real extra because it includes the optional APIs. CI and a real-profile test run `--check`, which fails if the checked-in reference differs. The generator works against either the source tree or an installed wheel.

The [coefficient-field tutorial](coefficient_fields.md) and [number-field guide](number_fields.md) have executable Python blocks. Record real-only examples in `examples/catalog.json` and mark real-only tutorial documents in `tests/test_documentation.py`; core documents remain collected in either profile.

## Distribution checks

Build the wheel and install it with matching backends in an environment outside the checkout. Run `python tools/check_installation.py --profile real` with source paths removed from `PYTHONPATH`. It rejects an import from `src`, executes every example script and all documented Python snippets, and reports the installed path. Use `--profile core` for the core workflows and a separate interpreter where semialg is unavailable. The source distribution also includes workflows, backend inputs, benchmark tools, release notes and the API generator. CI rebuilds a wheel from that source distribution, installs it and repeats the installed workflow checks. Follow [release checks](releasing.md) to run these steps locally.

The documentation link checks resolve relative links and this project's absolute GitHub source links against local files. They do not contact GitHub, establish that the proposed remote repository exists, or check third-party website uptime. Documentation and catalog text is decoded explicitly as UTF-8, and a regression test simulates Windows legacy text defaults with a non-ASCII executable tutorial.

The application figure is stored as SVG and can be regenerated with `python -m pip install ".[docs]"` followed by `python tools/plot_applications.py`. Plotting dependencies are separate from the computation and test dependencies.

## Performance measurements

The [benchmark guide](benchmarks.md) documents process isolation, warmups, sample distributions, structural dimensions, backend source fingerprints and per-case timeouts. Performance reports are separate from correctness tests. Larger field lattices and curve models can extend this corpus as the backend grows; include an independent mathematical expectation before adding a new measured workload.

## Installed examples, arithmetic floors and field policies

The example catalog in `examples/catalog.json` accounts for every runnable application. `tools/check_examples.py` rejects source imports and executes scripts in fresh interpreters and temporary directories, with assertions enabled and per-workflow timeouts. All alggeom Python tutorial fences execute. The backends mark executable tutorial fences explicitly; displayed signatures remain reference text.

CI exercises the supported effective arithmetic floor in `ci/minimum.txt` and independently resolved latest dependencies. The lower floor is SymPy 1.14.0 and python-flint 0.9.0, with matching backend inputs where needed. Python coefficient jobs test the conversion layer independently of FLINT ground coefficients. Local checks do not imply that remote matrix jobs have run.

[Certificate-field audits](certificate_audits.md) classify proof data and non-proof metadata explicitly and produce machine-readable outcomes. Coverage, disjointness, exact reconstruction, module map identities and known finite point sets supply independent mathematical expectations. Timing distributions and structural input sizes belong in benchmark reports, not correctness thresholds.

## External Macaulay2 oracle

`tests/test_macaulay2_primary_oracle.py` compares alggeom's absolute primary components with the independent `primaryDecomposition`, `associatedPrimes`, and `isPrimary` implementations in Macaulay2. The cases include nonprincipal ideals, embedded components, and splitting over the exact extension `QQ[a]/(a^2+1)`. The test checks primariness and ideal reconstruction as well as component matching, and does not infer correctness from reconstruction alone.

The test skips when the `M2` executable is unavailable. Set `M2_EXECUTABLE` to a Macaulay2 executable outside `PATH`; when running an extracted Macaulay2 distribution, set `M2_PREFIX` to its installation prefix and ensure the executable's shared libraries are discoverable. The test executes Macaulay2 through standard input and cleans up its process group after each run, so it does not require a writable installation. The especially expensive nonreduced geometric test is opt-in with `ALGGEOM_SLOW_ORACLE=1`. A skipped oracle is never counted as an agreement.

The companion Gröbner-cover tests check whole-parameter-space coverage, pairwise disjointness, replayed source and S-polynomial identities, and rejection of modified certificates. Metamorphic primary-decomposition tests check invariance under invertible linear coordinate changes and constant-field extension. These are complementary to independent computer-algebra oracles and do not replace general mathematical completeness proofs.
