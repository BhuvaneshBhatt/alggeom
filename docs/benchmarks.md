# Performance benchmarks

Use the benchmark runner to measure exact algebraic geometry and decomposition costs on your environment and to compare changes to the implementation. The suite records structural input sizes as well as elapsed time: a candidate count alone does not capture factorization, finite-algebra construction, or real elimination costs.

Install matching backends and the package, then run from the checkout:

```sh
python benchmarks/run.py --suite smoke --repeat 3 --output benchmark.json
python benchmarks/run.py --suite standard --real --repeat 5 --output benchmark-real.json
SYMPY_GROUND_TYPES=python python benchmarks/run.py --suite smoke --real --output benchmark-python.json
```

The real workload requires semialg. Without `--real`, it is listed under `omitted_cases`; no zero timing is substituted. The manual [benchmark workflow](../.github/workflows/benchmarks.yml) runs either suite on Linux with the selected coefficient backend and uploads its JSON report. CI also runs the smoke suite once to check that the workloads remain executable.

| Case | Measured operation | Smoke input | Standard input |
| --- | --- | --- | --- |
| `rational_search` | Difference factorization, candidate enumeration and recovery | `x⁸` | `x¹²` |
| `number_field_search` | Rational search over QQ(sqrt(2)) | `x²+2/x²`, three classes | Same field-sensitive input |
| `factor_discovery` | Construct the universal quadratic-factor algebra and its orbit fields | Cubic, algebra dimension 3 | Quartic, algebra dimension 6 |
| `reconstruction` | Coefficient-field construction, rationality and implicit components from a supplied factor | Quartic in a degree-two extension | Same reconstruction input |
| `certificate_replay` | Replay a factor discovery certificate, including coverage and inverse coordinates | Cubic quadratic factors | Same replay input |
| `rational_chain` | Build a recursive chain and leaf evidence | `x⁸`, three components | `x¹⁶`, four components |
| `absolute_monomial` | Absolute primary decomposition of a monomial ideal with an embedded component | `(x²,xy)` | Same input |
| `absolute_graph` | Geometric primality of a polynomial graph | `z=xy` | Same input |
| `groebner_generic` | Specialization-aware Gröbner cover with a linear exceptional branch | `(ax,x+y)` | Same input |
| `groebner_nested` | Gröbner cover with intersecting nonlinear exceptional graphs | `((a²+b²)x,ay)` | Same input |
| `primary_components` | Generic localization and embedded primary factors | `(x²,xy)`, two primary components | Same exact embedded-component problem |
| `graded_resolution` | Minimal graded syzygies and Betti shifts | Twisted cubic, shifts 0; 2,2,2; 3,3 | Same exact resolution |
| `local_tangent_cone` | Translate, saturate and specialize the complete local initial ideal | `(x²+y³,x²-y³)` at the origin | Same local length-six problem |
| `noether_projection` | Certified monic coordinate changes and finite module generators | `xy=1`, generic rank two | Same finite projection |
| `integral_closure` | Endomorphism-ring extensions and birational coordinates | Cusp `y²=x³` | Same conductor `(x,y)` problem |
| `finite_points` | Complete affine enumeration and equation checks over GF(9) | Parabola with 81 ambient candidates | Same exact nine-point set |
| `real_branches` | Totality and uniqueness of a selected square-root graph | `y²=x`, `y≥0`, domain `x≥0` | Same exact branch question |

## Measurement protocol

Each case runs in a fresh subprocess. Preparation is timed separately as `setup_seconds`, then the runner performs one unrecorded warmup and three recorded repetitions by default. Every result must pass its mathematical checks after the clock stops. The JSON contains all samples, their median and minimum, dimensions, package versions, Python version, operating system, processor information, arithmetic backend, and hashes of the benchmark code and installed alggeom and backend Python sources. Those hashes distinguish workload changes and patched sources that share a distribution version.

Repeated calls within a case retain normal SymPy and backend caches. Certificate replay reuses the prepared certificate; reconstruction reuses its ambient field. These timings measure repeated use, not interpreter startup or an entirely cold process for every sample. Setup can itself construct a sample search or certificate, so it should not be interpreted as import time. Real workloads measure exact solver execution, not sampled numerical comparisons. A package imported from an uninstalled source tree has a null distribution version; its source fingerprint is still recorded.

`--case factor_discovery` restricts the run to one case; the option can be repeated. `--warmup 0` disables the unrecorded repetitions. `--timeout 120` limits each complete case process, including preparation and verification. Timed-out or failed cases receive an explicit status, have no fabricated samples, and make the command exit unsuccessfully. These limits are separate from the package's degree and candidate budgets.

## Comparing results

Compare reports with the same suite, coefficient backend, interpreter, dependencies, source profile, machine and repetition settings. Keep the full sample distributions; a single fast or slow run is weak evidence. For a code change, repeat the same protocol on both revisions and inspect which structural operation changed. CI smoke results are an execution check, not a performance regression threshold.

The suites are small reproducible workloads, not an asymptotic complexity study or a general performance claim. Universal factor-algebra dimension grows as `binomial(n,d)`, so larger output degrees can become expensive before reconstruction begins. The benchmark runner is separate from the correctness suite, and the tests impose no wall-clock speed assertions.

## Comparing source changes

Run the same workload definitions on both source trees on one machine, with the same interpreter, dependencies, coefficient backend, warmup count, and repetitions. Record source fingerprints and all raw samples. Verify mathematical results outside each timed operation. Keep setup, native kernel work, exact verification, process startup, and total wall time distinct.

The benchmark runner writes JSON reports to the output location selected by the caller. The package does not ship historical measurements as current performance claims. Use `benchmarks/run.py` for geometry workflows and algroots' `benchmarks/bench_verified_backends.py` for verified Gröbner, quotient, and solver construction.
