# alggeom

Exact computational algebraic geometry for Python. Use it to turn parametrizations into equations, retain multiplicity in fibers, analyze singularities, compute finite projections and normalize affine models. Modules and graded resolutions describe relations; function decomposition exposes hidden intermediate quantities.

<!-- examples: executable -->

## Installation

Python 3.11 or later is required. Install the core package from PyPI:

```sh
python -m pip install alggeom
```

For real inequalities, selected branches and semialgebraic graph verification, install the optional real-geometry extra:

```sh
python -m pip install "alggeom[real]"
```

The core package requires algroots 0.4.0 or newer, SymPy 1.14 within the 1.x series, and python-flint 0.9 within the 0.x series. The real extra additionally requires semialg 1.5.0 or newer. Core algebraic geometry works without semialg.

The default Gröbner backend requires a native algroots installation. Compatible algroots wheels include the compiled kernels and their numerical libraries. If pip installs a portable algroots distribution, select the SymPy backend explicitly for ideal and quotient construction as shown below. See [backend installation and selection](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/backends.md) for source builds and dependency profiles.

## Quick start

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.schemes import AffineScheme
from alggeom.normalization import integral_closure

x, y = sp.symbols("x y")
S = PolynomialRing((x, y))
cusp = AffineScheme(S.ideal(y**2 - x**3))
assert cusp.coordinate_ring.krull_dimension == 1
assert cusp.multiplicity((0, 0)) == 2
assert not cusp.is_smooth((0, 0))
closure = integral_closure(cusp)
assert closure.mapping.is_injective and closure.verify()
assert closure.conductor == S.ideal(x, y)
```

The cusp has one-dimensional local geometry and a two-dimensional tangent space at the origin. Its normalization adds the integral function y/x. This is different from Noether normalization, which finds a polynomial subring over which the existing ring is finite.

## Package boundaries

[algroots](https://github.com/BhuvaneshBhatt/algroots) supplies Gröbner-basis construction, finite quotient and coefficient-field arithmetic. [semialg](https://github.com/BhuvaneshBhatt/semialg) supplies exact real decisions and branch conditions. alggeom owns general coordinate rings, schemes, normalization and decomposition workflows. See [package ownership](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/package_ownership.md).

## Capabilities and limits

| Task | Result | Scope |
| --- | --- | --- |
| Ideals and coordinate rings | Exact normal forms, ideal operations, elements and homomorphisms | Supported exact fields; Krull dimension is distinct from finite vector-space dimension |
| Schemes and maps | Affine/projective points, fibers, constructible sets, and exact existential images | Exact locally closed images of arbitrary affine polynomial and rational maps over QQ, with an explicit branch budget |
| Components and local geometry | Radicals, primary factors, absolute primary factors of hypersurfaces, selected finite algebras, separated complete intersections, and polynomial-graph prime ideals, tangent cones and local intersection multiplicity | Field-relative primary decomposition and certified QQ hypersurface absolute splitting; isolated proper local intersections |
| Modules and resolutions | Kernel/image/cokernel, tensor, Hom, Fitting ideals and Betti tables | Finite presentations; standard-graded polynomial resolutions |
| Noether normalization | Certified finite projection and monic integral relations | Spanning module generators need not be a free basis |
| Normalization | Affine finite birational maps, projective normalization patches and domain conductors | Reduced characteristic-zero affine algebras; regular projective-chart transitions and cocycle verification; no single global projective embedding |
| Finite-field points | Complete field-valued enumeration | Candidate budgets bound the ambient enumeration |
| Parameter families and morphisms | Univariate, hypersurface, affine-linear, and coupled nonlinear parameter strata; generic degree of affine maps | Exact constructible dimension and bounded geometric-point strata over QQ; specialization-aware Gröbner bases on supported parameter branches, including recursively split algebraic exceptional graphs; generic degrees for characteristic-zero affine-space targets; finiteness and relative-differential branch loci for finite affine morphisms |
| Function decomposition | Rational components, implicit reconstructions and chains | Declared characteristic-zero constants and documented completeness boundaries |
| Real verification | Domains, poles, selectors and graph conditions | Optional semialg extra; algebraic equality alone does not preserve real branches |

Choose the declared coefficient field before computing: a result over the rationals ℚ, written `QQ` in the API, need not remain valid after specialization or a change of field. Resource exhaustion raises or returns the incomplete search record documented by the interface. A bounded rational-point search cannot prove that a variety has no rational points.

The [capability matrix and mathematical scope](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/capabilities.md), [coefficient-field matrix](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/coefficient_fields_support.md) and [input limitations](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/known_issues.md) describe supported inputs, completeness boundaries and explicit refusals.

## Gröbner bases and exact verification

Ideals and coordinate rings delegate Gröbner-basis construction to algroots, using `groebner_backend="native"` by default. The native prime-field F4 kernel is `f4_tiled_finish`; prime-field calculations preserve the characteristic, and rational calculations use rational reconstruction. Accelerated bases pass independent exact verification before use in normal forms. Public timing includes that verification.

Use `ring.ideal(..., groebner_backend="sympy")` or `ring.quotient(..., groebner_backend="sympy")` for portable construction. Higher-level workflows such as function decomposition can also invoke native algroots operations and require the compiled kernels.

For a portable ideal calculation:

```python
import sympy as sp
from alggeom.ideals import PolynomialRing

x, y = sp.symbols("x y")
ring = PolynomialRing((x, y))
ideal = ring.ideal(x**2 - y, y**2 - 1, groebner_backend="sympy")
assert ideal.normal_form(x**4 - 1) == 0
```

Use `ideal.groebner_basis_from("native", limits=...)` for an explicitly requested calculation, or pass `groebner_limits` to the ideal or quotient factory. Availability errors, verification failures and exhausted budgets are reported to the caller. See [quotients and Gröbner backends](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/quotient_localization_and_groebner_backend.md) and the [algroots Gröbner guide](https://github.com/BhuvaneshBhatt/algroots/blob/main/docs/groebner-bases.md) for field support and limit semantics.

## Applications

| Problem | Executable application |
| --- | --- |
| Normalize a cusp, take a finite fiber and select its upper real branch | [Cusp and real branch](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/cusp_branch.py) |
| Separate quadratic, linear, empty and entire-line fibers | [Quadratic fibers](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/quadratic_fibers.py) |
| Count points over prime and extension fields | [Finite-field points](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/finite_points.py) |
| Derive the reachable disk of a planar two-link arm | [Linkage workspace](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/linkage_workspace.py) |
| Compare finite projection with integral closure | [Affine normalization](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/affine_normalization.py) |
| Compute kernels, Hom and tensor of torsion modules | [Module arithmetic](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/module_arithmetic.py) |

The [application tutorials](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/worked_applications.md) state the expected answers and assumptions. The [example catalog](https://github.com/BhuvaneshBhatt/alggeom/blob/main/examples/README.md) covers the complete set of runnable scripts.

## Documentation

Start with [getting started](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/quickstart.md), then choose a mathematical task in the [documentation index](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/index.md). The [API reference](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/api.md) gives signatures and result types.

| Topic | Guides |
| --- | --- |
| Geometric objects and conventions | [Geometry](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/geometry.md), [schemes](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/schemes.md), [conventions](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/conventions.md) |
| Exact images and parameter fibers | [Constructible images](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/constructive_chevalley.md), [nonlinear families](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/nonlinear_families.md), [coupled fibers](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/comprehensive_fibers.md) |
| Polynomial bases and parameter arithmetic | [Gröbner covers](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/groebner_covers.md), [prime parameter ideals](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/parameter_ideals.md), [quotient localizations](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/quotient_localization_and_groebner_backend.md) |
| Algebraic extensions and geometric components | [Constant fields](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/function_field_constants.md), [absolute factorization](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/absolute_factorization.md), [geometric constructions](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/geometric_constructions.md) |
| Function decomposition and real domains | [Decomposition guide](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/guide.md), [worked applications](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/worked_applications.md) |

The [mathematical references](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/references.md) and [bibliography](https://github.com/BhuvaneshBhatt/alggeom/blob/main/REFERENCES.md) explain the algorithms and conventions behind these interfaces.

## Development and testing

From a source checkout with matching backends installed:

```sh
python -m pip install ".[test,real]"
python -m pytest
ruff check .
ruff format --check .
python tools/generate_api.py --check
```

Use `python tools/install_backends.py --profile real` to install matching sibling backend checkouts when available, or add `--published` to test the declared released dependencies. For the core profile, install `".[test]"` and run `python -m pytest --without-real`.

The [testing guide](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/testing.md) covers installed-wheel examples, source-distribution rebuilding, optional dependency isolation and certificate checks. The [benchmark guide](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/benchmarks.md) separates performance measurements from correctness tests. Use the [release checks](https://github.com/BhuvaneshBhatt/alggeom/blob/main/docs/releasing.md) before publication.

## License

alggeom is licensed under the GNU General Public License, version 3 only (`GPL-3.0-only`). See [LICENSE](https://github.com/BhuvaneshBhatt/alggeom/blob/main/LICENSE).
