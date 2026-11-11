# Exact arithmetic and solver backends

alggeom delegates Gröbner-basis construction to algroots and uses it for finite quotient algebras, exact rational-univariate representations, finite-field arithmetic, function fields, and plane-curve calculations. These operations supply equality-only algebraic structure. The optional semialg adapter adds real quantifier elimination, inequalities, and branch conditions.

## Installation profiles

The core installation requires Python 3.11 or later and algroots 0.4.0 or newer:

```sh
python -m pip install alggeom
```

Install `"alggeom[real]"` to add semialg 1.5.0 or newer for real decisions and branch verification. Both profiles use SymPy 1.14 within the 1.x series and python-flint 0.9 within the 0.x series. The declared requirements are recorded in [package metadata](../pyproject.toml) and [backend inputs](../ci/backends.toml).

Compatible native algroots wheels include their compiled kernels and the required numerical libraries. If a native wheel is unavailable for the interpreter and platform, follow the [algroots installation guide](https://github.com/BhuvaneshBhatt/algroots/blob/main/README.md) to build the kernels, or select the portable backend for ideal and quotient construction. Native availability is required for the native calculation paths; installing semialg does not supply those kernels.

## Coefficient fields

Public rings retain the declared coefficient field and selected algebraic embedding. The shared `algroots.domains` layer converts complete algebraic coefficients inside that field, including rational functions over number fields. Its internal fraction representation normalizes denominators to control scalar growth. Existing SymPy classes and field instances are not patched.

Module calculations use SymPy's AGCA standard bases. Quotient-ring modules are lifted to polynomial modules with the defining ideal relations adjoined. Saturation uses exact elimination. These operations preserve coefficient characteristic and coordinate order.

`SYMPY_GROUND_TYPES=python` selects SymPy's coefficient representation. Gröbner backend selection is independent: native construction still requires the compiled algroots kernels.

## Gröbner selection and verification

`ring.ideal(..., groebner_backend="native", groebner_limits=...)` and `ring.quotient(...)` select verified accelerated construction, which is also their default. Use `groebner_backend="sympy"` for portable construction. `ideal.groebner_basis_from("native")` computes a separate explicitly requested basis. Prime-field requests run directly in their characteristic; QQ requests use rational reconstruction. The native kernel is fixed to `f4_tiled_finish`; `sympy` is the portable alternative.

Accelerated public results pass independent exact verification before entering geometric normal forms. Candidate-only calculation is an expert algroots kernel interface. Reference verification is included in public timing and can dominate total cost. See [quotients and Gröbner backends](quotient_localization_and_groebner_backend.md) for examples and limits.

The explicit SymPy policy covers ideal and quotient construction. Higher-level workflows, including function decomposition, can invoke additional native algroots operations and require its compiled kernels.

The proof interface `groebner_basis_proof_from` produces explicit ideal-membership evidence. `verify_reduced_groebner_basis` and `GroebnerBasisProof.verify` accept `normal_form_backend="flint"` over QQ for exact rational multivariate division. Both ideal inclusions, witness identities, reducedness, and Buchberger's criterion remain required. A kernel completion flag alone is insufficient.

## Install matching dependencies

[ci/backends.toml](../ci/backends.toml) declares dependency requirements. From an alggeom source checkout, install the released requirements explicitly:

```sh
python tools/install_backends.py --published --profile real
python -m pip install ".[test,real]"
```

For development against sibling source trees, omit `--published`. Use `--source-root /path/to/parent` to choose a combined checkout. The core profile installs algroots without the real adapter; run its tests with `python -m pytest --without-real`. Source installs of algroots need its documented native build configuration when native calculations are required.

CI installs the published requirements and checks the installed alggeom wheel outside its source checkout. It tests both core and real profiles, minimum and latest arithmetic dependencies, and FLINT and Python ground coefficients. Generated API documentation and executable examples validate the adapter contracts. Record backend versions and source fingerprints when comparing timings across dependency trees.
