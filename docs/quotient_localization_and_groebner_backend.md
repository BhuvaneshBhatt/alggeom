# Quotients and Gröbner backends

Use `ideal.groebner_basis_from("native", limits=...)` for an explicitly selected verified basis. It returns polynomial expressions in the declared coefficient field and coordinate order. Prime-field requests compute directly in that field; QQ requests use rational reconstruction. The cached `ideal.groebner_basis` property uses the policy stored on its ideal, defaulting to native. Select `ring.ideal(..., groebner_backend="sympy")` or `ring.quotient(..., groebner_backend="sympy")` for portable construction.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing

x, y = sp.symbols("x y")
ring = PolynomialRing((x, y), sp.GF(101), "grevlex")
ideal = ring.ideal(x**2 - y, y**2 - 1)
assert set(ideal.groebner_basis_from("sympy")) == set(ideal.groebner_basis)
```

The public method always verifies accelerated results and has no candidate-only switch. Independent exact comparison checks equality of the generated ideals, including both inclusions. Verification time is part of total API time. `groebner_basis_proof_from` is the separate expert interface for explicit ideal-membership proof objects, witness budgets, and change-of-order strategies.

`algroots.QuotientAlgebra.from_polynomials(..., groebner_backend="native", groebner_limits=...)` constructs a finite quotient using the same verified policy. Prime-field arithmetic preserves the characteristic; these algebras are not characteristic-zero solution sets. Quotient construction still requires a finite-dimensional quotient.

The compiled kernel must be installed for the default native construction. A pure Python installation requires explicit SymPy selection. Shared `GroebnerLimitError`, availability errors, and verification failures propagate to the caller without an automatic backend change. Backend budgets do not impose a wall-clock limit on independent reference verification.

The ideal and quotient factories accept `groebner_backend` and `groebner_limits` and store the selection on that defining ideal. Specify the selection when constructing a different ideal.
