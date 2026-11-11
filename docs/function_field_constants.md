# Algebraic constants and geometric primality

A prime ideal over a characteristic-zero field need not remain prime after extending the coefficient field. For a prime ideal \(P\subset \mathbb Q[x_1,\ldots,x_n]\), the relevant invariant is the algebraic closure of \(\mathbb Q\) **inside** the finitely generated function field \(K=\operatorname{Frac}(\mathbb Q[x]/P)\). When that constant field is exactly \(\mathbb Q\), the integral variety is geometrically integral. The issue is different from irreducibility of a polynomial after treating independent variables as rational-function coefficients.

## A nonlinear geometric-primality certificate

`certify_geometric_prime(ideal)` first proves field-relative primality. For a prime ideal with a smooth point rational over its coefficient field, it checks that the Jacobian rank at the point equals its codimension. This certifies geometric primality even for a nonlinear, positive-dimensional, nonprincipal ideal. An integral variety's geometric components form one Galois orbit, while a smooth rational point lies on exactly one geometric component; therefore there can be only one. Affine-linear and polynomial-graph prime ideals supply additional direct certificates. Smooth projective hypersurfaces of bounded degree also supply a sufficient certificate, including some curves without rational points.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.function_field_constants import certify_geometric_prime, algebraic_constant_field

x, y, z = sp.symbols("x y z")
ring = PolynomialRing((x, y, z), sp.QQ, "lex")
curve = ring.ideal(y**2 - x**3 + x, z - x * y)
assert certify_geometric_prime(curve)
assert algebraic_constant_field(curve).field == sp.QQ
```

The function returns `False` when the input is certifiably nonprime; otherwise failure to prove geometric primality raises `NotImplementedError`. A failure to find a rational smooth point does **not** imply that the variety is geometrically reducible.

## Complete algorithm for quadratic function-field presentations

For a prime hypersurface whose defining equation is quadratic in a coordinate `y` over the rational-function field `QQ(t_1,...,t_r)`, write `f=a*y**2+b*y+d` and `Delta=b**2-4*a*d`. Factoring the square class of `Delta` over the polynomial ring in the independent coordinates gives a **terminating and complete algorithm for this input class**:

- A positive-degree factor with odd valuation in `Delta` remains odd after extending the coefficient field to its algebraic closure (characteristic zero). Therefore `f` is absolutely irreducible and the relative algebraic constants are exactly `QQ`.
- Otherwise `Delta=c*s**2` with `c` a nonsquare rational number and `s` a rational function of the independent coordinates. The primitive constant `(2*a*y+b)/s` satisfies `T**2-c`, and because the extension has degree two, its constant field is exactly `QQ(sqrt(c))`. The geometric branches are the saturations of `(f, 2*a*y+b +/- s*sqrt(c))` by `a`.

The result includes a `QuadraticConstantCertificate`; calling `constants.verify(prime)` recomputes the discriminant independently using a separate irreducible-factorization path, checks the defining polynomial identity, reconstructs both saturated graph ideals, proves that each graph remains prime over an algebraic closure, and checks the intersection. A certified `QQ` answer has an independent odd-valuation witness.

```python
x, y = sp.symbols("x y")
prime = PolynomialRing((x, y)).ideal((x**7 + 1) ** 2 * y**2 - 2)
constants = algebraic_constant_field(prime, max_candidates=1)
assert constants.verify(prime)
assert constants.field.mod.degree() == 2
assert len(constants.geometric_primes) == 2
```

The fallback for other finitely generated function fields remains bounded and incomplete, not a general terminating algorithm. Do not infer that the constant field is `QQ` from an inconclusive search.

## Computing a full constant field with certificates

For a rational function \(c=p/q\in K\), alggeom adds a symbol \(T\), saturates the graph ideal \((P,Tq-p)\) by \(q\), and eliminates the coordinates. A nonzero elimination ideal gives the exact minimal polynomial of \(c\) over \(\mathbb Q\). The polynomial identity is checked directly modulo \(P\). The algorithm then adjoins its algebraic roots, constructs the conjugate graph components, certifies each as geometrically prime, and checks that their intersection equals the base-changed ideal. **Only after these checks is the entire constant field returned as computed.** Merely discovering an algebraic subfield is not sufficient to claim completeness.

```python
prime = ring.ideal(x**2 + y**2, z - y**2)
constants = algebraic_constant_field(prime)
assert sp.Poly(constants.minimal_polynomial).all_coeffs() == [1, 0, 1]
assert (constants.numerator, constants.denominator) == (x, y)
assert len(constants.geometric_primes) == 2
```

In this example \(x/y\) satisfies \(T^2+1\) in the function field, so the full constant field is \(\mathbb Q(i)\). A nonnormal example is \((x^3-2y^3,z-y^2)\): the constant field is \(\mathbb Q(\sqrt[3]{2})\), while the larger splitting field also contains the cube roots of unity. `AlgebraicConstants.field` describes the constant field; `splitting_field` describes the common coefficient field of all returned geometric components. Their degrees need not agree.

The associated-prime splitter can use these proven constant fields to extend the **original ideal** and recompute its primary decomposition, preserving embedded components. It accepts the result only after checking geometric primality of every extended associated prime, component primariness, and exact reconstruction.

## Scope and computation limits

- The constant-field search currently requires a positive-dimensional prime ideal over \(\mathbb Q\). Its bounded candidate set contains coordinate rational functions and, optionally, ratios of low-degree monomials. A primitive constant outside that set may not be found.
- `max_degree`, `numerator_degree`, `max_candidates`, and `max_points` bound portions of the search. Gröbner elimination and primary decomposition can still be expensive. Inconclusive cases raise `NotImplementedError`, not a false claim of a trivial constant field.
- The geometric-primality test is a collection of **sufficient, exact criteria**, not an unconditional decision procedure for arbitrary ideals. The returned constant field is complete *when returned*, because all conjugate factors receive geometric-primality certificates.
- Coefficient fields other than \(\mathbb Q\) are accepted by `certify_geometric_prime` when their ideal arithmetic is supported; constant-field discovery is currently over \(\mathbb Q\).

See [components](components.md), [prime parameter quotients](parameter_ideals.md), [capabilities](capabilities.md), and [mathematical references](references.md).
