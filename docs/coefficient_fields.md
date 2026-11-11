# Coefficient fields and resultant composition

Suppose an algebraic function is described by an irreducible equation `F(x,z)=0`, with output `x` and input `z`. A decomposition introduces a coordinate `t` and equations `G(x,t)=0` and `H(t,z)=0`. Eliminating `t` must recover the original relation, possibly with a constant multiplier and multiplicity:

\[ \operatorname{Res}_t(G,H)=cF^e,\qquad c\in K^\times,\quad e\ge1. \]

This construction is useful when neither component is an explicit rational function. The coefficient field identifies what the intermediate coordinate needs to encode. The resultant checks that the recovered components describe exactly the requested algebraic relation. It does not choose a real branch; that is a separate graph question.

This tutorial uses only the core dependencies. The complete runnable counterpart is [coefficient_fields.py](../examples/coefficient_fields.py). Read [number fields](number_fields.md) for the effect of changing `K`, and [result semantics](results.md) for the meaning of search completeness.

## 1. Separate constants, the input field, and factor coefficients

There are three fields to keep distinct:

| Field | Role |
| --- | --- |
| `K` | Declared exact constants, usually QQ or a number field |
| `K(z)` | Rational functions of the transcendental input coordinate |
| `L = K(z,a)` | A finite extension in which a useful factor of `F` is defined |

A quotient algebra presents `L` using polynomial equations in `a` over `K(z)`. Its finite dimension is the degree over `K(z)`; it does not say that `z` is algebraic over `K`.

```python
import sympy as sp
from algroots.quotient import QuotientAlgebra
from algroots.subfields import field_normal_form

x, z, a, t = sp.symbols("x z a t")
K = sp.QQ
base = K.frac_field(z)
ambient = QuotientAlgebra.from_polynomials([a**2 - z * a + z**3], (a,), domain=base)
assert ambient.is_field and ambient.dimension == 2
assert ambient.minimal_polynomial(a, t) == sp.Poly(t**2 - z * t + z**3, t, domain=base)
assert ambient.trace(a) == z
assert ambient.norm(a) == z**3
```

Consider the quartic and its quadratic factor in this field:

\[ F=x^4-zx^2(x+1)+z^3(x+1)^2,\qquad g=x^2-a(x+1). \]

The other conjugate is obtained by replacing `a` with `z-a`. Multiplying the two factors gives `F` because `a(z-a)=z³` in the quotient field.

```python
F = x**4 - z * x**2 * (x + 1) + z**3 * (x + 1) ** 2
g = x**2 - a * (x + 1)
conjugate = x**2 - (z - a) * (x + 1)
coefficients = sp.Poly(sp.expand(g * conjugate - F), x).all_coeffs()
assert all(field_normal_form(ambient, coefficient) == 0 for coefficient in coefficients)
```

## 2. Generate a field from the factor coefficients

The monic factor has coefficients `1,-a,-a`. The intermediate field used by reconstruction is `E=K(-a,-a)=K(a)`, generated over the constants. It is not formed by automatically adjoining `z`. In this example `z` is absent from `E`: it satisfies a cubic equation over `K(a)` and cannot be recovered as a rational function of `a`.

Adjoining `z` at this point would replace `E` with `L` and change the question. A decomposition should encode the factor coefficients through one intermediate coordinate, not require that this coordinate rationally recover the original input.

```python
from alggeom.decomposition.algebraic import reconstruct_algebraic

result = reconstruct_algebraic(F, (x, z), candidate=g, ambient=ambient, parameter=t)
E = result.coefficient_field
assert E.generators == (-a, -a)
assert E.contains(a) and not E.contains(z)
assert E.transcendence_degree == 1 and E.genus == 0
assert result.rationality.rational and result.parametrization.verify()
assert result.verify()
```

The backend builds a curve model of `E`, computes the genus of its normalization, and seeks a rational parameter over `K`. Genus zero alone is insufficient: the model also needs a `K`-rational point. A positive genus or a certified absence of a required ground-field point can reject a candidate. Resource limits and unsupported curve operations leave it unresolved. They cannot establish a negative decomposition theorem.

## 3. Recover the two relations

The parameter chosen here represents `-a`. Replacing the coefficients of `g` with their rational expressions in `t` produces the outer equation. The minimal polynomial of `-a` over `K(z)` produces the inner equation:

\[ G=x^2+t(x+1),\qquad H=t^2+zt+z^3. \]

```python
D = result.decomposition
assert D.outer == sp.Poly(x**2 + t * (x + 1), x, t, domain=K)
assert D.inner == sp.Poly(t**2 + z * t + z**3, t, z, domain=K)
assert D.exponent == 1 and D.unit == 1 and D.verify()
eliminated = sp.Poly(sp.resultant(D.outer.as_expr(), D.inner.as_expr(), t), x, z, domain=K)
assert eliminated == sp.Poly(F, x, z, domain=K)
assert not result.complete
```

A supplied factor proves this particular decomposition. It does not establish that all factors or all decompositions were considered. Intermediate generators are not canonical: changing `t` by a fractional-linear transformation can produce different printed components for the same field.

## 4. Discover factor coefficient extensions automatically

When no factor is supplied, `discover_factors` introduces the coefficients of a generic monic factor and sets the division remainder to zero. For a squarefree output-degree `n` polynomial, its universal degree-`d` factor algebra has dimension `binomial(n,d)`. A separating linear form converts this finite algebra into a cyclic presentation; its irreducible components describe Galois orbits of factors over `K(z)`.

```python
from alggeom.decomposition.factors import discover_factors

discovery = discover_factors(F, (x, z), degrees=(2,))
assert discovery.complete and discovery.verify()
(level,) = discovery.levels
assert level.geometric_factor_count == 6
assert [orbit.extension_degree for orbit in level.orbits] == [2, 4]
assert sum(orbit.extension_degree for orbit in level.orbits) == 6
```

There are six geometric quadratic factors, grouped into orbits of size two and four. The extension degree can exceed or fail to divide `n`: it is the degree of the field needed for a whole factor's coefficients, not the degree of one root of `F`. Completeness of this discovery result covers the requested factor degree only. `search_algebraic` additionally reconstructs candidates and accounts for their outcomes; its default can stop at the first successful decomposition.

## 5. Understand resultant powers

A resultant power greater than one is meaningful. Let `F₃=x³-3x-z` and let `a` be one of its roots. The quotient factor `x²+ax+a²-3` contains the other two roots. Its three conjugates contain each root of `F₃` twice, so their product is `F₃²`.

```python
cubic = x**3 - 3 * x - z
cubic_ambient = QuotientAlgebra.from_polynomials([a**3 - 3 * a - z], (a,), domain=base)
quadratic = x**2 + a * x + a**2 - 3
assert all(
    field_normal_form(cubic_ambient, c) == 0
    for c in sp.Poly((x - a) * quadratic - cubic, x).all_coeffs()
)
power_result = reconstruct_algebraic(
    cubic, (x, z), candidate=quadratic, ambient=cubic_ambient, parameter=t
)
assert power_result.verify() and power_result.decomposition.exponent == 2
G = x**2 + t * x + t**2 - 3
H = t**3 - 3 * t - z
assert sp.expand(sp.resultant(G, H, t) - cubic**2) == 0
```

The degree identities provide useful checks: `e degₓ(F) = degₓ(G) degₜ(H)` and `e deg_z(F) = degₜ(G) deg_z(H)`. Here they give `2·3=2·3` and `2·1=2·1`. They do not replace the full resultant identity, which also detects additional components and an incorrect scalar unit.

## 6. Interpret verification correctly

`AlgebraicDecomposition.verify()` checks irreducibility, declared coordinates and constants, nontrivial component degrees, the power, the nonzero constant unit, and the exact resultant. `AlgebraicReconstruction.verify()` also replays the candidate embedding, coefficient field and rational parametrization. Discovery certificates account for every covered orbit, including the inverse coordinates of the separating form. Chain verification additionally checks each split, leaf reason and tree budget.

Verification uses the same exact arithmetic backends as construction; it is not an independently implemented proof checker. An algebraic certificate says nothing about selecting a principal radical, preserving a visible pole, or making a relation single-valued. Use the optional real graph interfaces for those obligations. See [capabilities and limits](capabilities.md) for the supported algebraic search scope and current coefficient-field restrictions.
