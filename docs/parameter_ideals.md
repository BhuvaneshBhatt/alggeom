# Prime parameter quotients and recursive Gröbner specialization

Parameterized polynomial equations often acquire new solutions or different leading terms on exceptional algebraic varieties. A generic Gröbner basis over `QQ(a,b)` does not automatically specialize on the curve `a**2+b**2=0`. That curve has no polynomial parametrization over `QQ`, but its coordinate ring and function field are exact algebraic objects.

The `PrimeParameterField` class computes in `Frac(QQ[a_1,...,a_m]/P)` for a **certified prime** ideal `P`. Residue representatives are reduced by a Gröbner basis of `P`, and field elements are polynomial fractions with nonzero denominator classes. Equality is checked by cross multiplication modulo `P`. A nonprime ideal is rejected: inverting nonzero residues of a ring with zero divisors would not define a field.

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.parameter_fields import PrimeParameterField, groebner_over_prime_parameters

a, b, x, y = sp.symbols("a b x y")
P = PolynomialRing((a, b)).ideal(a**2 + b**2)
field = PrimeParameterField(P)
assert field.element(a**2) == field.element(-(b**2))

result = groebner_over_prime_parameters((a * x + b * y, y**2 - 1), (x, y), P)
assert result.verify()
assert a in result.denominator_guards
```

The Gröbner computation performs exact Buchberger arithmetic over the quotient function field. Its result contains a lexicographic basis, the polynomial nonvanishing guards needed for specialization, and membership identities expressing every basis polynomial in the original generators. The independent verifier checks these identities modulo `P`, reduces the source generators, checks all S-pairs, and rejects missing coefficient-denominator guards. A generic result is valid **only on the open subset** of `V(P)` on which all guards are nonzero. A proof of equality of ideals over `Frac(QQ[a]/P)` is not a proof that the basis specializes on the whole of `V(P)`.

## Recursive covers

`comprehensive_prime_cover(equations, variables, parameters)` starts with `P=(0)` and recursively analyzes the zero loci of inverted coefficients. It decomposes exceptional parameter ideals into field-relative minimal primes, computes a quotient-function-field basis over each prime, and partitions intersections of branches using exact equations and inequations. The output contains a replayable certificate for each piece and a symbolic coverage/disjointness check.

```python
from alggeom.parameter_cover import comprehensive_prime_cover

cover = comprehensive_prime_cover(((a**2 + b**2 + 1) * x, a * y), (x, y), (a, b))
assert cover.verify()
assert cover.specialize((0, sp.I)) == ()
```

This example has four strata, including the intersection `a=0`, `b**2+1=0` of two nonlinear exceptional conditions. The curve `a**2+b**2+1=0` is handled through its implicit quotient field, without solving for a polynomial graph.

## Geometric versus field-relative primality

A prime ideal over `QQ` can become reducible over an algebraic extension. The parameter quotient field only requires primality **over its specified coefficient field**: for example, `(a**2+b**2)` is prime over `QQ` and splits over `QQ(i)`. This is sufficient for exact generic coefficient arithmetic, but is not a geometric-primality certificate. Geometric splitting and primary decomposition remain separate operations. The [function-field constant guide](function_field_constants.md) describes exact sufficient geometric-primality certificates and bounded computation of complete constant fields for supported positive-dimensional prime ideals.

## Exact limits

- Coefficients are rational, parameters and fiber variables are distinct, and lexicographic fiber Gröbner bases are constructed.
- Primality is proved through sufficient exact criteria or the existing associated-prime implementation. Inputs without a supported certificate are rejected.
- Recursive exceptional ideals are decomposed using the current field-relative `minimal_primes` algorithm. If its resource limits or coefficient-domain restrictions are encountered, no partial cover is returned as complete.
- The verifier checks identities using exact arithmetic and constructible-set algebra; it does not establish that arbitrary unsupported positive-dimensional ideals admit automatic *absolute* primary decomposition.
- Branch count, output size, and polynomial-expression growth can be substantial. `max_nodes`, `max_pieces`, and `max_pairs` bound branching and critical pairs, but not wall-clock time.
