# Working with compositions

For polynomial ideals, coordinate rings, maps, graded invariants and modules, start with the [geometry tutorial](geometry.md).

## Common right components

Several rational functions can share one intermediate coordinate. `common_right_component` returns a SymPy Lambda generating their joint field when its transcendence degree is one. Degree-one components and linear outer maps are valid here. Constants have no distinguished nonconstant generator; higher-transcendence-degree inputs have no scalar intermediate. Both cases return `None`; the specialist certificate distinguishes them.

```python
import sympy as sp
from alggeom.decomposition import common_right_component
from alggeom.decomposition.common import common_components

u, v, q = sp.symbols("u v q")
shared = u**2 + v**2
functions = (shared**2, shared**3)
inner = common_right_component(functions, (u, v), parameter=q)
certificate = common_components(functions, (u, v), parameter=q)
assert certificate.verify() and inner == certificate.inner
assert all(
    sp.cancel(outer(inner(u, v)) - f) == 0 for outer, f in zip(certificate.outers, functions)
)
assert common_components((u, v), (u, v)).status == "higher_transcendence_degree"
```

This computes one generator of the joint field, with outer reconstructions and visible input poles. It does not enumerate all fields containing that field or every possible shared substitution.

For complete calculations, see the [application workflows](applications.md). For partial searches and stopping reasons, see [result interpretation](results.md).

## Choose the object you mean

A rational expression represents an element of a function field: canceling a common factor does not change that element. A real function also has a domain. The expression `x**4*(x-1)/(x-1)` can represent the same field element as `x**4` while retaining a missing input at `x=1`. To preserve such a hole, provide an unevaluated SymPy expression; information already erased by SymPy cannot be recovered.

An implicit polynomial `F(output,input)=0` describes a relation. It need not give a single output at each input. A real square-root branch additionally needs a sign selector and an input domain. Choose the rational search, implicit reconstruction, or real graph API according to the object under study.

## Rational components and fields

Use `decompose(f, variables)` for one complete chain or `decompositions(f, variables)` for every two-map class. The returned SymPy Lambdas bind their own coordinates, so evaluation does not depend on knowing an internally chosen intermediate symbol. For example, apply the last Lambda to the original inputs and work backwards through the tuple.

A fractional-linear substitution changes the representation of a split: replacing `h` by `(a*h+b)/(c*h+d)` and adjusting `g` gives the same intermediate field when `a*d-b*c != 0`. The search returns one representative per field. Exact representative expressions are not a stable choice of canonical form.

Pass `constants=sp.QQ.algebraic_field(sp.sqrt(2))` to work over a number field. Undeclared constants and floating-point coefficients are rejected. Decomposability over a larger constant field can differ, so an answer over QQ does not settle the question over its algebraic closure.

## Search evidence

`search_rational` returns factors of the cross-difference polynomial, a candidate count, the number examined, and verified decompositions. Its `complete` property means the entire candidate set has been examined. `status` can be `decomposable` even when `complete` is false: finding one split proves existence without exhausting the alternatives.

`rational_chain` chooses the first split at each node. Each leaf stores either a complete negative search or a stopping reason. `complete` concerns those leaves, so it can be true even when a parent search stopped after one successful candidate. The result does not enumerate all possible chains.

The decomposition entry points raise `IncompleteDecompositionError` when their promised complete result is unavailable. Inspect `error.result` for verified partial answers. `DegreeLimitError` is raised by a two-map search before processing an oversized input; chain builders retain this as a `degree_limit` leaf. A valid incomplete chain remains an exact composition identity.

## Implicit algebraic components

`search_algebraic(F, (output,input))` first constructs the algebra of monic factors of each requested output degree. It splits that finite algebra into coefficient fields and tries to reconstruct components from each orbit. `exhaustive=True` continues after successful reconstruction. `factor_degrees` restricts discovery; omitting a relevant degree prevents a complete decomposition search claim.

If an appropriate extension and factor are already known, call `reconstruct_algebraic(F, variables, candidate=g, ambient=A)`. The ambient algebra must be a certified field over `constants.frac_field(input)`. Coefficients of `g` generate the subfield to test. Adjoining the input coordinate to this generated field can erase the structure being sought, so reconstruction keeps the two roles distinct.

The resulting certificate stores both component polynomials, a positive resultant exponent, and a nonzero constant unit. `verify()` checks the full identity and degree constraints. A rejection can mean positive genus, no ground-field point, or a trivial or reducible reconstructed component. Rejection of one supplied candidate is not a global negative result.

## Real domains

`verify_decomposition` accepts one two-component certificate; `verify_chain` accepts a whole chain. Selectors are ordered outermost to innermost and may use only the corresponding component's input and output coordinates. Rational chains let you choose `output` and `intermediates`; implicit components already fix these coordinates.

The returned `equivalent` flag compares selected real graphs on `domain`. With `check_functions=True`, `function_checks` separately test totality and single-valuedness, and `function_verified` requires both graphs to satisfy both properties. Equal multivalued relations can have `equivalent=True` and `function_verified=False`. Empty input domains make these obligations vacuous.

Multivariate rational inputs remain free coordinates during graph comparison. For function checks, the last coordinate is the input and the preceding coordinates are treated as parameters; the domain quantifies over all of them.

## Parameter families

`polynomial_pair` gives the exact condition and normalized components for one ordered degree pair. The inner polynomial is monic with constant zero; this fixes the affine freedom. Its condition includes nonzero leading coefficients, residual equations, and coefficient denominator guards.

`polynomial_family` covers the declared real parameter domain with disjoint cases, including zero and constant polynomials, degree drops, indecomposable specializations, decomposable specializations, and undefined coefficient poles. `select` accepts exact assignments keyed by symbols or unambiguous parameter names. Check the selected case's status before using its components.
