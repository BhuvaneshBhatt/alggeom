# Exact affine and projective geometry

Use these interfaces when a parametrization needs implicit equations, unwanted components need to be removed from an elimination, or the dimension and singularities of a polynomial model matter. The same ideals also define finite module presentations, so Jacobian geometry and differential calculations share one coefficient field and coordinate convention.

The examples work over QQ (the rationals ℚ). Replace `domain` with `QQ.algebraic_field(alpha)` or a rational function field `K.frac_field(t)` for exact algebraic constants or generic parameters. A coefficient parameter is a scalar, not a geometric coordinate; its exceptional specializations require a separate calculation.

## Ideals and coordinate rings

`PolynomialRing` stores ordered coordinates, an exact coefficient field, and a global term order (`lex`, `grlex`, or `grevlex`). `Ideal` stores generators, computes a reduced Groebner basis, and supplies membership witnesses in the original generators. Ideal equality is mathematical equality in the same ambient ring, independent of term order. A zero generator is omitted and the empty generator list defines the zero ideal.

```python
import sympy as sp
from alggeom.ideals import CoordinateRing, PolynomialRing

x, y, z, t, h, n = sp.symbols("x y z t h n")
S = PolynomialRing((x, y))
I = S.ideal(x * y, x**2)
witness = I.membership(x**2 + y * x * y)
assert sp.expand(sum(c * g for c, g in zip(witness, I.generators)) - (x**2 + x * y**2)) == 0
assert I.normal_form(x**2 + y) == y
assert I.quotient(x) == S.ideal(x, y)
assert I.intersection(S.ideal(y)) == S.ideal(x * y)
assert I.saturate(x) == S.ideal(1)
```

`I.quotient(J)` means `I:J = {f : fJ ⊆ I}` and `I.saturate(J)` means `I:J^∞ = ⋃ₘ I:Jᵐ`. For a multi-generator ideal, saturation is the intersection of the saturations by its individual generators. Saturating by their product generally gives a different answer. For example, `(xy):(x,y)^∞=(xy)`, whereas `(xy):(xy)^∞=(1)`.

`CoordinateRing(I)` represents `S/I` without assuming that it is a field, reduced, irreducable, or zero-dimensional. Calling the ring constructs a parent-aware element, and `normal_form` returns a reduced polynomial expression. Elements support exact arithmetic and inversion of units; homomorphisms validate equations and support kernels, contractions, lifting and inversion. The [scheme tutorial](schemes.md) develops these operations and their contravariant geometric maps. General quotient-ring arithmetic stays in alggeom; the `finite_algebra()` bridge is available only when the quotient is finite over its coefficient field.

```python
A = CoordinateRing(S.ideal(y**2 - x**3))
assert A.krull_dimension == 1
assert A.vector_space_dimension is sp.oo
assert A.normal_form(x**3) == y**2

B = S.quotient(x**2, y**3)
assert B.krull_dimension == 0
assert B.vector_space_dimension == 6
assert B.finite_algebra().dimension == 6
```

Krull dimension measures chains of prime ideals. Vector-space dimension counts a basis over the coefficient field; it is infinite for a nonempty positive-dimensional affine algebra. `algroots.QuotientAlgebra.dimension` always has the latter, finite meaning. The zero ring has Krull dimension `-1` and vector-space dimension zero here. A polynomial ring with no coordinates is its coefficient field, with Krull dimension zero and vector-space dimension one.

## Graphs, implicitization and poles

A rational map has a source coordinate ring, rational coordinate expressions, and fresh target symbols. Graph construction clears denominators and then saturates by their product. Eliminating source coordinates gives the scheme-theoretic image closure. This is useful for removing a redundant parametrization and recovering the equations of an algebraic model.

```python
from alggeom.maps import RationalMap, implicitize, tangent_space, singular_locus

image = implicitize((t**2, t**3), t, targets=(x, y))
assert image == S.ideal(y**2 - x**3)
assert tangent_space(image, (0, 0)) == sp.eye(2)
assert singular_locus(image) == S.ideal(x**2, y)

reciprocal = RationalMap(PolynomialRing((t,)), (1 / t,), (x,))
assert reciprocal.graph_ideal == reciprocal.graph_ideal.ring.ideal(t * x - 1)
assert reciprocal.image_closure.is_zero
assert reciprocal.domain_guard == sp.Ne(t, 0)
```

The reciprocal map misses zero, although its image closure is the entire affine line. `domain_guard` records input exclusions; an image-closure ideal does not encode missing image points. Visible denominators are retained before cancellation, and map composition carries exclusions from both stages. Exclusions already simplified away before passing an expression cannot be recovered. For exact real images, inequalities, or selected branches, use `semialg`.

Tangent directions are columns in an immutable matrix over the declared coefficent field. Supply an exact point on the input scheme; extend the coefficient field explicitly for points outside it. `tangent_ideal` gives the translated affine tangent space. These operations work for nonreduced and mixed-dimensional ideals.

`singular_locus` returns the Jacobian singular subscheme for equidimensional inputs, including scheme multiplicity. It infers equidimensionality from minimal primes when necessary; `assume_equidimensional=True` skips that check when the hypothesis is already known. For a mixed-dimensional input, request `reduced=True` to combine dimension strata and return the reduced singular set. The default cusp singular ideal above is `(x²,y)`; `reduced=True` gives `(x,y)`. Scheme inputs return schemes of the same affine or projective kind. See [local geometry](local_geometry.md) for smoothness, tangent cones and multiplicities.

## Hilbert invariants and projective closure

Hilbert functions and series count graded pieces of `S/I`, assigning degree one to each coordinate. They require a homogeneous ideal and support the standard grading. The Hilbert polynomial is the eventual polynomial of graded-piece dimensions; it is not the cumulative Hilbert-Samuel polynomial.

```python
from alggeom.graded import (
    degree,
    hilbert_function,
    hilbert_polynomial,
    hilbert_series,
    projective_closure,
    projective_dimension,
)

affine = PolynomialRing((x, y, z)).ideal(y - x**2, z - x**3)
projective = projective_closure(affine, homogenizer=h)
assert projective == projective.ring.ideal(y * h - x**2, z * h - x * y, x * z - y**2)
assert projective.krull_dimension == 2
assert projective_dimension(projective) == 1
assert degree(projective) == 3
assert hilbert_polynomial(projective, parameter=n) == sp.Poly(3 * n + 1, n, domain=sp.QQ)
assert [hilbert_function(projective, d) for d in range(4)] == [1, 4, 7, 10]
assert sp.cancel(hilbert_series(projective, parameter=t) - (1 + 2 * t) / (1 - t) ** 2) == 0
```

Homogenising only the supplied affine generators can introduce components at infinity. `projective_closure` also saturates at the homogenising coordinate. The homogeneous coordinate ring of a projective curve has Krull dimension two; its projective dimension is one. An Artinian graded quotient has eventual Hilbert polynomial zero, empty Proj (dimension `-1`), and `degree` equal to its finite length. A unit ideal has degree zero.

Hilbert numerators are computed from the initial monomial ideal using least-common-multiple inclusion-exclusion. This can grow rapidly with many minimal monomial generators. Weighted gradings, multigraded series and numerical approximation are outside this interface.

## Presentations, syzygies and differentials

A matrix acts on column vectors. `syzygies(M, A)` returns column generators of its kernel over the polynomial or coordinate ring A. It need not return a minimal generating set. `Submodule` provides reduction, exact membership and coefficients in the original columns. `PresentedModule(A, M)` represents the cokernel `A^s → A^r → C → 0`, where M has r rows and s columns.

```python
from alggeom.modules import PresentedModule, Submodule, differential, differentials, syzygies

M = sp.Matrix([[x, y]])
kernel = syzygies(M, S)
assert M * kernel == sp.zeros(1, kernel.cols)
assert Submodule(S, kernel).contains([-y, x])

dual_numbers = PolynomialRing((x,)).quotient(x**2)
nilpotent_kernel = syzygies([[x]], dual_numbers)
assert Submodule(dual_numbers, nilpotent_kernel).contains([x])
assert syzygies([[x]], PolynomialRing((x,))).shape == (1, 0)
```

Relations in the base ring matter: multiplication by x is injective over `QQ[x]` but has kernel `(x)` over `QQ[x]/(x²)`. The implementation lifts module problems over `S/I` to S and adjoins each defining equation times each free basis vector. Projecting the lifted syzygies recovers the coordinate-ring kernel.

Fitting ideals detect where a module needs additional generators. `Fitt_i(C)` is generated by the `(r-i)`-minors of its presentation. Over `A=S/I`, `fitting_ideal` returns an ordinary ideal in S containing I, the inverse image of the ideal in A. This convention lets it be compared directly with geometric ideals.

```python
rank_model = PresentedModule(S, sp.diag(x, y))
assert rank_model.fitting_ideal(0) == S.ideal(x * y)
assert rank_model.fitting_ideal(1) == S.ideal(x, y)
assert rank_model.fitting_ideal(2).is_unit

omega = differentials(A)
assert omega.presentation == sp.Matrix([-3 * x**2, 2 * y])
assert omega.fitting_ideal(1) == singular_locus(A.ideal)
assert omega.normal_form(differential(y**2, A) - differential(x**3, A)) == sp.zeros(2, 1)
```

`differentials(A)` presents the relative Kahler module `Ω_(A/K)` using the transpose Jacobian. Coefficient-field parameters are constants in this differential. The columns returned by `differential(f,A)` represent df modulo the module relations. Presentations are not assumed minimal. The [resolution interface](resolutions.md) computes minimal standard-graded resolutions and Betti tables over polynomial rings. Ext and Tor are not computed.

## Schemes, components and curves

The [scheme tutorial](schemes.md) connects coordinate-ring elements and homomorphisms to affine/projective schemes, exact extension-field points, projective charts and scheme-theoretic fibers. The [component guide](components.md) separates radicals, minimal primes, associated primes and primary components, including embedded multiplicity. The [curve guide](curves.md) exposes plane-curve genus, normalization analysis, rational parametrization and projective normalization maps, with the supported cases stated explicitly.

## Mathematical references

The elimination constructions for intersections, quotients, saturation and projective closure follow Cox, Little and O'Shea, *Ideals, Varieties, and Algorithms*, fourth edition, Chapters 4, 8 and 9. Hilbert invariants, presentations, differentials and Fitting ideals use the commutative-algebra conventions described by Eisenbud, *Commutative Algebra with a View Toward Algebraic Geometry*. Bosch and Hartshorne supply the affine/projective geometric context. Basu, Pollack and Roy's *Algorithms in Real Algebraic Geometry* concerns the real-locus layer provided by semialg; no inequality solver is added to these algebraic ideals.
