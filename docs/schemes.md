# Coordinate rings, schemes and fibers

Use coordinate-ring elements when formulas should be equal modulo equations, and use morphisms when equations must be respected by a map. Affine and projective schemes retain multiplicity, so a double point differs from its reduced point even though both have the same field-valued points.

## Arithmetic and contravariance

```python
import sympy as sp
from alggeom.ideals import PolynomialRing

x, y, t = sp.symbols("x y t")
A = PolynomialRing((x, y)).quotient(x * y - 1)
a, b = A.gens
assert a * b == A.one
assert a.inverse() == b
assert (a + b) ** 2 == A(x * x + 2 + y * y)

cusp = PolynomialRing((x, y)).quotient(y * y - x**3)
line = PolynomialRing((t,)).quotient()
pullback = cusp.hom((t * t, t**3), line)
assert pullback.kernel == cusp.ideal
assert pullback.is_injective and not pullback.is_surjective
assert pullback.lift(t * t) == cusp(x)
assert pullback.preimage(line.ambient.ideal(t)) == cusp.ambient.ideal(x, y)
```

A ring element has a parent and a canonical polynomial representative accessible through `as_expr()`. Elements from different parents require an explicit map; they are unhashable because equality is mathematical equality modulo the defining ideal. Inversion succeeds exactly for units, including units in positive-dimensional rings such as the hyperbola above. A nonunit raises `ZeroDivisionError`.

A homomorphism fixes the declared coefficient field and sends source coordinates to target elements. Its constructor checks all source equations. `kernel` and `preimage` are ideals in the source ambient polynomial ring containing its defining ideal. `extend` is the corresponding ambient ideal in the target. `lift` solves polynomial image membership; it raises if the target element has no polynomial preimage. `inverse` requires an isomorphism. Kernels and contractions use graph elimination with fresh symbols, so shared coordinate names do not cause variable capture.

The cusp normalization illustrates contravariance: a geometric map from the line to the cusp induces a ring map from the cusp ring into the line ring. It is injective but omits the function t from its image.

## Points and scheme-theoretic fibers

```python
from alggeom.schemes import AffineScheme

source = AffineScheme(PolynomialRing((t,)))
target = AffineScheme(PolynomialRing((x,)))
square = source.morphism(target, (t * t,))
double_point = square.fiber(target.point((0,)))
assert double_point.coordinate_ring.vector_space_dimension == 2
assert double_point.coordinate_ring(t) != 0
assert double_point.coordinate_ring(t) ** 2 == 0
assert [p.coordinates for p in double_point.rational_points()] == [(0,)]

field = sp.QQ.algebraic_field(sp.sqrt(2))
fiber = square.fiber(target.point((2,), field=field))
assert fiber.domain == field
assert {p.coordinates for p in fiber.rational_points()} == {(-sp.sqrt(2),), (sp.sqrt(2),)}
```

An exact point belongs to its declared scheme and carries a coefficient field containing the scheme's base field. Equations and coordinates are validated in that field. Fibers extend the source coefficients to the point's field, then impose the target-coordinate equations. Fibers retain nilpotents; taking only their points would lose ramification multiplicity.

`rational_points(field=K)` enumerates every distinct K-valued point of a zero-dimensional scheme by its finite-algebra rational univariate representation. Over finite fields, `rational_points` enumerates every field-valued point in any dimension, within an ambient candidate budget. `point_search(height=H)` performs a bounded-height search over QQ in any dimension. Neither operation automatically returns all algebraic points. See [finite fields and point searches](finite_fields.md) for height conventions and search limits. A nonempty scheme can have no rational points: `Spec(QQ[x]/(x*x+1))` is one example. Neither scheme emptiness nor rational-point enumeration is a real-locus decision.

Affine schemes also support closed subschemes, scheme intersections and unions, products, fiber products, and principal open subsets. `principal_open(f)` returns the affine presentation with a new coordinate u and relation `u*f=1`. Base change extends a field embedding; it does not specialize a transcendental coefficient parameter.

## Projective charts and morphisms

```python
from alggeom.schemes import ProjectiveScheme

z = sp.Symbol("z")
conic = ProjectiveScheme(PolynomialRing((x, y, z)).ideal(x * z - y * y))
chart = conic.affine_chart(2)
assert chart.scheme.ideal == chart.scheme.ambient.ideal(x - y * y)
point = chart.scheme.point((4, 2))
assert chart.to_projective(point).coordinates == (1, sp.Rational(1, 2), sp.Rational(1, 4))
assert chart.from_projective(chart.to_projective(point)) == point
assert chart.overlap_guard(conic.affine_chart(0)) == sp.Ne(x, 0)
```

`ProjectiveScheme` requires homogeneous equations and saturates at the irrelevant ideal. Its homogeneous coordinate ring has one more Krull dimension than its nonempty projective scheme. Projective points are normalized by their first nonzero coordinate, so scalar multiples compare equal. `affine_chart(i)` sets the ith projective coordinate to one; chart transitions return rational coordinates together with a separate overlap guard.

Projective morphisms use same-degree homogeneous forms for a projective source, check the target equations, and reject a common base locus. Maps from an affine chart to a projective target use regular coordinates without a common zero. Their image closures and fibers preserve scheme structure. A projective graph ideal describes the product through source and fresh target coordinates; it is saturated in each projective factor, not treated as an ordinary singly graded projective space.

## Scope

These are closed affine and projective schemes over supported exact fields, with explicit polynomial morphisms and projective charts. There is no general scheme gluing, structure-sheaf implementation, arbitrary base-ring scheme, projective bundle, or a general algorithm for decomposing constructible images into finite locally closed pieces. `AffineMorphism.constructible_image` supplies exact geometric-fiber membership and closure; see [geometric constructions](geometric_constructions.md). Rational maps and domain guards remain in `alggeom.maps`; inequalities and real branches remain in semialg. See [components](components.md), [curve models](curves.md), and the [capability matrix](capabilities.md).
