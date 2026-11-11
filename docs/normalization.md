# Finite projections and affine normalization

Two constructions called normalization solve different problems. Noether normalization finds a polynomial subring over which the given coordinate ring is finite. Integral closure adds missing integral functions in the same total function field, producing a normal ring and a finite birational morphism. Neither construction means discarding nilpotents without an explicit request.

## Noether normalization

```python
import sympy as sp
from alggeom.ideals import PolynomialRing
from alggeom.noether import noether_normalization

x, y, z = sp.symbols("x y z")
A = PolynomialRing((x, y)).quotient(x * y - 1)
projection = noether_normalization(A)
assert projection.verify()
assert projection.base.krull_dimension == 1
assert projection.generic_rank == 2
assert projection.base_embedding.is_injective
assert projection.inverse_map.compose(projection.coordinate_map).images == A.gens
```

The result records an isomorphic transformed coordinate ring, coordinate maps in both directions, an embedded polynomial base, triangular monic integral relations and module spanning generators. The geometric `.projection` is a finite morphism to affine space. Merely finding algebraically independent coordinates would establish generic finiteness, which is weaker: the original x-coordinate on xy=1 is not a finite projection to the affine line.

The polynomial base uses the final coordinates of the transformed presentation. Verification binds its embedding to those coordinate classes and to the transformed ring. Changing the embedding to a power of a coordinate changes the projection and its generic rank, so injectivity alone is insufficient to verify the stored projection.

The finite module generators are a spanning set, not a claim of freeness. In a non-equidimensional or nonreduced ring, generic rank can be smaller than their number. Successive monic elimination gives constructive finite coordinates. Over small finite fields, a bounded linear search falls back to triangular polynomial changes; no field-size hypothesis is imposed. `max_transform_degree` and `max_generators` bound the resulting presentation. The zero ring is excluded. See the [Noether normalization conventions](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/NoetherNormalization/html/index.html).

## Integral closure

```python
from alggeom.normalization import integral_closure, normalize
from alggeom.schemes import AffineScheme

S = PolynomialRing((x, y))
cusp = S.quotient(y**2 - x**3)
closure = integral_closure(cusp)
assert closure.verify() and closure.mapping.is_injective
assert closure.ring.is_normal and not cusp.is_normal
assert closure.conductor == S.ideal(x, y)
assert normalize(cusp).is_normal
map_to_cusp = AffineScheme(cusp).normalization_map()
assert map_to_cusp.target.coordinate_ring == cusp
assert map_to_cusp.source.is_normal
```

`normalize(A)` returns the mathematical coordinate ring directly. `normalization_map(A)` returns its ring homomorphism A→Abar; for an affine scheme it returns the contravariant morphism Spec(Abar)→Spec(A). `integral_closure(A)` additionally returns replayable evidence. For domains, `.inverse_functions` express normalized coordinates in the original function field, `.module_generators` span Abar over A, and `.conductor` is the ambient inverse-image ideal of {a in A : a Abar ⊆ A}.

The Grauert–Remmert algorithm uses a radical Jacobian test ideal J and nonzero g in J. It computes End_A(J)=(gJ:J)/g, adds its module generators as integral coordinates, and repeats until the endomorphism ring equals A. Each extension is presented by its saturated graph. `verify()` replays the ideal operations, finite module evidence, rational inverse coordinates and final normality test. It may cost as much as construction. See the [normalization algorithm](https://agag-mschulze.math.rptu.de/download/nor.pdf).

Step verification also requires the same coefficient field, the retained source coordinates, and a regular denominator: I:g must equal I. A nonzero zero divisor is insufficient, since saturating its graph can delete a component instead of producing a birational extension.

```python
surface = PolynomialRing((x, y, z)).quotient(x**2 - y**2 * z)
normal_surface = integral_closure(surface)
assert normal_surface.verify() and normal_surface.ring.krull_dimension == 2
assert normal_surface.conductor == surface.ambient.ideal(x, y)

crossing = S.quotient(x * y)
separated = integral_closure(crossing)
assert len(separated.components) == 2 and separated.verify()
assert separated.mapping.is_injective
```

Reduced reducible rings are normalized componentwise and assembled as a product with explicit orthogonal idempotents. Components that intersect show why a reduced ring may still fail to be normal. Disjoint normal components form a normal ring. The zero ring is normal by the empty-scheme convention. Nonreduced input raises; call `.radical()` and form the reduced quotient explicitly when that is the desired operation. The product conductor is not currently exposed.

General affine integral closure supports the declared characteristic-zero exact fields and bounded ideal computations, including positive-genus curves, space curves and higher-dimensional models. It does not depend on rational parametrization. `max_steps`, `max_new_variables`, `max_module_generators` and `max_components` raise on exhaustion; an intermediate extension is never returned as a completed normalization.

`is_normal(ProjectiveScheme)` checks affine charts, but assembling a normalized projective model and gluing its maps are not implemented. The [plane-curve interface](curves.md) separately computes smooth projective genus and rationality certificates and gives P^1 normalization maps for rational models. Its `affine_normalization()` method connects to the general affine construction. Finite-characteristic integral closure and automatic preservation of real branch conditions remain unsupported.

## Normalization of integral projective curves

`projective_normalization_atlas(curve)` normalizes the standard projective affine charts by integral closure. `atlas.transition(i, j)` returns `NormalizationTransition(guard, coordinates)` for the map from the normalized source patch `i` to normalized target patch `j`, defined on `D(guard)`. Every coordinate is represented in the localization of the source normal ring: an exact ideal-membership calculation finds a numerator over a power of the overlap guard. The transition is checked against the target equations and the original projective coordinate change.

`atlas.verify_cocycle()` checks identity and triple-overlap composition in the exact fraction fields. This supplies explicit affine gluing for integral projective curves whenever all required chart integral closures complete. It does not produce one homogeneous projective embedding. The normalization is therefore an abstract projective scheme presented by glued affine normal charts.
