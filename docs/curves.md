# Plane curves and normalization

Use `alggeom.curves` to ask whether a singular plane model has a rational parameter, compute the genus of its smooth projective normalization, and connect a parametrization to affine or projective scheme maps. The public interface uses algroots' exact curve algorithms and exposes their certificates for inspecting the underlying mathematical evidence.

```python
import sympy as sp
from alggeom.curves import PlaneCurve, genus, parametrize

x, y, t, u, v = sp.symbols("x y t u v")
curve = PlaneCurve(y * y - x**3, (x, y))
assert genus(curve) == 0
analysis = curve.normalization()
assert (analysis.arithmetic_genus, analysis.delta_sum, analysis.genus) == (1, 1, 0)
assert analysis.verify()
affine_map = parametrize(curve, t)
assert affine_map.functions == (t * t, t**3)
assert affine_map.image_closure == curve.affine_scheme.ideal
projective_map = curve.normalization_map(parameters=(u, v))
assert projective_map.coordinates == (u * u * v, u**3, v**3)
assert projective_map.image_closure.ideal == curve.projective_scheme.ideal
```

The cusp has a degree-three plane equation, arithmetic genus one, and a singularity contribution of one. Its normalization has genus zero. The affine map uses one parameter chart; the projective normalization morphism also includes the point at infinity. A generic rational inverse may fail at a singular point even though the normalization map reaches that point.

## Genus zero and the field of definition

```python
conic = PlaneCurve(x * x + y * y - 3, (x, y))
result = conic.rationality(t)
assert conic.genus == 0
assert result.status == "no_ground_field_point"
assert result.verify()
extended = conic.parametrize(t, extend_constants=True)
assert extended.source.domain == sp.QQ.algebraic_field(sp.sqrt(3))
assert sp.cancel(sum(f * f for f in extended.functions) - 3) == 0
```

Genus zero alone does not guarantee a rational parameter over the declared field. A ground-field point is also needed. Over QQ the backend decides the conic obstruction exactly. Over other number fields it requires a supplied point witness, or an explicit request to extend constants; it does not claim that a failed point search proves nonrationality. For a point witness use `curve.rationality(t, point=(a,b))` or pass the same option to `parametrize`. A chosen point must satisfy the backend's smooth-point requirements.

The specialist rationality result has statuses `rational`, `rational_over_extension`, `positive_genus`, and `no_ground_field_point`. When available, `result.parametrization` stores its parameter, constant field, projective and affine coordinates, generic inverse, guards and exceptional parameter polynomial. It supplies `verify()` and homogeneous-coordinate conversion. Any splitting field and conic model are explicit. `parametrize` returns the mathematical `RationalMap` directly; absence of a requested parametrization raises with the status.

## Supported cases and the meaning of normalization

| Operation | Supported models | Result and boundary |
| --- | --- | --- |
| Plane model | Reduced irreducible affine equations in two coordinates over QQ or a declared number field | Affine scheme and saturated projective closure; no automatic space-curve projection |
| Genus and normalization analysis | Geometrically integral plane curves within exact computation limits | Smooth projective genus, singularity resolution and adjoints; geometric reducibility is rejected |
| Rationality over QQ | Genus-zero models reducible by the backend's adjoint and conic constructions | Exact rational parameter or a conic point obstruction; positive genus is an obstruction |
| Rationality over number fields | Supplied suitable point witness or explicit constant extension | Verified maps over the resulting declared field; no general number-field point decision |
| Projective normalization morphism | Rational models with a computed parametrization | A morphism from P^1 to the projective plane model, including infinity |
| Affine integral closure | Reduced affine characteristic-zero algebras through `normalization` | Explicit normalization rings and finite birational maps, also for positive-genus curves and higher-dimensional models |

`normalization()` returns algroots' `CurveNormalization` analysis certificate: `genus`, `arithmetic_genus`, `delta_sum`, `singularities`, `canonical_adjoints`, and replay verification. The name refers to the smooth projective normalization whose invariants are computed. For an actual affine integral-closure scheme use `curve.affine_normalization()` or `alggeom.normalization`; these are separate from this projective analysis certificate. `normalization_map()` supplies an actual normalized source only in the rational case, when P^1 is available.

Nonreduced and reducible equations must first be separated with the [component interface](components.md). Rational-function constant fields, finite characteristic, general space curves, and automatic preservation of real branch ranges are unsupported by the plane-curve analysis interface. The separate [affine normalization interface](normalization.md) supports general reduced affine characteristic-zero models. `curve.projective_normalization_atlas()` normalizes every nonempty projective affine chart and records its overlap nonvanishing functions, including for positive-genus integral curves. It computes exact regular transitions for integral generators on principal overlaps and verifies their cocycle identities. The affine atlas defines the normalized scheme by gluing, but a single global projective embedding and homogeneous coordinate description are not constructed. See [geometric constructions](geometric_constructions.md). `max_steps` and `max_dimension` are exact backend resource limits; exhaustion raises without a genus or rationality decision. For real loci and branch conditions, retain semialg guards separately.
