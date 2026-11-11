"""The cusp normalizes to P^1; a genus-zero conic can need larger constants."""

import sympy as sp

from alggeom.curves import PlaneCurve

x, y, t, u, v = sp.symbols("x y t u v")
cusp = PlaneCurve(y * y - x**3, (x, y))
normalization = cusp.normalization()
assert (normalization.arithmetic_genus, normalization.delta_sum, normalization.genus) == (1, 1, 0)
assert normalization.verify()
mapping = cusp.normalization_map(parameters=(u, v))
assert mapping.coordinates == (u * u * v, u**3, v**3)
assert mapping.image_closure.ideal == cusp.projective_scheme.ideal
conic = PlaneCurve(x * x + y * y - 3, (x, y))
assert conic.genus == 0
assert conic.rationality(t).status == "no_ground_field_point"
extended = conic.parametrize(t, extend_constants=True)
assert extended.source.domain == sp.QQ.algebraic_field(sp.sqrt(3))
assert sp.cancel(sum(f * f for f in extended.functions) - 3) == 0
print("Cusp normalization genus:", normalization.genus, "map:", mapping.coordinates)
print("Conic parameter field:", extended.source.domain, "coordinates:", extended.functions)
