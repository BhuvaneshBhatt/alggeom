"""Prepared workloads; setup and certificate checks are outside measured intervals."""

from collections.abc import Callable
from dataclasses import dataclass

import sympy as sp
from algroots.quotient import QuotientAlgebra

from alggeom.components import primary_decomposition
from alggeom.decomposition.algebraic import reconstruct_algebraic
from alggeom.decomposition.chains import rational_chain
from alggeom.decomposition.factors import discover_factors
from alggeom.decomposition.rational import search_rational
from alggeom.fields import finite_field
from alggeom.ideals import PolynomialRing
from alggeom.local import tangent_cone
from alggeom.noether import noether_normalization
from alggeom.normalization import integral_closure
from alggeom.resolutions import free_resolution
from alggeom.schemes import AffineScheme

CASES = (
    "rational_search",
    "number_field_search",
    "factor_discovery",
    "reconstruction",
    "certificate_replay",
    "rational_chain",
    "real_branches",
    "primary_components",
    "absolute_monomial",
    "absolute_graph",
    "groebner_generic",
    "groebner_nested",
    "graded_resolution",
    "local_tangent_cone",
    "noether_projection",
    "integral_closure",
    "finite_points",
    "normalization_family",
    "finite_product_family",
)
REAL_CASES = frozenset({"real_branches"})


@dataclass(frozen=True)
class Workload:
    """A timed operation, an independent result check and structural input measurements."""

    operation: Callable[[], object]
    validate: Callable[[object], bool]
    dimensions: dict[str, int]


def prepare(name, suite):
    """Build one small or standard workload without importing optional real code in core cases."""
    x, z, a, t = sp.symbols("x z a t")
    degree = 8 if suite == "smoke" else 12
    if name == "local_tangent_cone":
        y = sp.Symbol("y")
        ring = PolynomialRing((x, y))
        ideal = ring.ideal(x * x + y**3, x * x - y**3)
        return Workload(
            lambda: tangent_cone(ideal, (0, 0)),
            lambda result: result == ring.ideal(x * x, y**3),
            {"coordinates": 2, "relations": 2, "local_length": 6},
        )
    if name == "noether_projection":
        y = sp.Symbol("y")
        ideal = PolynomialRing((x, y)).ideal(x * y - 1)
        return Workload(
            lambda: noether_normalization(ideal),
            lambda result: result.verify() and result.generic_rank == 2,
            {"coordinates": 2, "dimension": 1, "generic_rank": 2},
        )
    if name == "integral_closure":
        y = sp.Symbol("y")
        ring = PolynomialRing((x, y))
        source = ring.quotient(y * y - x**3)
        return Workload(
            lambda: integral_closure(source),
            lambda result: result.verify() and result.conductor == ring.ideal(x, y),
            {"coordinates": 2, "dimension": 1, "conductor_codimension": 2},
        )
    if name == "finite_points":
        y = sp.Symbol("y")
        field = finite_field(3, 2)
        curve = AffineScheme(PolynomialRing((x, y), field).ideal(y - x * x))
        return Workload(
            curve.rational_points,
            lambda result: len(result) == 9 and len({p.coordinates for p in result}) == 9,
            {"coordinates": 2, "field_order": 9, "ambient_candidates": 81},
        )
    if name == "normalization_family":
        y = sp.Symbol("y")
        n = 2 if suite == "smoke" else 3
        ring = PolynomialRing((x, y))
        source = ring.quotient(y * y - x ** (2 * n + 1))
        return Workload(
            lambda: integral_closure(source),
            lambda value: value.verify() and value.conductor == ring.ideal(x**n, y),
            {"coordinates": 2, "semigroup_generator": 2 * n + 1, "conductor_exponent": 2 * n},
        )
    if name == "finite_product_family":
        from algroots.finite_decomposition import (
            primary_decomposition as finite_primary_decomposition,
        )

        y = sp.Symbol("y")
        m, n = (2, 3) if suite == "smoke" else (4, 5)
        algebra = QuotientAlgebra.from_polynomials(
            (x**m * (x - 1) ** n, y * y), (x, y), domain=finite_field(3, 2)
        )
        return Workload(
            lambda: finite_primary_decomposition(algebra),
            lambda value: (
                sorted(c.dimension for c, k in value) == [2 * m, 2 * n]
                and all(k.dimension == 1 and k.is_field for c, k in value)
            ),
            {
                "coordinates": 2,
                "field_order": 9,
                "algebra_dimension": 2 * (m + n),
                "primary_factors": 2,
            },
        )
    if name == "absolute_monomial":
        from alggeom.geometric import absolute_primary_decomposition

        y = sp.Symbol("y")
        ideal = PolynomialRing((x, y)).ideal(x**2, x * y)
        return Workload(
            lambda: absolute_primary_decomposition(ideal),
            lambda result: (
                len(result.components) == 2
                and ideal == result.components[0].intersection(result.components[1])
            ),
            {"coordinates": 2, "generators": 2, "embedded_primes": 1},
        )
    if name == "absolute_graph":
        from alggeom.geometric import absolute_primary_decomposition

        y = sp.Symbol("y")
        ideal = PolynomialRing((x, y, z)).ideal(z - x * y)
        return Workload(
            lambda: absolute_primary_decomposition(ideal),
            lambda result: len(result.components) == 1 and result.components[0] == ideal,
            {"coordinates": 3, "generators": 1, "dimension": 2},
        )
    if name == "groebner_generic":
        from alggeom.groebner_cover import comprehensive_groebner_cover

        y = sp.Symbol("y")
        return Workload(
            lambda: comprehensive_groebner_cover((a * x, x + y), (x, y), (a,)),
            lambda result: len(result.pieces) == 2 and result.verify_partition(),
            {"fiber_variables": 2, "parameters": 1, "parameter_branches": 2},
        )
    if name == "groebner_nested":
        from alggeom.groebner_cover import comprehensive_groebner_cover

        b, y = sp.symbols("b y")
        return Workload(
            lambda: comprehensive_groebner_cover(((a**2 + b**2) * x, a * y), (x, y), (a, b)),
            lambda result: len(result.pieces) >= 4 and result.verify_partition(),
            {"fiber_variables": 2, "parameters": 2, "nonlinear_exceptional_degree": 2},
        )
    if name == "primary_components":
        y = sp.Symbol("y")
        ring = PolynomialRing((x, y))
        ideal = ring.ideal(x * x, x * y)

        def validate(result):
            return (
                len(result) == 2
                and result[0].intersection(result[1]) == ideal
                and any(q == ring.ideal(x) for q in result)
                and all(
                    q == ring.ideal(x) or (q.radical_contains(x) and q.radical_contains(y))
                    for q in result
                )
            )

        return Workload(
            lambda: primary_decomposition(ideal),
            validate,
            {"coordinates": 2, "expected_components": 2, "embedded_primes": 1},
        )
    if name == "graded_resolution":
        y, w = sp.symbols("y w")
        ideal = PolynomialRing((x, y, z, w)).ideal(x * z - y * y, x * w - y * z, y * w - z * z)

        def validate(result):
            return (
                result.complete and result.shifts == ((0,), (2, 2, 2), (3, 3)) and result.verify()
            )

        return Workload(
            lambda: free_resolution(ideal),
            validate,
            {"coordinates": 4, "relations": 3, "resolution_length": 2},
        )
    if name == "rational_search":
        sample = search_rational(x**degree, x, parameter=t)

        def validate(result):
            return (
                result.complete
                and result.verify()
                and len(result.decompositions) == (2 if degree == 8 else 4)
            )

        return Workload(
            lambda: search_rational(x**degree, x, parameter=t),
            validate,
            {"input_degree": degree, "candidate_count": sample.candidate_count},
        )
    if name == "number_field_search":
        constants = sp.QQ.algebraic_field(sp.sqrt(2))

        def validate(result):
            return result.complete and result.verify() and len(result.decompositions) == 3

        return Workload(
            lambda: search_rational(x**2 + 2 / x**2, x, parameter=t, constants=constants),
            validate,
            {"input_degree": 4, "constant_field_degree": 2, "expected_classes": 3},
        )
    if name == "factor_discovery":
        relation = (
            x**3 - 3 * x - z
            if suite == "smoke"
            else x**4 - z * x**2 * (x + 1) + z**3 * (x + 1) ** 2
        )
        dimension = 3 if suite == "smoke" else 6

        def validate(result):
            return (
                result.complete
                and result.verify()
                and sum(o.extension_degree for o in result.orbits) == dimension
            )

        return Workload(
            lambda: discover_factors(relation, (x, z), degrees=(2,)),
            validate,
            {
                "output_degree": 3 if suite == "smoke" else 4,
                "factor_degree": 2,
                "factor_algebra_dimension": dimension,
            },
        )
    if name == "reconstruction":
        relation = x**4 - z * x**2 * (x + 1) + z**3 * (x + 1) ** 2
        ambient = QuotientAlgebra.from_polynomials(
            [a**2 - z * a + z**3], (a,), domain=sp.QQ.frac_field(z)
        )

        def validate(result):
            return (
                result.verify()
                and result.decomposition.exponent == 1
                and not result.coefficient_field.contains(z)
            )

        return Workload(
            lambda: reconstruct_algebraic(
                relation, (x, z), candidate=x**2 - a * (x + 1), ambient=ambient, parameter=t
            ),
            validate,
            {"output_degree": 4, "extension_degree": 2, "candidate_degree": 2},
        )
    if name == "certificate_replay":
        discovery = discover_factors(x**3 - 3 * x - z, (x, z), degrees=(2,))
        return Workload(
            discovery.verify,
            lambda result: result is True,
            {"output_degree": 3, "factor_degree": 2, "factor_algebra_dimension": 3},
        )
    if name == "rational_chain":
        degree = 8 if suite == "smoke" else 16

        def validate(result):
            return (
                result.complete
                and result.verify()
                and result.compose() == x**degree
                and len(result.components) == (3 if degree == 8 else 4)
            )

        return Workload(lambda: rational_chain(x**degree, x), validate, {"input_degree": degree})
    if name == "real_branches":
        from semialg.function_verification import verify_real_function

        x, y = sp.symbols("x y", real=True)
        graph = sp.And(sp.Eq(y**2, x), y >= 0)

        def validate(result):
            return result.total and result.single_valued and result.certified

        return Workload(
            lambda: verify_real_function(graph, input=x, output=y, domain=x >= 0),
            validate,
            {"variables": 2, "graph_degree": 2, "branch_inequalities": 1},
        )
    raise ValueError(f"unknown benchmark: {name}")
