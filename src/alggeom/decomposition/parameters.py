"""Polynomial decomposition on exact real parameter strata."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

import sympy as sp
from semialg.conditional import (
    ConditionalBranch,
    ParameterStratifiedResult,
    verify_parameter_stratification,
)
from semialg.decision import implies, is_satisfiable
from semialg.simplify.atoms import normalize_atoms

from ._validation import validate_real_symbols


def _parameter_assignments(assignments, parameters):
    resolved = {}
    for key, value in assignments.items():
        if isinstance(key, str):
            matches = [p for p in parameters if p.name == key]
            if len(matches) != 1:
                raise ValueError(f"unknown or ambiguous parameter name {key!r}")
            key = matches[0]
        resolved[key] = value
    return resolved


def _equations(coefficients):
    equations, guards = ([], [])
    for coefficient in coefficients:
        numerator, denominator = sp.fraction(sp.together(sp.cancel(coefficient)))
        equations.append(sp.Eq(numerator, 0))
        if denominator != 1:
            guards.append(sp.Ne(denominator, 0))
    return (tuple(equations), tuple(guards))


@dataclass(frozen=True)
class PolynomialSystem:
    """Normalized coefficient equations for fixed outer and inner degrees.

    ``inner`` is monic with constant zero. ``outer`` has the original leading
    and constant coefficients. This fixes affine changes of the middle variable.
    ``formula`` is existential in ``unknowns`` and free in ``parameters``.
    """

    polynomial: sp.Poly
    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    degrees: tuple[int, int]
    component_variable: sp.Symbol
    outer: sp.Expr
    inner: sp.Expr
    unknowns: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    guard: sp.Expr

    @property
    def formula(self):
        """Exact normalized coefficient relation, including degree and pole guards."""
        return sp.And(self.guard, *self.equations)


@dataclass(frozen=True)
class PolynomialPair:
    """Necessary and sufficient condition and explicit normalized components.

    Components are valid only on ``condition``. The coefficient construction is
    deterministic in characteristic zero; no parameter sampling is used.
    """

    system: PolynomialSystem
    condition: sp.Expr
    outer: sp.Expr
    inner: sp.Expr
    coefficient_values: tuple[tuple[sp.Symbol, sp.Expr], ...]
    method: str = "normalized_coefficient_reconstruction"

    @property
    def composition(self):
        """Candidate composition, to be interpreted under ``condition``."""
        return self.outer.subs(self.system.component_variable, self.inner)


def polynomial_system(
    polynomial,
    variable: sp.Symbol,
    *,
    degrees: tuple[int, int],
    parameters: Sequence[sp.Symbol] = (),
    component_variable: sp.Symbol | None = None,
) -> PolynomialSystem:
    """Build exact coefficient equations for ``polynomial = outer(inner(x))``.

    ``degrees=(outer_degree, inner_degree)`` must contain integers at least two.
    Coefficients may be rational functions of the declared real parameters with
    exact real algebraic constants. Leading-coefficient and pole guards exclude
    specializations with smaller degree or undefined input coefficients.
    """
    parameters = validate_real_symbols(parameters, "parameters")
    validate_real_symbols((variable, *parameters), "variable and parameters")
    if len(degrees) != 2 or any(
        isinstance(d, bool) or not isinstance(d, int) or d < 2 for d in degrees
    ):
        raise ValueError("degrees must be a pair of integers at least two")
    outer_degree, inner_degree = degrees
    middle = (
        sp.Dummy("decomposition_t", real=True) if component_variable is None else component_variable
    )
    validate_real_symbols((variable, middle, *parameters), "variables and parameters")
    expression = polynomial.as_expr() if isinstance(polynomial, sp.Poly) else sp.sympify(polynomial)
    if expression.has(sp.Float) or expression.free_symbols - {variable, *parameters}:
        raise ValueError("polynomial must have exact coefficients and only declared parameters")
    numerator, denominator = sp.fraction(sp.together(expression))
    if variable in denominator.free_symbols:
        raise ValueError("input must be polynomial in the decomposition variable")
    try:
        _, options = sp.parallel_poly_from_expr(
            (numerator, denominator), variable, *parameters, extension=True
        )
        ground = options["domain"]
        if ground not in (sp.ZZ, sp.QQ) and (
            not (
                getattr(ground, "is_AlgebraicField", False) and ground.ext.as_expr().is_real is True
            )
        ):
            raise ValueError("coefficients must be exact real algebraic numbers")
        ground = ground.get_field()
        domain = ground.frac_field(*parameters) if parameters else ground
        target = sp.Poly(expression, variable, domain=domain)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed) as exc:
        raise ValueError("input must be polynomial over an exact real coefficient field") from exc
    if target.degree() != outer_degree * inner_degree:
        raise ValueError("polynomial degree must equal the product of the requested degrees")
    inner_unknowns = tuple(sp.Dummy(f"inner_{i}", real=True) for i in range(1, inner_degree))
    outer_unknowns = tuple(sp.Dummy(f"outer_{j}", real=True) for j in range(1, outer_degree))
    inner = variable**inner_degree + sum((c * variable**i for i, c in enumerate(inner_unknowns, 1)))
    outer = (
        target.LC() * middle**outer_degree
        + target.TC()
        + sum((c * middle**j for j, c in enumerate(outer_unknowns, 1)))
    )
    residual = sp.Poly(sp.expand(outer.subs(middle, inner) - expression), variable)
    equations, guards = _equations(residual.all_coeffs())
    leading_numerator = sp.fraction(sp.together(target.LC()))[0]
    guard = sp.And(sp.Ne(leading_numerator, 0), sp.Ne(denominator, 0), *guards)
    return PolynomialSystem(
        target,
        variable,
        parameters,
        tuple(degrees),
        middle,
        outer,
        inner,
        (*inner_unknowns, *outer_unknowns),
        equations,
        guard,
    )


def polynomial_pair(
    polynomial,
    variable: sp.Symbol,
    *,
    degrees: tuple[int, int],
    parameters: Sequence[sp.Symbol] = (),
    component_variable: sp.Symbol | None = None,
) -> PolynomialPair:
    """Return a real parameter condition and explicit fixed-degree components.

    Normalize the inner polynomial to be monic with constant zero. Its remaining
    coefficients are uniquely determined, top down, by the highest coefficients
    of the input. Then recover the outer coefficients top down. Vanishing of the
    residual coefficients is necessary and sufficient on the degree/pole guard.
    This avoids an unnormalized nonlinear QE search and works in characteristic
    zero. It does not search rational or general implicit algebraic components.
    """
    system = polynomial_system(
        polynomial,
        variable,
        degrees=degrees,
        parameters=parameters,
        component_variable=component_variable,
    )
    target = system.polynomial
    domain = target.domain
    m, n = system.degrees
    middle = system.component_variable
    inner = variable**n
    # The highest remaining term determines one inner coefficient without solving a nonlinear system.
    for i in range(n - 1, 0, -1):
        residual = target - sp.Poly(target.LC() * inner**m, variable, domain=domain)
        coefficient = sp.cancel(residual.nth((m - 1) * n + i) / (m * target.LC()))
        inner += coefficient * variable**i
    outer = target.LC() * middle**m + target.TC()
    for j in range(m - 1, 0, -1):
        residual = target - sp.Poly(outer.subs(middle, inner), variable, domain=domain)
        outer += residual.nth(n * j) * middle**j
    inner = sp.Poly(inner, variable, domain=domain).as_expr()
    outer_poly = sp.Poly(outer, middle, domain=domain)
    values = tuple(
        zip(
            system.unknowns[: n - 1],
            (sp.Poly(inner, variable, domain=domain).nth(i) for i in range(1, n)),
            strict=True,
        )
    )
    values += tuple(
        zip(system.unknowns[n - 1 :], (outer_poly.nth(j) for j in range(1, m)), strict=True)
    )
    residual = target - sp.Poly(outer_poly.as_expr().subs(middle, inner), variable, domain=domain)
    equations, guards = _equations(residual.all_coeffs())
    _, component_guards = _equations(
        (*sp.Poly(inner, variable, domain=domain).all_coeffs(), *outer_poly.all_coeffs())
    )
    condition = sp.And(system.guard, *guards, *component_guards, *equations)
    return PolynomialPair(system, condition, outer_poly.as_expr(), inner, values)


def _family_polynomial(polynomial, variable, parameters):
    """Exact rational-parameter polynomial, including constants and zero."""
    expression = polynomial.as_expr() if isinstance(polynomial, sp.Poly) else sp.sympify(polynomial)
    if expression.has(sp.Float) or expression.free_symbols - {variable, *parameters}:
        raise ValueError("polynomial must have exact coefficients and only declared parameters")
    numerator, denominator = sp.fraction(sp.together(expression))
    if variable in denominator.free_symbols:
        raise ValueError("input must be polynomial in the decomposition variable")
    _, options = sp.parallel_poly_from_expr(
        (numerator, denominator), variable, *parameters, extension=True
    )
    ground = options["domain"]
    if ground not in (sp.QQ, sp.ZZ) and (
        not (getattr(ground, "is_AlgebraicField", False) and ground.ext.as_expr().is_real is True)
    ):
        raise ValueError("coefficients must be exact real algebraic numbers")
    ground = ground.get_field()
    domain = ground.frac_field(*parameters) if parameters else ground
    return (sp.Poly(expression, variable, domain=domain), denominator, ground)


def _equalities(guard):
    return tuple(
        sp.fraction(sp.together(atom.lhs - atom.rhs))[0]
        for atom in sp.And.make_args(normalize_atoms(guard))
        if isinstance(atom, sp.Equality)
    )


def _reduce_on_guard(expression, guard, parameters, ground):
    """Reduce an expression modulo the conjunction of parameter equalities."""
    if not parameters:
        return sp.cancel(expression)
    equations = _equalities(guard)
    if not equations:
        return sp.cancel(expression)
    basis = sp.groebner(equations, *parameters, domain=ground, order="grevlex")
    numerator, denominator = sp.fraction(sp.together(sp.cancel(expression)))
    numerator = basis.reduce(numerator)[1]
    denominator = basis.reduce(denominator)[1]
    if denominator == 0:
        if not is_satisfiable(guard, parameters):
            return sp.Integer(0)
        raise ValueError("coefficient expression has a pole on its parameter guard")
    return sp.cancel(numerator / denominator)


def _prove_guard(guard, conclusion, parameters):
    """Prove an implication after equality reduction, using exact decisions if needed."""
    guard, conclusion = (normalize_atoms(guard), normalize_atoms(conclusion))
    if isinstance(conclusion, (sp.Equality, sp.Unequality)):
        numerator, denominator = sp.fraction(sp.together(conclusion.lhs - conclusion.rhs))
        if denominator != 1:
            return _prove_guard(guard, sp.Ne(denominator, 0), parameters) and _prove_guard(
                guard, conclusion.func(numerator, 0), parameters
            )
    if guard is sp.false or conclusion is sp.true or sp.And(guard, sp.Not(conclusion)) is sp.false:
        return True
    return bool(implies(guard, conclusion, parameters))


@dataclass(frozen=True)
class PolynomialCase:
    """All normalized fixed-degree components on one guarded degree stratum."""

    degree: int
    polynomial: sp.Expr
    components: tuple[PolynomialPair, ...]
    status: str

    def subs(self, assignments):
        """Specialize symbolic components after selecting their parameter guard."""
        components = []
        for component in self.components:
            system = component.system
            system = replace(
                system,
                polynomial=sp.Poly(system.polynomial.as_expr().subs(assignments), system.variable),
                parameters=tuple(p for p in system.parameters if p not in assignments),
                outer=system.outer.subs(assignments),
                inner=system.inner.subs(assignments),
                equations=tuple(eq.subs(assignments) for eq in system.equations),
                guard=system.guard.subs(assignments),
            )
            components.append(
                replace(
                    component,
                    system=system,
                    condition=component.condition.subs(assignments),
                    outer=component.outer.subs(assignments),
                    inner=component.inner.subs(assignments),
                    coefficient_values=tuple(
                        ((k, v.subs(assignments)) for k, v in component.coefficient_values)
                    ),
                )
            )
        return replace(
            self, polynomial=self.polynomial.subs(assignments), components=tuple(components)
        )


@dataclass(frozen=True)
class PolynomialFamily:
    """A complete guarded partition, including exceptional degrees and poles."""

    polynomial: sp.Expr
    variable: sp.Symbol
    parameters: tuple[sp.Symbol, ...]
    component_variable: sp.Symbol
    parameter_domain: sp.Expr
    stratification: ParameterStratifiedResult

    @property
    def branches(self):
        """Disjoint guarded cases covering the declared parameter domain."""
        return self.stratification.branches

    def select(self, assignments):
        """Specialize the unique case containing an exact parameter assignment."""
        assignments = _parameter_assignments(assignments, self.parameters)
        if set(assignments) - set(self.parameters) or any(
            sp.sympify(v).has(sp.Float) for v in assignments.values()
        ):
            raise ValueError("selection requires exact values for declared parameters")
        return self.stratification.select(assignments)

    def verify(self):
        """Replay the parameter partition and all guarded component identities."""
        return verify_family(self)


def polynomial_family(
    polynomial,
    variable: sp.Symbol,
    *,
    parameters: Sequence[sp.Symbol] = (),
    parameter_domain=sp.true,
    component_variable: sp.Symbol | None = None,
):
    """Reconstruct normalized components across every exact degree-drop stratum.

    Branches include decomposable, indecomposable, degree-below-two, and undefined
    coefficient-pole cases. Every ordered nontrivial factorization of the actual
    degree is considered. Guards are exact formulas, never sample substitutions.
    Empty guards may remain symbolically; they assert no result at any real point.
    """
    parameters = validate_real_symbols(parameters, "parameters")
    validate_real_symbols((variable, *parameters), "variable and parameters")
    middle = (
        sp.Dummy("parametric_decomposition_t", real=True)
        if component_variable is None
        else component_variable
    )
    validate_real_symbols((variable, middle, *parameters), "variables and parameters")
    parameter_domain = sp.sympify(parameter_domain)
    if parameter_domain.free_symbols - set(parameters):
        raise ValueError("parameter_domain may contain only declared parameters")
    target, denominator, ground = _family_polynomial(polynomial, variable, parameters)
    branches = []
    undefined = normalize_atoms(sp.And(parameter_domain, sp.Eq(denominator, 0)))
    if undefined is not sp.false:
        branches.append(ConditionalBranch(undefined, PolynomialCase(-1, sp.nan, (), "undefined")))
    guard = normalize_atoms(sp.And(parameter_domain, sp.Ne(denominator, 0)))
    upper = -1 if target.is_zero else int(target.degree())
    for degree in range(upper, -1, -1):
        if guard is sp.false:
            break
        leading = _reduce_on_guard(target.nth(degree), guard, parameters, ground)
        numerator = sp.fraction(sp.together(leading))[0]
        current = normalize_atoms(sp.And(guard, sp.Ne(numerator, 0)))
        if current is not sp.false:
            coefficients = [
                _reduce_on_guard(target.nth(i), guard, parameters, ground)
                for i in range(degree + 1)
            ]
            reduced = sum((c * variable**i for i, c in enumerate(coefficients)))
            candidates = tuple(
                polynomial_pair(
                    reduced,
                    variable,
                    degrees=(degree // n, n),
                    parameters=parameters,
                    component_variable=middle,
                )
                for n in range(2, degree)
                if degree % n == 0
            )
            cases = [(current, ())]
            for candidate in candidates:
                split = []
                for case_guard, components in cases:
                    yes = normalize_atoms(sp.And(case_guard, candidate.condition))
                    no = normalize_atoms(sp.And(case_guard, sp.Not(candidate.condition)))
                    if yes is not sp.false:
                        split.append((yes, (*components, candidate)))
                    if no is not sp.false:
                        split.append((no, components))
                cases = split
            for condition, components in cases:
                status = (
                    "decomposable"
                    if components
                    else "degree_below_two"
                    if degree < 2
                    else "indecomposable"
                )
                branches.append(
                    ConditionalBranch(
                        condition, PolynomialCase(degree, reduced, components, status)
                    )
                )
        guard = normalize_atoms(sp.And(guard, sp.Eq(numerator, 0)))
    if guard is not sp.false:
        branches.append(
            ConditionalBranch(guard, PolynomialCase(-1, sp.Integer(0), (), "degree_below_two"))
        )
    stratification = ParameterStratifiedResult(
        parameters,
        tuple(branches),
        parameter_domain,
        complete=True,
        disjoint=True,
        certified=True,
        method="coefficient_degree_stratification",
    )
    return PolynomialFamily(
        target.as_expr(), variable, parameters, middle, parameter_domain, stratification
    )


def verify_family(result: PolynomialFamily) -> bool:
    """Replay branch coverage, actual degrees, component guards and identities."""
    if result.stratification.parameters != result.parameters:
        return False
    if not verify_parameter_stratification(
        result.stratification, parameter_domain=result.parameter_domain
    ).verified:
        return False
    target, denominator, ground = _family_polynomial(
        result.polynomial, result.variable, result.parameters
    )
    for branch in result.branches:
        guard, case = (branch.condition, branch.value)
        if case.status == "undefined":
            if not _prove_guard(guard, sp.Eq(denominator, 0), result.parameters):
                return False
            continue
        if not _prove_guard(guard, sp.Ne(denominator, 0), result.parameters):
            return False
        residual = sp.Poly(target.as_expr() - case.polynomial, result.variable)
        for coefficient in residual.all_coeffs():
            reduced = _reduce_on_guard(coefficient, guard, result.parameters, ground)
            if reduced != 0 and (
                not _prove_guard(
                    guard, sp.Eq(sp.fraction(sp.together(reduced))[0], 0), result.parameters
                )
            ):
                return False
        actual = sp.Poly(case.polynomial, result.variable)
        if case.degree < 0:
            if not actual.is_zero or case.components or case.status != "degree_below_two":
                return False
            continue
        if actual.degree() != case.degree or not _prove_guard(
            guard, sp.Ne(actual.LC(), 0), result.parameters
        ):
            return False
        expected_pairs = tuple(
            (case.degree // n, n) for n in range(2, case.degree) if case.degree % n == 0
        )
        supplied = {c.system.degrees: c for c in case.components}
        if len(supplied) != len(case.components) or set(supplied) - set(expected_pairs):
            return False
        for pair in expected_pairs:
            candidate = polynomial_pair(
                case.polynomial,
                result.variable,
                degrees=pair,
                parameters=result.parameters,
                component_variable=result.component_variable,
            )
            if pair not in supplied:
                if not _prove_guard(guard, sp.Not(candidate.condition), result.parameters):
                    return False
                continue
            component = supplied[pair]
            if (
                component.system.variable != result.variable
                or component.system.component_variable != result.component_variable
            ):
                return False
            inner = sp.Poly(component.inner, result.variable)
            outer = sp.Poly(component.outer, result.component_variable)
            if (
                inner.degree() != pair[1]
                or outer.degree() != pair[0]
                or inner.LC() != 1
                or (inner.nth(0) != 0)
            ):
                return False
            if not _prove_guard(guard, component.condition, result.parameters):
                return False
            difference = sp.Poly(component.composition - result.polynomial, result.variable)
            for coefficient in difference.all_coeffs():
                reduced = _reduce_on_guard(coefficient, guard, result.parameters, ground)
                if reduced != 0 and (
                    not _prove_guard(
                        guard, sp.Eq(sp.fraction(sp.together(reduced))[0], 0), result.parameters
                    )
                ):
                    return False
        expected_status = (
            "decomposable"
            if supplied
            else "degree_below_two"
            if case.degree < 2
            else "indecomposable"
        )
        if case.status != expected_status:
            return False
    return True


__all__ = [
    "PolynomialSystem",
    "PolynomialPair",
    "PolynomialCase",
    "PolynomialFamily",
    "polynomial_system",
    "polynomial_pair",
    "polynomial_family",
    "verify_family",
]
