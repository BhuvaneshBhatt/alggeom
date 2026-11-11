"""Parameter-sensitive invariants of a univariate polynomial fiber family."""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class FiberInvariants:
    """Geometric dimension and distinct-point count of one fiber over an algebraic closure."""

    dimension: int
    geometric_points: object


@dataclass(frozen=True)
class FiberStratum:
    """A parameter condition with exact geometric fiber invariants."""

    condition: sp.Expr
    invariants: FiberInvariants


@dataclass(frozen=True)
class FiberStratification:
    """Generic invariants and a complete finite exceptional locus in A^1.

    The complement of ``exceptional_polynomial=0`` is a certified open stratum.
    ``strata()`` isolates its distinct algebraic zeros and computes their
    exact fiber invariants, covering all geometric parameter points.
    """

    equation: sp.Expr
    variable: sp.Symbol
    parameter: sp.Symbol
    exceptional_polynomial: sp.Expr
    generic: FiberInvariants

    @property
    def generic_condition(self):
        """A sufficient Zariski-open condition for the generic invariants."""
        return sp.Ne(self.exceptional_polynomial, 0)

    def strata(self, *, max_exceptional_roots=32, max_root_degree=2):
        """Enumerate a full stratification of the one-dimensional parameter line.

        Over an algebraic closure the exceptional polynomial has finitely many
        roots. Each distinct root forms one closed stratum, with its fiber
        invariants computed by exact specialization. The complement is one
        open stratum. A declared root-count limit prevents accidental large
        algebraic isolations. The default degree limit excludes higher-degree
        RootOf specializations that can be prohibitively expensive.
        """
        exceptional = sp.Poly(self.exceptional_polynomial, self.parameter, domain=sp.QQ)
        if exceptional.degree() > max_exceptional_roots:
            raise ValueError("exceptional polynomial exceeds max_exceptional_roots")
        if exceptional.degree() <= 0:
            return (FiberStratum(sp.true, self.generic),)
        factors = sp.factor_list(exceptional)[1]
        if any(factor.degree() > max_root_degree for factor, _ in factors):
            raise NotImplementedError("exceptional root degree exceeds max_root_degree")
        roots = tuple(dict.fromkeys(exceptional.all_roots()))
        return (
            FiberStratum(self.generic_condition, self.generic),
            *(FiberStratum(sp.Eq(self.parameter, root), self.specialize(root)) for root in roots),
        )

    def specialize(self, value):
        """Compute the geometric dimension and number of distinct points exactly."""
        value = sp.sympify(value)
        if value.has(sp.Float) or value.is_algebraic is not True:
            raise ValueError("specialization requires an exact algebraic parameter value")
        # RootOf objects commonly bind a symbol named x. Using a fresh
        # polynomial generator prevents capture by the bound root variable.
        fiber_coordinate = sp.Dummy("fiber_coordinate")
        expression = sp.expand(self.equation.subs(self.parameter, value))
        polynomial = sp.Poly(
            expression.xreplace({self.variable: fiber_coordinate}), fiber_coordinate, extension=True
        )
        if polynomial.is_zero:
            return FiberInvariants(1, sp.oo)
        if polynomial.degree() == 0:
            return FiberInvariants(-1, 0)
        return FiberInvariants(0, int(polynomial.sqf_part().degree()))


def fiber_stratification(equation, variable, parameter):
    """Stratify univariate polynomial or multidimensional affine-linear fibers.

    A sequence of equations and sequences of fiber/parameter symbols select
    the determinantal algorithm, which supports any finite number of
    variables and polynomial parameters but requires linear fiber equations.
    Scalar arguments select the discriminant-based univariate algorithm.

    The discriminant of the generic squarefree polynomial, its coefficient
    denominators and leading coefficient, and the content of the original
    family determine a certified open set with constant geometric root count.
    This avoids confusing repeated generic roots with specialization failures.
    """
    if any(isinstance(arg, (tuple, list)) for arg in (equation, variable, parameter)):
        from .linear_strata import linear_fiber_stratification

        variables = (variable,) if isinstance(variable, sp.Symbol) else tuple(variable)
        parameters = (parameter,) if isinstance(parameter, sp.Symbol) else tuple(parameter)
        equations = tuple(equation) if isinstance(equation, (tuple, list)) else (equation,)
        if len(equations) == 1:
            from .nonlinear_strata import hypersurface_fiber_stratification

            return hypersurface_fiber_stratification(equations[0], variables, parameters)
        try:
            return linear_fiber_stratification(equations, variables, parameters)
        except NotImplementedError:
            from .nonlinear_strata import separable_fiber_stratification

            try:
                return separable_fiber_stratification(equations, variables, parameters)
            except NotImplementedError:
                from .comprehensive_strata import comprehensive_fiber_stratification

                return comprehensive_fiber_stratification(
                    equations, variables, parameters, count_points=False
                )
    variable, parameter = sp.sympify(variable), sp.sympify(parameter)
    if (
        not isinstance(variable, sp.Symbol)
        or not isinstance(parameter, sp.Symbol)
        or variable == parameter
    ):
        raise ValueError("one distinct variable and parameter Symbol are required")
    expression = sp.sympify(equation)
    if expression.has(sp.Float) or expression.free_symbols - {variable, parameter}:
        raise ValueError("a polynomial over QQ in exactly the declared symbols is required")
    polynomial = sp.Poly(expression, variable, parameter, domain=sp.QQ)
    if polynomial.is_zero:
        return FiberStratification(
            expression, variable, parameter, sp.S.One, FiberInvariants(1, sp.oo)
        )
    generic_poly = sp.Poly(expression, variable, domain=sp.QQ.frac_field(parameter))
    if generic_poly.degree() == 0:
        generic = FiberInvariants(-1, 0)
    else:
        generic = FiberInvariants(0, int(generic_poly.sqf_part().degree()))
    candidates = []
    univariate = sp.Poly(expression, variable)
    coefficients = univariate.all_coeffs()
    content = sp.gcd_list(coefficients)
    candidates.append(content)
    candidates.append(univariate.LC())
    squarefree = generic_poly.sqf_part()
    if squarefree.degree() > 0:
        candidates.append(squarefree.LC())
    if squarefree.degree() > 1:
        candidates.append(sp.discriminant(squarefree.as_expr(), variable))
    for coefficient in squarefree.all_coeffs():
        candidates.append(sp.together(coefficient).as_numer_denom()[1])
    factors = []
    for candidate in candidates:
        numerator, denominator = sp.together(candidate).as_numer_denom()
        factors.extend((numerator, denominator))
    product = sp.prod(
        sp.Poly(factor, parameter, domain=sp.QQ).sqf_part().as_expr() for factor in factors
    )
    if product == 0:
        raise RuntimeError("unexpected zero generic discriminant factor")
    return FiberStratification(expression, variable, parameter, sp.factor(product), generic)


__all__ = ["FiberInvariants", "FiberStratum", "FiberStratification", "fiber_stratification"]
