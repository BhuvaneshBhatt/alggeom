"""Implicit composition certificates and supplied-factor reconstruction."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from algroots.curves import CurveRationalityResult
from algroots.function_fields import (
    AlgebraicCoefficientField,
    CoefficientFieldParametrization,
    coefficient_field,
)
from algroots.quotient import QuotientAlgebra
from algroots.subfields import field_normal_form
from sympy.polys.polyerrors import BasePolynomialError

from ._validation import (
    DecompositionError,
    constant_field,
    normalize_relation,
    parse_relation,
    validate_integer,
    validate_symbols,
)


@dataclass(frozen=True)
class AlgebraicDecomposition:
    """An irreducible resultant decomposition of an implicit relation."""

    relation: sp.Poly
    outer: sp.Poly
    inner: sp.Poly
    exponent: int
    unit: sp.Expr

    @property
    def output(self):
        """Output coordinate of the implicit relation."""
        return self.relation.gens[0]

    @property
    def input(self):
        """Input coordinate of the implicit relation."""
        return self.relation.gens[1]

    @property
    def intermediate(self):
        """Coordinate shared by the two components."""
        return self.outer.gens[1]

    @property
    def constants(self):
        """Declared constant field of the implicit relation."""
        return self.relation.domain

    def verify(self):
        """Check all polynomials, nontriviality, degrees, power and scalar unit."""
        try:
            x, z, y = (self.output, self.input, self.intermediate)
            validate_symbols((x, y, z))
            if constant_field(self.constants) != self.constants:
                return False
            validate_integer(self.exponent, "exponent", 1)
            if self.outer.gens != (x, y) or self.inner.gens != (y, z):
                return False
            if any(
                p.domain != self.constants
                or not p.is_irreducible
                or any(p.degree(v) <= 0 for v in p.gens)
                for p in (self.relation, self.outer, self.inner)
            ):
                return False
            if any(max(p.degree(v) for v in p.gens) <= 1 for p in (self.outer, self.inner)):
                return False
            if sp.sympify(self.unit).has(sp.Float):
                return False
            unit = self.constants.convert(self.unit)
            if (
                not unit
                or self.exponent * self.relation.degree(x)
                != self.outer.degree(x) * self.inner.degree(y)
                or self.exponent * self.relation.degree(z)
                != self.outer.degree(y) * self.inner.degree(z)
            ):
                return False
            resultant = sp.Poly(
                sp.resultant(self.outer.as_expr(), self.inner.as_expr(), y),
                x,
                z,
                domain=self.constants,
            )
            return resultant == self.relation**self.exponent * self.constants.to_sympy(unit)
        except (
            ValueError,
            TypeError,
            AttributeError,
            IndexError,
            BasePolynomialError,
            ZeroDivisionError,
        ):
            return False


def normalize_factor(expression, variable, ambient):
    """Reduce coefficients in the ambient field and make the factor monic."""
    try:
        polynomial = sp.Poly(expression, variable)
        coefficients = [field_normal_form(ambient, c) for c in polynomial.all_coeffs()]
        while coefficients and coefficients[0] == 0:
            coefficients.pop(0)
        if len(coefficients) <= 1:
            raise DecompositionError("candidate must have positive degree")
        inverse = ambient.inverse(coefficients[0])
        coefficients = tuple(field_normal_form(ambient, c * inverse) for c in coefficients)
    except (ValueError, TypeError, BasePolynomialError, ZeroDivisionError) as exc:
        raise DecompositionError("candidate must be a polynomial over the ambient field") from exc
    return sp.Add(
        *(c * variable ** (len(coefficients) - i - 1) for i, c in enumerate(coefficients))
    )


def factor_remainder(relation, candidate, variable, ambient):
    """Divide by a monic candidate using ambient-field coefficient reduction."""
    left = list(sp.Poly(relation, variable).all_coeffs())
    right = sp.Poly(candidate, variable).all_coeffs()
    while len(left) >= len(right):
        lead = field_normal_form(ambient, left[0])
        for i, c in enumerate(right):
            left[i] = field_normal_form(ambient, left[i] - lead * c)
        left.pop(0)
    return tuple(field_normal_form(ambient, c) for c in left)


@dataclass(frozen=True)
class AlgebraicReconstruction:
    """Reconstruction from one supplied factor; never a complete discovery claim."""

    relation: sp.Poly
    ambient: QuotientAlgebra
    candidate: sp.Expr
    parameter: sp.Symbol
    coefficient_field: AlgebraicCoefficientField
    rationality: CurveRationalityResult
    parametrization: CoefficientFieldParametrization | None
    decomposition: AlgebraicDecomposition | None
    reason: str | None

    @property
    def complete(self):
        """False: one supplied factor does not establish discovery completeness."""
        return False

    @property
    def status(self):
        """Whether the supplied factor yields a nontrivial composition."""
        return "decomposable" if self.decomposition is not None else "candidate_rejected"

    @property
    def decompositions(self):
        """The reconstructed composition, or an empty tuple for a rejected candidate."""
        return () if self.decomposition is None else (self.decomposition,)

    def verify(self):
        """Replay the supplied-factor embedding, rational maps and final identity."""
        try:
            x, z = self.relation.gens
            if (
                parse_relation(self.relation.as_expr(), (x, z), self.relation.domain)
                != self.relation
            ):
                return False
            validate_ambient(self.ambient, self.relation.domain, z, x, self.parameter)
            if normalize_factor(self.candidate, x, self.ambient) != self.candidate or any(
                factor_remainder(self.relation.as_expr(), self.candidate, x, self.ambient)
            ):
                return False
            E = self.coefficient_field
            coefficients = tuple(
                field_normal_form(self.ambient, c)
                for c in sp.Poly(self.candidate, x).all_coeffs()[1:]
            )
            if (
                E.ambient is not self.ambient
                or E.constants != self.relation.domain
                or E.generators != coefficients
                or (not E.verify())
            ):
                return False
            model = E.curve_model()
            if self.rationality.normalization.curve != model.curve or not self.rationality.verify():
                return False
            if not self.rationality.rational:
                return (
                    self.parametrization is None
                    and self.decomposition is None
                    and (self.reason == self.rationality.status)
                )
            P = self.parametrization
            if (
                P.model != model
                or P.parametrization != self.rationality.parametrization
                or P.parameter != self.parameter
                or (not P.verify())
            ):
                return False
            G, H = reconstruct_components(self.relation, P, self.parameter, self.ambient)
            rejection = component_rejection(G, H)
            if rejection is not None:
                return self.decomposition is None and self.reason == rejection
            D = self.decomposition
            return (
                D is not None
                and D.relation == self.relation
                and (D.outer == G)
                and (D.inner == H)
                and (self.reason is None)
                and D.verify()
            )
        except (ValueError, TypeError, AttributeError, BasePolynomialError, ZeroDivisionError):
            return False


def validate_ambient(ambient, constants, z, x, y):
    """Require a certified finite field over constants(input) with fresh coordinates."""
    if (
        not isinstance(ambient, QuotientAlgebra)
        or ambient.domain != constants.frac_field(z)
        or (not ambient.is_field)
    ):
        raise DecompositionError("ambient must be a certified quotient field over constants(input)")
    validate_symbols((x, y, z, *ambient.variables))


def reconstruct_components(relation, parametrization, parameter, ambient):
    """Recover the outer factor and the minimal relation of its coefficient generator."""
    x, z = relation.gens
    images = parametrization.coefficient_images
    expression = x ** len(images) + sum(
        (c * x ** (len(images) - i - 1) for i, c in enumerate(images))
    )
    G = normalize_relation(expression, (x, parameter), relation.domain)
    H = normalize_relation(
        ambient.minimal_polynomial(parametrization.generator, parameter).as_expr(),
        (parameter, z),
        relation.domain,
    )
    return (G, H)


def component_rejection(G, H):
    """Identify reducible or degree-one reconstructed components."""
    if any(not p.is_irreducible for p in (G, H)):
        return "reducible_candidate"
    if any(
        any(p.degree(v) <= 0 for v in p.gens) or max(p.degree(v) for v in p.gens) <= 1
        for p in (G, H)
    ):
        return "trivial_component"
    return None


def reconstruct_algebraic(
    relation,
    variables,
    *,
    candidate,
    ambient,
    parameter=None,
    constants=sp.QQ,
    point=None,
    max_steps=256,
    max_dimension=256,
):
    """Reconstruct an implicit composition from a factor in a declared finite field."""
    x, z = tuple(variables)
    parameter = sp.Dummy("intermediate") if parameter is None else parameter
    constants = constant_field(constants)
    if not isinstance(ambient, QuotientAlgebra):
        raise DecompositionError("ambient must be a QuotientAlgebra field")
    validate_ambient(ambient, constants, z, x, parameter)
    F = parse_relation(relation, (x, z), constants)
    g = normalize_factor(candidate, x, ambient)
    if any(factor_remainder(F.as_expr(), g, x, ambient)):
        raise DecompositionError("candidate does not divide the input in the ambient field")
    coefficients = sp.Poly(g, x).all_coeffs()[1:]
    # Do not adjoin the input coordinate: that can destroy a proper coefficient subfield.
    E = coefficient_field(coefficients, ambient=ambient, constants=constants)
    if E.transcendence_degree != 1:
        raise DecompositionError("candidate coefficients must generate a function field")
    R = E.rationality(parameter, point=point, max_steps=max_steps, max_dimension=max_dimension)
    if not R.rational:
        return AlgebraicReconstruction(F, ambient, g, parameter, E, R, None, None, R.status)
    P = CoefficientFieldParametrization(E.curve_model(), R.parametrization)
    G, H = reconstruct_components(F, P, parameter, ambient)
    reason = component_rejection(G, H)
    D = None
    if reason is None:
        resultant = sp.Poly(
            sp.resultant(G.as_expr(), H.as_expr(), parameter), x, z, domain=constants
        )
        exponent, remainder = divmod(G.degree(x) * H.degree(parameter), F.degree(x))
        if remainder or exponent < 1:
            raise DecompositionError("component degrees do not give an integral resultant power")
        unit = constants.to_sympy(constants.convert(resultant.LC() / (F**exponent).LC()))
        D = AlgebraicDecomposition(F, G, H, exponent, unit)
        if not D.verify():
            raise DecompositionError(
                "candidate produced extra resultant components or failed exact reconstruction"
            )
    return AlgebraicReconstruction(F, ambient, g, parameter, E, R, P, D, reason)


__all__ = ["AlgebraicDecomposition", "AlgebraicReconstruction", "reconstruct_algebraic"]
