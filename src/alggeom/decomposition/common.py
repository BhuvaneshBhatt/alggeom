"""A shared rational intermediate coordinate for a tuple of functions."""

from dataclasses import dataclass

import sympy as sp
from algroots.function_fields import LurothGeneratorResult, luroth_generator
from sympy.polys.polyerrors import BasePolynomialError

from ._validation import constant_field, input_variables, visible_denominators
from .rational import normalize_rational


@dataclass(frozen=True)
class CommonDecomposition:
    """Constructive equality K(f_1,...,f_m)=K(h), including negative field evidence.

    Individual outer maps may be linear. The inner map generates exactly the
    joint field, so a degree-one inner map is a valid answer. Constant inputs
    have no distinguished nonconstant generator; higher transcendence degree
    rules out a common scalar rational intermediate.
    """

    functions: tuple[sp.Expr, ...]
    field: LurothGeneratorResult

    @property
    def status(self):
        """Rational, constant, or higher_transcendence_degree."""
        return self.field.status

    @property
    def inner(self):
        """A SymPy Lambda generating the joint field, or None when none is selected."""
        return (
            sp.Lambda(self.field.variables, self.field.generator)
            if self.status == "rational"
            else None
        )

    @property
    def outers(self):
        """The univariate Lambda maps f_i=g_i(h), or an empty tuple without a generator."""
        return (
            tuple(sp.Lambda(self.field.parameter, image) for image in self.field.coefficient_images)
            if self.status == "rational"
            else ()
        )

    @property
    def input_guard(self):
        """Denominator exclusions visible in the original input functions."""
        return sp.And(*(sp.Ne(p, 0) for f in self.functions for p in visible_denominators(f)))

    def verify(self):
        """Replay field generation, inverse coordinates and each original function."""
        try:
            values = tuple(
                normalize_rational(f, self.field.variables, self.field.constants)
                for f in self.functions
            )
            return values == self.field.generators and self.field.verify()
        except (ValueError, TypeError, AttributeError, BasePolynomialError, ZeroDivisionError):
            return False


def common_components(functions, variables, *, parameter=None, constants=sp.QQ):
    """Compute all outer maps through one generator of the joint rational function field.

    Equality is in the function field; input_guard records original affine
    poles. This is a generator of the joint field, not an enumeration of all
    intermediate fields or all possible shared substitutions.
    """
    variables = input_variables(variables)
    constants = constant_field(constants)
    functions = tuple(f.as_expr() if isinstance(f, sp.Poly) else sp.sympify(f) for f in functions)
    if not functions:
        raise ValueError("at least one rational function is required")
    field = luroth_generator(
        functions, variables=variables, constants=constants, parameter=parameter
    )
    return CommonDecomposition(functions, field)


__all__ = ["CommonDecomposition", "common_components"]
