"""Prime and finite-extension coefficient fields for exact algebraic geometry."""

from algroots import AlgebraicFunctionField
from algroots.finite_fields import field_order as _field_order
from algroots.finite_fields import finite_field as _finite_field


def finite_field(characteristic, degree=1, *, name="a", modulus=None):
    """Return GF(p) or GF(p^degree) with a validated irreducible defining polynomial.

    A supplied modulus determines the extension; otherwise python-flint chooses
    one. The returned object is a SymPy coefficient domain. GF(p^degree) is
    never interpreted as the composite residue ring Z/(p^degree).
    """
    return _finite_field(characteristic, degree, name=name, modulus=modulus)


def field_order(field):
    """Return the cardinality of a supported finite coefficient field."""
    return _field_order(field)


__all__ = ["AlgebraicFunctionField", "finite_field", "field_order"]
