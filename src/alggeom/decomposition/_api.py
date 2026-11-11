"""Mathematical entry points for functional decomposition"""

import sympy as sp

from ._validation import IncompleteDecompositionError
from .chains import algebraic_chain, rational_chain
from .common import common_components
from .rational import search_rational


def common_right_component(functions, variables, *, constants=sp.QQ, parameter=None):
    """Return a Lambda generating the joint field of the supplied rational functions.

    Return None for constant-only inputs or a joint
    field of transcendence degree greater than one. A degree-one Lambda is a
    valid generator. Use common.common_components for outer maps and evidence.
    """
    return common_components(functions, variables, constants=constants, parameter=parameter).inner


def decompose(expression, variables, *, constants=sp.QQ, **limits):
    """Return an outer-to-inner tuple of SymPy Lambdas for one resolved chain.

    ``variables`` is a Symbol or an ordered sequence of input coordinates. Only
    the last map may have several arguments. Constant and indecomposible inputs
    give a one-element tuple. Search limits are those of ``rational_chain``.
    An unfinished chain raises ``IncompleteDecompositionError`` with its
    certificate in ``error.result``.
    """
    chain = rational_chain(expression, variables, constants=constants, **limits)
    if not chain.complete:
        raise IncompleteDecompositionError(chain)
    return tuple(
        sp.Lambda(component.variables, component.expression) for component in chain.components
    )


def decompositions(expression, variables, *, parameter=None, constants=sp.QQ, **limits):
    """Return every nontrivial two-map pair over the declared constant field.

    Each pair is ``(outer, inner)`` with SymPy Lambda maps. Answers are distinct
    up to a fractional-linear change of the intermediate coordinate. An empty
    tuple proves that no such pair exists in this scope. Incomplete searches
    raise ``IncompleteDecompositionError``; their partial answers remain in
    ``error.result``. Limits are those of ``rational.search_rational``.
    """
    result = search_rational(
        expression, variables, parameter=parameter, constants=constants, **limits
    )
    if not result.complete:
        raise IncompleteDecompositionError(result)
    return tuple(
        (sp.Lambda(pair.parameter, pair.outer), sp.Lambda(pair.variables, pair.inner))
        for pair in result.decompositions
    )


def decompose_relation(relation, variables, *, constants=sp.QQ, **limits):
    """Return an ordered tuple of implicit component polynomials.

    ``variables=(output,input)`` fixes the orientation. Adjacent polynomials
    share an intermediate coordinate. Their iterated resultant is a nonzero
    constant times a positive power of the original relation. Completeness
    concerns the selected chain under coefficient-field and rational-graph
    search. An unresolved leaf raises ``IncompleteDecompositionError``; the
    exact partial chain remains available in ``error.result``.
    """
    chain = algebraic_chain(relation, variables, constants=constants, **limits)
    if not chain.complete:
        raise IncompleteDecompositionError(chain)
    return chain.components


def decomposition_condition(polynomial, variable, *, degrees, parameters=()):
    """Return the exact real parameter condition for fixed polynomial degrees.

    ``degrees=(outer_degree,inner_degree)`` fixes the requested composition.
    The condition excludes coefficient poles and degree drops. For explicit
    components or all degree strata use ``parameters.polynomial_pair`` or
    ``parameters.polynomial_family``. Requires the real-domain dependency.
    """
    from .parameters import polynomial_pair

    return polynomial_pair(polynomial, variable, degrees=degrees, parameters=parameters).condition


def real_equivalent(chain, **options):
    """Decide equality of a chain's selected real graph on a declared domain.

    ``options`` are those of ``real.verify_chain``. Equality of graphs is
    independent of totality and uniqueness; use that specialist interface for
    their individual results. Unsupported exact decisions raise their solver
    error. Requires the real-domain dependency.
    """
    from .real import verify_chain

    return verify_chain(chain, **options).equivalent
