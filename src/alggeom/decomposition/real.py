"""Exact real graph comparison, branch selection, totality and uniqueness."""

from __future__ import annotations

import sympy as sp
from semialg.function_verification import compare_real_graphs
from semialg.qe import QuantifierEliminationResult, quantifier_eliminate
from semialg.simplify.atoms import normalize_atoms

from ._validation import validate_real_symbols, visible_denominators
from .algebraic import AlgebraicDecomposition
from .chains import AlgebraicChain, RationalChain, algebraic_chain, rational_chain
from .rational import RationalDecomposition


def _rational_chain_projection(graphs, components, inputs, middle, output):
    """Eliminate uniquely determined rational intermediates with their guards."""
    substitutions = {}
    for index in range(len(components) - 1, 0, -1):
        sources = inputs if index == len(components) - 1 else (substitutions[middle[index]],)
        renaming = dict(zip(components[index].variables, sources, strict=True))
        substitutions[middle[index - 1]] = sp.cancel(
            components[index].expression.xreplace(renaming)
        )
    formula = sp.And(*graphs).xreplace(substitutions)
    replacements, guards = ({}, [])
    for atom in formula.atoms(sp.core.relational.Relational):
        numerator, denominator = sp.cancel(atom.lhs - atom.rhs).as_numer_denom()
        guards.append(sp.Ne(denominator, 0))
        # Multiplication by denominator squared preserves the sign on the pole-free graph.
        polynomial = (
            numerator if isinstance(atom, (sp.Equality, sp.Unequality)) else numerator * denominator
        )
        replacements[atom] = atom.func(polynomial, 0)
    formula = normalize_atoms(sp.And(formula.xreplace(replacements), *guards))
    return QuantifierEliminationResult(
        formula=formula,
        free_variables=(*inputs, output),
        quantified_variables=middle,
        quantifiers=tuple(("exists", v) for v in middle),
        method="guarded-rational-chain-substitution",
        presolve_substitutions=tuple(substitutions.items()),
        notes=("each intermediate is uniquely determined by its guarded rational component",),
    )


def verify_chain(
    chain,
    *,
    output=None,
    intermediates=None,
    domain=sp.true,
    target_selector=sp.true,
    component_selectors=None,
    check_functions=True,
):
    """Compare a decomposition chain with its selected affine real composition.

    Components and selectors are ordered outermost to innermost. Rational
    chains may name their intermediate coordinates explicitly; implicit chains
    already fix those coordinates in ``chain.components``. Each selector may
    use only its component's input and output coordinates. Chain search need
    not be complete, but its exact composition certificate must verify.
    """
    if not isinstance(chain, (RationalChain, AlgebraicChain)):
        raise ValueError("a rational or algebraic decomposition chain is required")
    if not chain.verify():
        raise ValueError("the decomposition chain failed verification")
    if chain.constants.is_AlgebraicField and chain.constants.ext.as_expr().is_real is not True:
        raise ValueError("real graph verification requires a real embedding of the constants")
    rational = isinstance(chain, RationalChain)
    components = chain.components
    count = len(components)
    if rational:
        inputs = chain.variables
        z = sp.Dummy("chain_output", real=True) if output is None else output
        middle = (
            tuple(sp.Dummy(f"chain_middle_{i}", real=True) for i in range(count - 1))
            if intermediates is None
            else tuple(intermediates)
        )
        if len(middle) != count - 1:
            raise ValueError("one intermediate coordinate per adjacent pair is required")
    else:
        z, x = chain.variables
        inputs = (x,)
        middle = tuple(p.gens[1] for p in components[:-1])
        if output is not None and output != z:
            raise ValueError("an implicit chain already fixes its output variable")
        if intermediates is not None and tuple(intermediates) != middle:
            raise ValueError("an implicit chain already fixes its intermediate variables")
    validate_real_symbols((*inputs, z, *middle), "chain variables")
    x, parameters = (inputs[-1], inputs[:-1])
    selectors = (sp.true,) * count if component_selectors is None else tuple(component_selectors)
    if len(selectors) != count:
        raise ValueError("one selector per component is required")
    target_selector, domain = (sp.sympify(target_selector), sp.sympify(domain))
    if target_selector.free_symbols - {*inputs, z}:
        raise ValueError("target selector contains undeclared symbols")
    if domain.free_symbols - set(inputs):
        raise ValueError("domain may contain only the input coordinates")

    def graph(expression, target):
        numerator, denominator = sp.cancel(expression).as_numer_denom()
        guards = [sp.Ne(denominator, 0)]
        guards.extend(sp.Ne(p, 0) for p in visible_denominators(expression))
        return sp.And(sp.Eq(target * denominator - numerator, 0), *guards)

    target = (
        sp.And(graph(chain.original_expression, z), chain.input_guard)
        if rational
        else sp.Eq(chain.relation.as_expr(), 0)
    )
    graphs = []
    for index, (component, selector) in enumerate(zip(components, selectors, strict=True)):
        target_coordinate = z if index == 0 else middle[index - 1]
        source_coordinates = inputs if index == count - 1 else (middle[index],)
        selector = sp.sympify(selector)
        if selector.free_symbols - {target_coordinate, *source_coordinates}:
            raise ValueError("component selector contains undeclared symbols")
        if rational:
            renaming = dict(zip(component.variables, source_coordinates, strict=True))
            formula = graph(component.expression.xreplace(renaming), target_coordinate)
        else:
            formula = sp.Eq(component.as_expr(), 0)
        graphs.append(sp.And(formula, selector))
    renaming = {v: sp.Dummy("bound_chain_middle", real=True) for v in middle}
    bound = tuple(renaming.values())
    variables = (*inputs, z)
    target_qe = quantifier_eliminate(
        sp.And(target, target_selector), [], variables=variables, return_result=True
    )
    composition_qe = (
        _rational_chain_projection(graphs, components, inputs, middle, z)
        if rational
        else quantifier_eliminate(
            sp.And(*graphs).xreplace(renaming),
            [("exists", v) for v in bound],
            variables=(*variables, *bound),
            return_result=True,
        )
    )
    return compare_real_graphs(
        target_qe,
        composition_qe,
        input=x,
        output=z,
        parameters=parameters,
        domain=domain,
        check_functions=check_functions,
    )


def verify_decomposition(
    decomposition,
    *,
    output=None,
    domain=sp.true,
    target_selector=sp.true,
    outer_selector=sp.true,
    inner_selector=sp.true,
    check_functions=True,
):
    """Compare a two-component certificate on a selected real input domain.

    Rational certificates require a separate output coordinate when selectors
    mention it. Their intermediate coordinate is ``decomposition.parameter``.
    Implicit certificates fix all three coordinates in their polynomials.
    Selectors use only the input and output of the corresponding component.
    """
    if not isinstance(decomposition, (RationalDecomposition, AlgebraicDecomposition)):
        raise ValueError("a rational or algebraic decomposition certificate is required")
    if not decomposition.verify():
        raise ValueError("the decomposition certificate failed verification")
    if isinstance(decomposition, RationalDecomposition):
        chain = rational_chain(
            decomposition.original_expression,
            decomposition.variables,
            constants=decomposition.constants,
            initial=decomposition,
            max_depth=1,
        )
        middle = (decomposition.parameter,)
    else:
        chain = algebraic_chain(
            decomposition.relation,
            decomposition.relation.gens,
            constants=decomposition.constants,
            initial=decomposition,
            max_depth=1,
        )
        middle = (decomposition.intermediate,)
    return verify_chain(
        chain,
        output=output,
        intermediates=middle,
        domain=domain,
        target_selector=target_selector,
        component_selectors=(outer_selector, inner_selector),
        check_functions=check_functions,
    )


__all__ = [
    "verify_chain",
    "verify_decomposition",
]
