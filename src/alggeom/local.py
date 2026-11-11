"""Tangent cones, local dimensions, multiplicities and isolated intersection lengths."""

from itertools import combinations

import sympy as sp
from algroots.domains import domain_matrix, field_element

from ._algebra import integer
from .graded import degree
from .ideals import CoordinateRing, PolynomialRing
from .properties import defining_ideal


def _affine_point(value, point):
    """Choose an affine chart and extend coefficients before a local calculation."""
    from .schemes import AffinePoint, AffineScheme, ProjectivePoint, ProjectiveScheme

    if isinstance(value, ProjectiveScheme):
        if not isinstance(point, ProjectivePoint) or point.scheme != value:
            raise ValueError("a point of this projective scheme is required")
        chart = value.affine_chart(next(i for i, c in enumerate(point.coordinates) if c != 0))
        return _affine_point(chart.scheme, chart.from_projective(point))
    scheme = value if isinstance(value, AffineScheme) else AffineScheme(defining_ideal(value))
    if not isinstance(point, AffinePoint):
        point = scheme.point(point)
    else:
        if point.scheme.variables != scheme.variables:
            raise ValueError("the point must have the same ambient coordinates")
        point = scheme.point(point.coordinates, field=point.field)
    changed = scheme.base_change(point.field)
    return changed.ideal, point.coordinates


def tangent_cone(value, point):
    """Return the homogeneous tangent-cone ideal at an exact point.

    Dilating coordinates about the point and saturating the deformation
    parameter computes initial forms of the entire ideal, including forms
    created by cancellation. Taking only initial terms of input generators
    would miss such equations. Projective points use their first nonzero chart.
    """
    ideal, values = _affine_point(value, point)
    h = sp.Dummy("deformation")
    ambient = ideal.ring
    ring = PolynomialRing((h, *ambient.variables), ambient.domain, "lex")
    substitution = {v: p + h * v for v, p in zip(ambient.variables, values, strict=True)}
    deformation = ring.ideal(
        *(g.subs(substitution, simultaneous=True) for g in ideal.generators)
    ).saturate(h)
    return ambient.ideal(*(g.subs(h, 0) for g in deformation.groebner_basis))


def local_dimension(value, point):
    """Return local dimension after extending to the point's coefficient field."""
    return tangent_cone(value, point).krull_dimension


def embedding_dimension(value, point):
    """Return the dimension of the Zariski tangent space at the point."""
    from .maps import jacobian

    ideal, values = _affine_point(value, point)
    matrix = jacobian(ideal).subs(
        dict(zip(ideal.ring.variables, values, strict=True)), simultaneous=True
    )
    return len(ideal.ring.variables) - domain_matrix(matrix, ideal.ring.domain).rank()


def multiplicity(value, point):
    """Return Hilbert-Samuel multiplicity from the graded tangent cone."""
    return degree(tangent_cone(value, point))


def module_length(module, point=None, *, max_dimension=256, max_vector_dimension=4096):
    """Return the length of a finite-support module, globally or at an exact point.

    Global length is vector-space dimension over the coefficient field.
    At a field-valued point the module is base changed first and its primary
    local factor is selected. Nonzero-dimensional support is rejected.
    """
    from .components import primary_decomposition
    from .modules import PresentedModule, annihilator

    if not isinstance(module, PresentedModule):
        raise TypeError("a PresentedModule is required")
    module = module.restrict_scalars()
    if module.is_zero:
        return 0
    if point is not None:
        support, values = _affine_point(module.ring.ideal(), point)
        ambient = support.ring
        module = PresentedModule(ambient, module.presentation)
    else:
        ambient = module.ring
    ideal = annihilator(module)
    if ideal.is_unit:
        return 0
    if ideal.krull_dimension != 0:
        raise ValueError("module length requires zero-dimensional support")
    if point is not None:
        substitution = dict(zip(ambient.variables, values, strict=True))
        factors = [
            q
            for q in primary_decomposition(ideal)
            if all(
                field_element(g.subs(substitution, simultaneous=True), ambient.domain)
                == ambient.domain.zero
                for g in q.generators
            )
        ]
        if not factors:
            return 0
        ideal = factors[0]
    algebra = CoordinateRing(ideal).finite_algebra(max_dimension=max_dimension)
    size = algebra.dimension * module.generator_count
    if size > integer(max_vector_dimension, "max_vector_dimension", 1):
        raise ValueError("module vector-space dimension bound exceeded max_vector_dimension")
    matrix = sp.zeros(size, algebra.dimension * module.relation_count)
    for i in range(module.generator_count):
        for j in range(module.relation_count):
            matrix[
                i * algebra.dimension : (i + 1) * algebra.dimension,
                j * algebra.dimension : (j + 1) * algebra.dimension,
            ] = algebra.multiplication_matrix(module.presentation[i, j])
    return int(size - domain_matrix(matrix, ambient.domain).rank())


def intersection_multiplicity(left, right, point, *, max_dimension=256):
    """Compute a proper isolated local intersection using a Koszul Tor complex.

    One input must be a complete intersection in its supplied equations.
    Higher Tor lengths are included, so the answer is not assumed to equal
    the length of the sum ideal. General intersections needing a non-Koszul
    free resolution are outside this interface.
    """
    from .modules import ModuleHomomorphism, PresentedModule

    first, values = _affine_point(left, point)
    original = defining_ideal(right)
    if original.ring.variables != first.ring.variables:
        raise ValueError("the inputs must have the same ambient coordinates")
    second = original.change_ring(first.ring)
    from .schemes import AffineScheme

    AffineScheme(second).point(values)
    total = first + second
    if total.krull_dimension > 0:
        # Components not passing through the point become units in its local ring.
        # Their removal leaves the local Tor calculation unchanged.
        from .components import primary_decomposition

        substitution = dict(zip(first.ring.variables, values, strict=True))

        def local_components(ideal):
            factors = primary_decomposition(ideal)
            selected = [
                component
                for component in factors
                if all(
                    field_element(g.subs(substitution, simultaneous=True), first.ring.domain)
                    == first.ring.domain.zero
                    for g in component.generators
                )
            ]
            if not selected:
                raise ValueError("the intersection point is not in the local support")
            local_ideal = first.ring.ideal(1)
            for component in selected:
                local_ideal = local_ideal.intersection(component)
            return local_ideal

        first, second = local_components(first), local_components(second)
        total = first + second
    if total.krull_dimension != 0 or local_dimension(first, values) + local_dimension(
        second, values
    ) != len(first.ring.variables):
        raise ValueError("the intersection must be proper and finite at the point")
    if len(first.generators) != len(first.ring.variables) - first.krull_dimension:
        first, second = second, first
    if len(first.generators) != len(first.ring.variables) - first.krull_dimension:
        raise NotImplementedError("one input must have a complete-intersection presentation")
    ring = first.ring
    equations = first.generators
    target = PresentedModule(ring, sp.ImmutableMatrix(1, len(second.generators), second.generators))
    ranks = [len(tuple(combinations(range(len(equations)), i))) for i in range(len(equations) + 1)]
    terms = [
        PresentedModule(ring, sp.kronecker_product(sp.eye(rank), target.presentation))
        for rank in ranks
    ]
    differentials = []
    for i in range(1, len(terms)):
        rows = list(combinations(range(len(equations)), i - 1))
        columns = list(combinations(range(len(equations)), i))
        matrix = sp.zeros(len(rows), len(columns))
        for column, indices in enumerate(columns):
            for index, variable in enumerate(indices):
                row = rows.index(indices[:index] + indices[index + 1 :])
                matrix[row, column] = (-1) ** index * equations[variable]
        differentials.append(ModuleHomomorphism(terms[i], terms[i - 1], matrix))
    result = 0
    zero = PresentedModule(ring, sp.zeros(0, 0))
    for i, term in enumerate(terms):
        outgoing = (
            differentials[i - 1]
            if i
            else ModuleHomomorphism(term, zero, sp.zeros(0, term.generator_count))
        )
        kernel = outgoing.kernel
        if i < len(differentials):
            incoming = differentials[i].matrix
            columns = [
                outgoing.kernel_map.lift(term.element(incoming[:, j])).vector
                for j in range(incoming.cols)
            ]
            image = sp.ImmutableMatrix(
                kernel.generator_count,
                len(columns),
                [columns[c][r] for r in range(kernel.generator_count) for c in range(len(columns))],
            )
        else:
            image = sp.zeros(kernel.generator_count, 0)
        homology = PresentedModule(ring, kernel.presentation.row_join(image))
        result += (-1) ** i * module_length(homology, values, max_dimension=max_dimension)
    return result


__all__ = [
    "tangent_cone",
    "local_dimension",
    "embedding_dimension",
    "multiplicity",
    "module_length",
    "intersection_multiplicity",
]
