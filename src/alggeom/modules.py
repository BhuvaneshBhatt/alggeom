"""Finite module presentations, polynomial syzygies, and Kahler differentials."""

from dataclasses import dataclass
from functools import cached_property
from itertools import combinations

import sympy as sp
from algroots.domains import polynomial
from sympy.polys.distributedmodules import sdm_nf_buchberger_reduced

from ._algebra import integer
from .ideals import CoordinateRing, Ideal, PolynomialRing


def _base(ring):
    if not isinstance(ring, (PolynomialRing, CoordinateRing)):
        raise TypeError("a polynomial ring or coordinate ring is required")
    return ring


def _ambient(ring):
    return ring.ambient if isinstance(ring, CoordinateRing) else ring


def _convert(ring, expression):
    return (
        ring.normal_form(expression)
        if isinstance(ring, CoordinateRing)
        else ring.convert(expression)
    )


def _matrix(ring, values):
    return sp.ImmutableMatrix(values).applyfunc(lambda value: _convert(ring, value))


def _vector(ring, values, rank):
    vector = _matrix(ring, values)
    if vector.shape != (rank, 1):
        if rank == 0 and vector.shape == (0, 0):
            return sp.ImmutableMatrix(0, 1, [])
        raise ValueError(f"a column vector with {rank} entries is required")
    return vector


def _module_data(ring, matrix):
    """Lift columns in (S/I)^r to S^r, adjoining I times each basis vector.

    Membership in the augmented S-module is membership modulo I. Projecting
    its syzygies onto the original columns generates the kernel over S/I.
    A coefficient field uses one auxiliary variable killed by a relation,
    so the same algorithm also handles polynomial rings with no coordinates.
    """
    ambient = _ambient(ring)
    auxiliary = None if ambient.variables else sp.Dummy("module_coordinate")
    variables = ambient.variables or (auxiliary,)
    old = ambient._poly_ring.domain.old_poly_ring(*variables, order=ambient.order)
    relations = ring.ideal.generators if isinstance(ring, CoordinateRing) else ()
    if auxiliary is not None:
        relations = (*relations, auxiliary)
    columns = [list(matrix[:, j]) for j in range(matrix.cols)]
    columns.extend(
        [relation if i == j else sp.S.Zero for i in range(matrix.rows)]
        for relation in relations
        for j in range(matrix.rows)
    )
    columns_in_ring = [
        [old.convert(polynomial(value, variables, old.dom).rep) for value in column]
        for column in columns
    ]
    # SymPy's module Groebner routine cannot reduce a list containing a zero
    # generator. Retain the column positions separately for witnesses/kernels.
    indices = tuple(i for i, column in enumerate(columns_in_ring) if any(column))
    module = old.free_module(matrix.rows).submodule(*(columns_in_ring[i] for i in indices))
    return module, auxiliary, indices


def _expression(ring, old, value, auxiliary):
    expression = old.to_sympy(value)
    if auxiliary is not None:
        expression = expression.subs(auxiliary, 0)
    return _convert(ring, expression)


@dataclass(frozen=True)
class Submodule:
    """The column span of a matrix in a finite free module over S or S/I.

    Columns are retained in their supplied order, including zero columns, so
    membership coefficients and syzygy coordinates refer to the input matrix.
    """

    ring: PolynomialRing | CoordinateRing
    generators: sp.ImmutableMatrix

    def __post_init__(self):
        _base(self.ring)
        object.__setattr__(self, "generators", _matrix(self.ring, self.generators))

    @property
    def ambient_rank(self):
        """The number of rows of the generator matrix, not the module's rank."""
        return self.generators.rows

    @cached_property
    def _data(self):
        return _module_data(self.ring, self.generators)

    def _module_vector(self, vector):
        module, auxiliary, _ = self._data
        variables = _ambient(self.ring).variables or (auxiliary,)
        return module.container.convert(
            [
                module.ring.convert(polynomial(value, variables, module.ring.dom).rep)
                for value in vector
            ]
        )

    def normal_form(self, vector):
        """Reduce a column modulo the submodule using its module Groebner basis."""
        vector = _vector(self.ring, vector, self.ambient_rank)
        if not self.ambient_rank:
            return vector
        module, auxiliary, _ = self._data
        reduced = module.reduce_element(self._module_vector(vector), NF=sdm_nf_buchberger_reduced)
        return sp.ImmutableMatrix(
            self.ambient_rank,
            1,
            [_expression(self.ring, module.ring, value, auxiliary) for value in reduced],
        )

    def contains(self, vector):
        """Decide whether a column belongs to this submodule."""
        return self.normal_form(vector).is_zero_matrix is True

    def coordinates(self, vector):
        """Return coefficients in the input columns, or None for a nonmember.

        Over S/I the resulting identity holds modulo I in each row. Auxiliary
        coefficients for the lifted ideal relations are omitted.
        """
        vector = _vector(self.ring, vector, self.ambient_rank)
        if not self.contains(vector):
            return None
        if not self.ambient_rank:
            return (sp.S.Zero,) * self.generators.cols
        module, auxiliary, indices = self._data
        coefficients = module.in_terms_of_generators(self._module_vector(vector))
        result = [sp.S.Zero] * self.generators.cols
        for index, value in zip(indices, coefficients, strict=True):
            if index < len(result):
                result[index] = _expression(self.ring, module.ring, value, auxiliary)
        return tuple(result)

    def syzygies(self):
        """Return a matrix whose columns generate the relations among the input columns."""
        return syzygies(self.generators, self.ring)


def syzygies(matrix, ring):
    """Return column generators of the exact kernel of a polynomial matrix.

    The matrix acts on column vectors. For a coordinate ring S/I, relations
    are computed modulo I by lifting to a free S-module with the I-relations
    included. The output need not be a minimal generating set or a basis.
    """
    ring = _base(ring)
    matrix = _matrix(ring, matrix)
    if not matrix.rows:
        return sp.ImmutableMatrix.eye(matrix.cols).applyfunc(lambda value: _convert(ring, value))
    if not matrix.cols:
        return sp.ImmutableMatrix(0, 0, [])
    module, auxiliary, indices = _module_data(ring, matrix)
    relations = module.syzygy_module()
    columns = []
    for relation in relations.gens:
        values = [sp.S.Zero] * matrix.cols
        for index, value in zip(indices, relation, strict=True):
            if index < matrix.cols:
                values[index] = _expression(ring, module.ring, value, auxiliary)
        column = tuple(values)
        if any(column) and column not in columns:
            columns.append(column)
    for index in range(matrix.cols):
        if index not in indices:
            column = tuple(_convert(ring, int(i == index)) for i in range(matrix.cols))
            if any(column):
                columns.append(column)
    return sp.ImmutableMatrix(matrix.cols, len(columns), lambda i, j: columns[j][i])


@dataclass(frozen=True)
class PresentedModule:
    """The cokernel of a matrix A^s -> A^r over a polynomial or coordinate ring.

    Rows index the r module generators and columns the s defining relations.
    Presentations are finite; they are not assumed minimal or free resolutions.
    """

    ring: PolynomialRing | CoordinateRing
    presentation: sp.ImmutableMatrix

    def __post_init__(self):
        _base(self.ring)
        object.__setattr__(self, "presentation", _matrix(self.ring, self.presentation))

    def element(self, vector):
        """Construct a residue class in this presented module."""
        if isinstance(vector, ModuleElement):
            if vector.module != self:
                raise TypeError("module elements must have the same parent")
            return vector
        return ModuleElement(self, vector)

    def __call__(self, vector):
        return self.element(vector)

    @cached_property
    def gens(self):
        """The classes of the ordered presentation generators."""
        return tuple(
            self.element(sp.eye(self.generator_count)[:, i]) for i in range(self.generator_count)
        )

    @property
    def zero(self):
        """The additive identity in this module."""
        return self.element(sp.zeros(self.generator_count, 1))

    def hom(self, matrix, target):
        """Define a linear map by columns of generator images."""
        return ModuleHomomorphism(self, target, matrix)

    def tensor_product(self, other):
        """Return the tensor product over the current base ring."""
        return tensor_product(self, other)

    def hom_module(self, target):
        """Return the module of linear maps to target."""
        return HomModule(self, target)

    def annihilator(self):
        """Return the annihilator in the ambient polynomial ring."""
        return annihilator(self)

    @property
    def generator_count(self):
        """The number of generators in this presentation."""
        return self.presentation.rows

    @property
    def relation_count(self):
        """The number of relation columns in this presentation."""
        return self.presentation.cols

    @cached_property
    def relations(self):
        """The submodule of the free module generated by the presentation columns."""
        return Submodule(self.ring, self.presentation)

    def normal_form(self, vector):
        """Return a representative of the class of a column in this cokernel."""
        return self.relations.normal_form(vector)

    @property
    def is_zero(self):
        """Decide whether the presented module is zero."""
        return all(
            self.relations.contains(sp.eye(self.generator_count)[:, i])
            for i in range(self.generator_count)
        )

    def syzygies(self):
        """Return generators for the relations among the presentation's relation columns."""
        return syzygies(self.presentation, self.ring)

    def fitting_ideal(self, index):
        """Return the indexed Fitting ideal, as an ideal in the ambient polynomial ring."""
        return fitting_ideal(self, index)

    def restrict_scalars(self):
        """Present this module over its ambient polynomial ring S.

        For a module over S/I, append I times every free generator to its
        relation matrix. This changes the base of a subsequent resolution:
        an S-resolution is finite, whereas an S/I-resolution may be infinite.
        """
        if isinstance(self.ring, PolynomialRing):
            return self
        rank = self.generator_count
        relations = sp.ImmutableMatrix(
            rank,
            rank * len(self.ring.ideal.generators),
            lambda i, j: self.ring.ideal.generators[j // rank] if i == j % rank else 0,
        )
        return PresentedModule(self.ring.ambient, self.presentation.row_join(relations))

    def free_resolution(self, *, generator_degrees=None, max_length=None):
        """Compute a minimal standard-graded resolution over a polynomial ring."""
        from .resolutions import free_resolution

        return free_resolution(self, generator_degrees=generator_degrees, max_length=max_length)


def fitting_ideal(module, index):
    """Return Fitt_index(M), generated by minors of size generator_count-index.

    For A=S/I the result is the inverse-image ideal in S and contains I.
    Minors of size <=0 generate the unit ideal; a size larger than either
    matrix dimension gives the zero ideal of A, whose inverse image is I.
    """
    if not isinstance(module, PresentedModule):
        raise TypeError("a PresentedModule is required")
    index = integer(index, "index")
    matrix = module.presentation
    size = matrix.rows - index
    ambient = _ambient(module.ring)
    defining = module.ring.ideal.generators if isinstance(module.ring, CoordinateRing) else ()
    if size <= 0:
        return ambient.ideal(1)
    minors = [
        matrix.extract(rows, columns).det()
        for rows in combinations(range(matrix.rows), size)
        for columns in combinations(range(matrix.cols), size)
    ]
    return ambient.ideal(*defining, *minors)


def _coordinate_ring(ring):
    if isinstance(ring, Ideal):
        return CoordinateRing(ring)
    if isinstance(ring, PolynomialRing):
        return ring.quotient()
    if not isinstance(ring, CoordinateRing):
        raise TypeError("an ideal, polynomial ring, or coordinate ring is required")
    return ring


def differentials(ring):
    """Present Omega_(A/K) by the transposed Jacobian of A=S/I.

    The generators are the ordered ambient coordinate differentials. Both
    defining equations and their differential relations are respected over A.
    """
    ring = _coordinate_ring(ring)
    matrix = sp.ImmutableMatrix(
        len(ring.variables),
        len(ring.ideal.generators),
        lambda i, j: sp.diff(ring.ideal.generators[j], ring.variables[i]),
    )
    return PresentedModule(ring, matrix)


def differential(expression, ring):
    """Return the coordinate column representing df in Omega_(A/K).

    Coefficient-field parameters are constants for this relative differential.
    Equal representatives in A give equal classes in the differential module.
    """
    ring = _coordinate_ring(ring)
    expression = ring.normal_form(expression)
    return sp.ImmutableMatrix(
        len(ring.variables),
        1,
        [ring.normal_form(sp.diff(expression, variable)) for variable in ring.variables],
    )


@dataclass(frozen=True, eq=False)
class ModuleElement:
    """A parent-aware residue class in a finitely presented module."""

    module: PresentedModule
    vector: sp.ImmutableMatrix

    __hash__ = None
    _op_priority = 1000

    def __post_init__(self):
        if not isinstance(self.module, PresentedModule):
            raise TypeError("a PresentedModule parent is required")
        object.__setattr__(self, "vector", self.module.normal_form(self.vector))

    def __eq__(self, other):
        return (
            isinstance(other, ModuleElement)
            and self.module == other.module
            and self.module.normal_form(self.vector - other.vector).is_zero_matrix is True
        )

    def __bool__(self):
        return self.vector.is_zero_matrix is not True

    def __add__(self, other):
        other = self.module.element(other)
        return self.module.element(self.vector + other.vector)

    def __neg__(self):
        return self.module.element(-self.vector)

    def __sub__(self, other):
        return self + (-self.module.element(other))

    def __mul__(self, scalar):
        scalar = scalar.as_expr() if hasattr(scalar, "as_expr") else scalar
        return self.module.element(self.vector * _convert(self.module.ring, scalar))

    __rmul__ = __mul__


def _compatible(source, target):
    if not isinstance(source, PresentedModule) or not isinstance(target, PresentedModule):
        raise TypeError("presented modules are required")
    if source.ring != target.ring:
        raise ValueError("the modules must have the same base ring")


def _relative_kernel(matrix, relations, ring):
    """Lift the kernel of a map into a presented module to its free source."""
    lifted = syzygies(matrix.row_join(-relations), ring)
    return sp.ImmutableMatrix(lifted[: matrix.cols, :])


@dataclass(frozen=True)
class ModuleHomomorphism:
    """A linear map of presented modules, with columns the images of source generators.

    Relations are checked in the target cokernel. Kernel and image objects
    retain presentations; kernel_map and image_map supply their natural
    inclusions. Cokernel_map is the natural target quotient.
    """

    source: PresentedModule
    target: PresentedModule
    matrix: sp.ImmutableMatrix

    def __post_init__(self):
        _compatible(self.source, self.target)
        matrix = _matrix(self.source.ring, self.matrix)
        if matrix.shape != (self.target.generator_count, self.source.generator_count):
            raise ValueError("the map needs one target column per source generator")
        if self.source.relation_count:
            relation_images = matrix * self.source.presentation
            if not all(
                self.target.relations.contains(relation_images[:, j])
                for j in range(self.source.relation_count)
            ):
                raise ValueError("the generator images do not respect the source relations")
        object.__setattr__(self, "matrix", matrix)

    def __call__(self, element):
        element = self.source.element(element)
        return self.target.element(self.matrix * element.vector)

    def compose(self, inner):
        """Return self after a map into its source module."""
        if not isinstance(inner, ModuleHomomorphism) or inner.target != self.source:
            raise ValueError("the intermediate modules must agree")
        return ModuleHomomorphism(inner.source, self.target, self.matrix * inner.matrix)

    @cached_property
    def _kernel_lifts(self):
        return _relative_kernel(self.matrix, self.target.presentation, self.source.ring)

    @cached_property
    def kernel(self):
        """The kernel as a finitely presented module over the same base."""
        relations = _relative_kernel(self._kernel_lifts, self.source.presentation, self.source.ring)
        return PresentedModule(self.source.ring, relations)

    @property
    def kernel_map(self):
        """The natural inclusion of kernel into source."""
        return ModuleHomomorphism(self.kernel, self.source, self._kernel_lifts)

    @cached_property
    def image(self):
        """The image, presented by relations among its supplied generator images."""
        return PresentedModule(self.source.ring, self._kernel_lifts)

    @property
    def image_map(self):
        """The natural inclusion of the presented image into target."""
        return ModuleHomomorphism(self.image, self.target, self.matrix)

    @cached_property
    def cokernel(self):
        """The target divided by this map's image."""
        return PresentedModule(self.source.ring, self.target.presentation.row_join(self.matrix))

    @property
    def cokernel_map(self):
        """The natural quotient from target to cokernel."""
        return ModuleHomomorphism(self.target, self.cokernel, sp.eye(self.target.generator_count))

    @property
    def is_injective(self):
        """Whether the kernel is zero."""
        return self.kernel.is_zero

    @property
    def is_surjective(self):
        """Whether the cokernel is zero."""
        return self.cokernel.is_zero

    def lift(self, element):
        """Return a source preimage of a target element, raising for a nonmember."""
        element = self.target.element(element)
        span = Submodule(self.source.ring, self.matrix.row_join(self.target.presentation))
        coordinates = span.coordinates(element.vector)
        if coordinates is None:
            raise ValueError("the element is not in the map's image")
        return self.source.element(coordinates[: self.source.generator_count])

    def inverse(self):
        """Return the inverse of an isomorphism; nonisomorphisms raise."""
        if not self.is_injective or not self.is_surjective:
            raise ValueError("the module map is not an isomorphism")
        columns = [self.lift(e).vector for e in self.target.gens]
        matrix = sp.ImmutableMatrix(
            self.source.generator_count, self.target.generator_count, lambda i, j: columns[j][i]
        )
        return ModuleHomomorphism(self.target, self.source, matrix)


def direct_sum(*modules):
    """Return the direct sum of one or more modules over the same ring."""
    if not modules:
        raise ValueError("at least one module is required")
    for module in modules:
        _compatible(modules[0], module)
    return PresentedModule(modules[0].ring, sp.diag(*(m.presentation for m in modules)))


def _kronecker(left, right):
    """Form an explicit tensor matrix, preserving empty row and column shapes."""
    return sp.ImmutableMatrix(
        left.rows * right.rows,
        left.cols * right.cols,
        lambda i, j: left[i // right.rows, j // right.cols] * right[i % right.rows, j % right.cols],
    )


def tensor_product(left, right):
    """Present a tensor product with basis order e_i tensor f_j, j varying fastest."""
    _compatible(left, right)
    relations = _kronecker(left.presentation, sp.eye(right.generator_count)).row_join(
        _kronecker(sp.eye(left.generator_count), right.presentation)
    )
    return PresentedModule(left.ring, relations)


def annihilator(module):
    """Return the annihilator, lifted to the ambient polynomial ring for a quotient base."""
    if not isinstance(module, PresentedModule):
        raise TypeError("a PresentedModule is required")
    ambient = _ambient(module.ring)
    defining = module.ring.ideal.generators if isinstance(module.ring, CoordinateRing) else ()
    result = ambient.ideal(1)
    for i in range(module.generator_count):
        basis = sp.eye(module.generator_count)[:, i]
        relations = syzygies(module.presentation.row_join(-basis), module.ring)
        result = result.intersection(ambient.ideal(*defining, *relations[-1, :]))
    return result


@dataclass(frozen=True, init=False)
class HomModule(PresentedModule):
    """The module Hom(M,N), with conversions between elements and actual linear maps."""

    source: PresentedModule
    target: PresentedModule

    def __init__(self, source, target):
        _compatible(source, target)
        rows, columns = target.generator_count, source.generator_count
        # Column-major vectorization turns Q*P_M into (P_M^T tensor I)*vec(Q).
        condition = _kronecker(source.presentation.T, sp.eye(rows))
        target_relations = _kronecker(sp.eye(source.relation_count), target.presentation)
        lifts = _relative_kernel(condition, target_relations, source.ring)
        zero_maps = _kronecker(sp.eye(columns), target.presentation)
        relations = _relative_kernel(lifts, zero_maps, source.ring)
        PresentedModule.__init__(self, source.ring, relations)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "_lifts", lifts)
        object.__setattr__(self, "_zero_maps", zero_maps)

    def homomorphism(self, element):
        """Interpret a Hom-module element as a map M -> N."""
        vector = self._lifts * self.element(element).vector
        rows, columns = self.target.generator_count, self.source.generator_count
        matrix = sp.ImmutableMatrix(rows, columns, lambda i, j: vector[i + j * rows])
        return ModuleHomomorphism(self.source, self.target, matrix)

    def from_homomorphism(self, mapping):
        """Express a map in this Hom presentation, modulo target relations."""
        if (
            not isinstance(mapping, ModuleHomomorphism)
            or mapping.source != self.source
            or mapping.target != self.target
        ):
            raise ValueError("a map with this Hom module's source and target is required")
        vector = sp.ImmutableMatrix(
            [
                mapping.matrix[i, j]
                for j in range(mapping.matrix.cols)
                for i in range(mapping.matrix.rows)
            ]
        )
        span = Submodule(self.ring, self._lifts.row_join(self._zero_maps))
        coefficients = span.coordinates(vector)
        return self.element(coefficients[: self.generator_count])


__all__ = [
    "Submodule",
    "PresentedModule",
    "ModuleElement",
    "ModuleHomomorphism",
    "HomModule",
    "syzygies",
    "fitting_ideal",
    "differentials",
    "differential",
    "direct_sum",
    "tensor_product",
    "annihilator",
]
