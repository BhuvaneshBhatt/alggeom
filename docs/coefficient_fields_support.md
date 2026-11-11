# Coefficient fields and mathematical scope

A coefficient parameter is a scalar in a declared field. A geometric coordinate is a generator of the polynomial ring. Specializing a parameter can change dimension or degree and must use original equations, not only a generic Groebner basis.

| Operation | QQ/number fields | K(t), characteristic zero | GF(p^d) | Positive-characteristic rational-function fields |
| --- | --- | --- | --- | --- |
| Polynomial/coordinate/module arithmetic | Supported | Supported | Supported | Arithmetic where the underlying exact domain supports it |
| Finite quotient trace/norm/minimal polynomials | Supported | Supported | Supported | Supported arithmetic; separability is not assumed |
| Finite quotient radicals/primary factors | Trace/separating methods | Supported within finite bounds | Frobenius methods | General algorithms over imperfect fields are not supplied |
| General components and affine normalization | alggeom, characteristic zero | alggeom with stated limits | Radicals in any dimension; primary factors in dimension zero; normalization unsupported | Unsupported |
| Numerical root extraction | Supported solver inputs | Specialize constants before numerical extraction | Unsupported | Unsupported |
| Real inequalities and branch checks | Declared real embeddings | Parameters belong in semialg formulas and strata | No real ordering | No real ordering |

Krull dimension measures prime chains. Finite algebra dimension is length as a vector space over the coefficient field. Geometric point count discards nilpotent multiplicity and is different from either dimension. Global finite-support module length includes residue-field degrees; local length at a field-valued point is computed over that point's field.

Construct GF(p^d) with `algroots.finite_fields.finite_field(p, d)` or `alggeom.fields.finite_field`. `sympy.GF(p**d)` is a composite-modulus residue ring when d>1 and is not that extension field. Different finite-extension presentations need explicit embeddings. Number-field embeddings also carry mathematical choices; the backend does not identify different parents merely because printed coefficients resemble one another.
