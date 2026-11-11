"""Independent rational-specialization checks through Singular's standard bases."""

import os
import shutil
import subprocess

import pytest
import sympy as sp

from alggeom.groebner_cover import comprehensive_groebner_cover


def _singular_polynomial(expression, variables):
    polynomial = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
    return str(polynomial.as_expr()).replace("**", "^")


def _compare_ideals(first, second, variables):
    executable = os.environ.get("SINGULAR_BIN") or shutil.which("Singular")
    if not executable:
        pytest.skip("Singular executable is unavailable")
    ring = ",".join(str(variable) for variable in variables)
    one = ",".join(_singular_polynomial(f, variables) for f in first) or "0"
    two = ",".join(_singular_polynomial(f, variables) for f in second) or "0"
    script = f"""ring r=0,({ring}),lp;
ideal initial={one};
ideal candidate={two};
ideal initial_basis=std(initial);
ideal candidate_basis=std(candidate);
int equal=1;
for(int i=1;i<=size(initial);i++){{if(reduce(initial[i],candidate_basis)!=0){{equal=0;}}}}
for(int j=1;j<=size(candidate);j++){{if(reduce(candidate[j],initial_basis)!=0){{equal=0;}}}}
print("ALGGEOM_EQUAL="+string(equal));
quit;
"""
    completed = subprocess.run(
        [executable, "-q", "-t"],
        input=script,
        capture_output=True,
        text=True,
        timeout=25,
        check=False,
    )
    output = completed.stdout + completed.stderr
    assert completed.returncode == 0 and "ALGGEOM_EQUAL=1" in output, output


@pytest.mark.oracle
def test_linear_exception_oracle():
    a, x, y = sp.symbols("a x y")
    equations = (a * x, x + y)
    cover = comprehensive_groebner_cover(equations, (x, y), (a,))
    for value in (0, 1, -3):
        specialized = tuple(sp.expand(eq.subs(a, value)) for eq in equations)
        _compare_ideals(specialized, cover.specialize((value,)), (x, y))


@pytest.mark.oracle
def test_intersecting_exceptions_oracle():
    a, b, x, y = sp.symbols("a b x y")
    equations = ((a**2 + b**2) * x, a * y)
    cover = comprehensive_groebner_cover(equations, (x, y), (a, b))
    for point in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 3)):
        substitutions = dict(zip((a, b), point, strict=True))
        specialized = tuple(sp.expand(eq.subs(substitutions)) for eq in equations)
        _compare_ideals(specialized, cover.specialize(point), (x, y))


@pytest.mark.oracle
@pytest.mark.parametrize("degree", [3, 4])
def test_implicit_high_degree_oracle(degree):
    """Compare every rational branch type with Singular, not just a generic point."""
    from alggeom.parameter_cover import comprehensive_prime_cover

    a, b, x, y = sp.symbols("a b x y")
    exceptional = a**degree + b**3 + 1
    equations = (exceptional * x, a * y)
    cover = comprehensive_prime_cover(equations, (x, y), (a, b))
    assert cover.verify()
    for point in ((0, 0), (0, -1), (-1, 0), (1, 0), (2, 1)):
        substitutions = dict(zip((a, b), point, strict=True))
        specialized = tuple(sp.expand(f.subs(substitutions)) for f in equations)
        _compare_ideals(specialized, cover.specialize(point), (x, y))


@pytest.mark.oracle
@pytest.mark.parametrize("degree", [3, 4])
def test_implicit_quotient_basis_oracle(degree):
    """Verify the generic quotient basis at a point satisfying its coefficient guards."""
    from alggeom.ideals import PolynomialRing
    from alggeom.parameter_fields import groebner_over_prime_parameters

    a, b, x, y = sp.symbols("a b x y")
    prime = PolynomialRing((a, b), sp.QQ).ideal(a**degree + b**3 + 1)
    equations = (a * x + b * y, y - 1)
    proof = groebner_over_prime_parameters(equations, (x, y), prime)
    assert proof.verify()
    point = {a: -1, b: 0}
    assert all(sp.simplify(guard.subs(point)) != 0 for guard in proof.denominator_guards)
    specialized = tuple(sp.expand(f.subs(point)) for f in equations)
    basis = tuple(sp.cancel(g.subs(point)) for g in proof.basis)
    _compare_ideals(specialized, basis, (x, y))
