"""Ideal and quotient adapters to the algroots Gröbner implementation."""

import json
from pathlib import Path

import pytest
import sympy as sp
from algroots import NativeF4Limits
from algroots.native_f4 import compiled_available

from alggeom.ideals import PolynomialRing


@pytest.mark.skipif(not compiled_available(), reason="compiled kernel is not installed")
@pytest.mark.parametrize("prime", [2, 101, 1009])
def test_prime_routing(prime, monkeypatch):
    from algroots import groebner as native

    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y), domain=sp.GF(prime))
    ideal = ring.ideal(x**2 - y, y**2 - 1)

    def reject_lift(*args, **kwargs):
        raise AssertionError("prime-field request reached rational lifting")

    monkeypatch.setattr(native, "native_modular_basis", reject_lift)
    assert set(ideal.groebner_basis_from("native")) == set(ideal.groebner_basis)


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_quotient_selection(monkeypatch):
    from algroots import groebner

    x, y = sp.symbols("x y")
    calls = []
    original = groebner.compute_groebner_basis

    def record(inputs, ring, **kwargs):
        calls.append(kwargs["backend"])
        return original(inputs, ring, **kwargs)

    monkeypatch.setattr(groebner, "compute_groebner_basis", record)
    quotient = PolynomialRing((x, y), sp.GF(101)).quotient(
        x * x - y, y * y - 1, groebner_backend="native"
    )
    assert quotient.normal_form(x**4 - 1) == 0
    assert quotient.normal_form(x) == x
    assert calls == ["native"]


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_bound_coefficient_symbols(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 extension unavailable")
    x = sp.Symbol("x")
    alpha = sp.CRootOf(x**3 - x - 1, 0)
    field = sp.QQ.algebraic_field(alpha)
    ring = PolynomialRing((x,), field)
    ideal = ring.ideal(x - alpha, groebner_backend=backend)
    assert (
        ideal.groebner_basis_from(backend)
        == ring.ideal(x - alpha, groebner_backend="sympy").groebner_basis
    )
    assert ideal.normal_form(x) == alpha


@pytest.mark.parametrize("domain", [sp.QQ, sp.GF(2), sp.GF(101)])
@pytest.mark.parametrize("order", ["lex", "grevlex"])
@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_finite_conversion(domain, order, backend):
    from algroots import NativeF4Limits, groebner_backends

    if not groebner_backends()[backend].available:
        pytest.skip(f"{backend} is not installed")
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y), domain, order)
    coordinate = ring.quotient(x * x - y, y * y - 1, groebner_backend=backend)
    algebra = coordinate.finite_algebra()
    assert algebra.dimension == coordinate.vector_space_dimension == 4
    assert algebra.domain == domain
    assert algebra.normal_form(x**4 - 1) == coordinate.normal_form(x**4 - 1) == 0
    if backend == "native":
        from algroots import GroebnerLimitError

        bounded = ring.quotient(
            x * x - y,
            y * y - 1,
            groebner_backend=backend,
            groebner_limits=NativeF4Limits(max_terms=1),
        )
        with pytest.raises(GroebnerLimitError):
            bounded.finite_algebra()


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_derived_policy(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y))
    source = ring.ideal(x * y, groebner_backend=backend)
    other = ring.ideal(x, groebner_backend=backend)
    derived = [
        source + other,
        source * other,
        source.intersection(other),
        source.quotient(x),
        source.saturate(x),
        source.eliminate(y),
        source.change_ring(ring.with_order("lex")),
    ]
    expected = [(x,), (x * x * y,), (x * y,), (y,), (y,), (), (x * y,)]
    for ideal, basis in zip(derived, expected, strict=True):
        assert ideal.groebner_backend == backend
        assert ideal.groebner_basis == basis
    further = ring.quotient(x * y, groebner_backend=backend).quotient(x, groebner_backend=backend)
    assert further.ideal.groebner_backend == backend
    assert further.normal_form(x) == 0
    assert further.normal_form(y) == y


def test_factory_defaults():
    import inspect

    from alggeom.ideals import CoordinateRing, Ideal

    for constructor in (PolynomialRing.ideal, PolynomialRing.quotient, CoordinateRing.quotient):
        assert inspect.signature(constructor).parameters["groebner_backend"].default == "native"
    assert inspect.signature(Ideal).parameters["groebner_backend"].default == "native"
    assert inspect.signature(Ideal.groebner_basis_from).parameters["backend"].default == "native"


def _rational_systems():
    x, y, z, w = sp.symbols("x y z w")
    return (
        ((x, y), (x - 1, x + 1)),
        ((x, y), (x * x - y, y * y - x, 2 * x * x - 2 * y)),
        ((x, y, z), (x * x - y, y * y - z, z * z - 1)),
        (
            (x, y, z, w),
            (
                3 * w * w - 8 * z,
                w * z - 6 * y,
                z * z - 36 * x,
                w * y - 16 * x,
                z * y - 6 * w * x,
                3 * y * y - 8 * z * x,
            ),
        ),
    )


@pytest.mark.skipif(not compiled_available(), reason="native algroots kernel unavailable")
@pytest.mark.parametrize("index", range(4))
def test_rational_basis(index):
    variables, equations = _rational_systems()[index]
    ideal = PolynomialRing(variables, sp.QQ, "grevlex").ideal(*equations)
    observed = ideal.groebner_basis_from(
        "native", limits=NativeF4Limits(max_basis=512, max_reduction_steps=2000000)
    )
    reference = sp.groebner(equations, *variables, domain=sp.QQ, order="grevlex")
    assert observed == tuple(poly.as_expr() for poly in reference.polys)


@pytest.mark.skipif(not compiled_available(), reason="native algroots kernel unavailable")
@pytest.mark.parametrize("index", range(2))
def test_application_basis(index):
    case = json.loads(
        (Path(__file__).parent / "fixtures/groebner_examples.json").read_text(encoding="utf-8")
    )[index]
    variables = sp.symbols(f"x0:{case['nvars']}")
    equations = tuple(
        sum(c * sp.prod(x**degree for x, degree in zip(variables, m, strict=True)) for m, c in p)
        for p in case["polynomials"]
    )
    ideal = PolynomialRing(variables, sp.QQ, "grevlex").ideal(*equations)
    observed = ideal.groebner_basis_from(
        "native", limits=NativeF4Limits(max_basis=512, max_reduction_steps=2000000)
    )
    reference = sp.groebner(equations, *variables, domain=sp.QQ, order="grevlex")
    assert observed == tuple(poly.as_expr() for poly in reference.polys)
