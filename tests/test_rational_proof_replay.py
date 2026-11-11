"""Rational certificate replay using shared algroots/FLINT arithmetic."""

import sympy as sp
from sympy.polys.rings import PolyRing

from alggeom.groebner_certificates import GroebnerBasisProof, verify_reduced_groebner_basis


def test_exact_rational_replay():
    x, y = sp.symbols("x y")
    ring = PolyRing((x, y), sp.QQ, "lex")
    generators = tuple(ring.from_expr(f) for f in (x**2 - y, y**2 - 2))
    proof = GroebnerBasisProof(generators, ((ring.one, ring.zero), (ring.zero, ring.one)))
    assert proof.verify(generators, ring, normal_form_backend="flint")
    assert verify_reduced_groebner_basis(
        generators,
        generators,
        ring,
        membership_witnesses=proof.memberships,
        normal_form_backend="sympy",
    )
    bad = (ring.from_expr(x**2 - y + 1), generators[1])
    assert not verify_reduced_groebner_basis(
        generators, bad, ring, membership_witnesses=proof.memberships, normal_form_backend="flint"
    )


def test_sparse_witness_validation():
    x, y = sp.symbols("x y")
    ring = PolyRing((x, y), sp.QQ, "lex")
    generators = tuple(ring.from_expr(f) for f in (x**2 - y, y**2 - 2))
    rows = (({(0, 0): sp.S.One}, {}), ({}, {(0, 0): sp.S.One}))
    assert verify_reduced_groebner_basis(
        generators, generators, ring, membership_witnesses=rows, normal_form_backend="flint"
    )
    assert not verify_reduced_groebner_basis(
        generators,
        generators,
        ring,
        membership_witnesses=((({(-1, 0): sp.S.One}), {}), rows[1]),
        normal_form_backend="flint",
    )
    assert not verify_reduced_groebner_basis(
        generators,
        generators,
        ring,
        membership_witnesses=((({(0, 0): sp.S(2)}), {}), rows[1]),
        normal_form_backend="flint",
    )
