"""Replayable proof certificates independently of the candidate backend."""

from dataclasses import replace

import pytest
import sympy as sp
from sympy.polys.rings import PolyRing

from alggeom.groebner_certificates import (
    construct_membership_witnesses,
    verify_reduced_groebner_basis,
)
from alggeom.ideals import PolynomialRing


def test_witness_replay(monkeypatch):
    x, y = sp.symbols("x y")
    ideal = PolynomialRing((x, y), sp.QQ, "lex").ideal(x**2 - y, x * y - 1)
    proof = ideal.groebner_basis_proof_from("sympy")
    import alggeom.groebner_certificates as certificates

    monkeypatch.setattr(
        certificates, "reference_groebner", lambda *args: pytest.fail("reference basis called")
    )
    assert proof.verify(ideal.generators, ideal.ring._poly_ring)
    assert len(proof.memberships) == len(proof.basis)
    assert len(proof.memberships[0]) == len(ideal.generators)
    bad = replace(proof, memberships=((sp.S.One,) * len(ideal.generators),) + proof.memberships[1:])
    assert not bad.verify(ideal.generators, ideal.ring._poly_ring)
    assert not replace(proof, basis=proof.basis[:-1]).verify(
        ideal.generators, ideal.ring._poly_ring
    )
    assert not replace(proof, memberships=proof.memberships[:-1]).verify(
        ideal.generators, ideal.ring._poly_ring
    )
    assert not proof.verify((x, y), ideal.ring._poly_ring)


def test_proof_generator_and_zero_unit_ideal():
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z), sp.QQ, "lex")
    for equations in [
        (),
        (x, y, z),
        (x * y - z, x * x - y, y * y - z),
        (1,),
        (x**2 / 3 + y / 5, x * y / 7 - 1),
    ]:
        ideal = ring.ideal(*equations)
        proof = ideal.groebner_basis_proof_from("sympy")
        assert proof.verify(ideal.generators, ideal.ring._poly_ring)
        assert set(proof.basis) == set(ideal.groebner_basis)
    with pytest.raises(ValueError, match="max_basis"):
        construct_membership_witnesses((x,), (x,), ring._poly_ring, max_basis=0)
    with pytest.raises((ArithmeticError, ValueError), match="basis count"):
        construct_membership_witnesses(
            (x * y - z, x * x - y, y * y - z), (z**2,), ring._poly_ring, max_basis=1
        )


def test_nonmember_rejection():
    x, y, z = sp.symbols("x y z")
    ring = PolynomialRing((x, y, z), sp.QQ, "lex")
    ideal = ring.ideal(x * y - z, x * x - y, y * y - z)
    proof = ideal.groebner_basis_proof_from("sympy")
    assert proof.verify(ideal.generators, ring._poly_ring)
    assert not verify_reduced_groebner_basis(
        (x,), (y,), PolyRing((x, y), sp.QQ, "lex"), membership_witnesses=((1,),)
    )
    with pytest.raises(ValueError, match="not in the input ideal"):
        construct_membership_witnesses((x,), (y,), PolyRing((x, y), sp.QQ, "lex"))


def test_witness_pair_budget_and_malformed_proof():
    x, y = sp.symbols("x y")
    ring = PolyRing((x, y), sp.QQ, "lex")
    with pytest.raises(ValueError, match="max_pairs"):
        construct_membership_witnesses((x,), (x,), ring, max_pairs=0)
    assert not verify_reduced_groebner_basis((x,), (x,), ring, membership_witnesses=((1 / x,),))


def test_candidate_witnesses():
    x, y = sp.symbols("x y")
    ideal = PolynomialRing((x, y), sp.QQ, "lex").ideal(x * x - y, x * y - 1)
    assert ideal.groebner_basis_from("sympy") == ideal.groebner_basis
    proof = ideal.groebner_basis_proof_from("sympy")
    assert proof.verify(ideal.generators, ideal.ring._poly_ring)
    assert not replace(proof, memberships=((sp.S.Zero, sp.S.Zero),) + proof.memberships[1:]).verify(
        ideal.generators, ideal.ring._poly_ring
    )


def test_targeted_witness_budget():
    x, y = sp.symbols("x y")
    ring = PolynomialRing((x, y), sp.QQ, "lex")._poly_ring
    generators = (x * x - y, x * y - 1)
    # Not reducible directly by the generators, but a single new Buchberger
    # polynomial suffices to express x-y**2 in the ideal.
    witness = construct_membership_witnesses(
        generators, (x - y * y,), ring, max_pairs=1, max_basis=4
    )
    reconstructed = sum(ring.from_expr(f) * a for f, a in zip(generators, witness[0], strict=True))
    assert reconstructed == ring.from_expr(x - y * y)
