"""Universal Gröbner-cover proof replay and exact stratum partition checks."""

from dataclasses import replace

import sympy as sp

from alggeom.groebner_cover import comprehensive_groebner_cover


def test_linear_cover_certificate():
    x, y, a = sp.symbols("x y a")
    cover = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,))
    proof = cover.certificate()
    assert len(proof.pieces) == 2
    assert proof.verify()
    assert cover.verify_partition()
    assert cover.verify_at((0,))
    assert cover.verify_at((2,))


def test_tampered_reduction_is_rejected():
    x, y, a = sp.symbols("x y a")
    proof = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,)).certificate()
    piece = proof.pieces[0]
    identity = piece.basis_membership[0]
    altered = replace(identity, remainder=sp.Integer(1))
    modified = replace(piece, basis_membership=(altered, *piece.basis_membership[1:]))
    assert not replace(proof, pieces=(modified, *proof.pieces[1:])).verify()


def test_nonlinear_partition_covers_intersections():
    x, y, a, b = sp.symbols("x y a b")
    cover = comprehensive_groebner_cover(((a**2 + b**2) * x, a * y), (x, y), (a, b))
    assert len(cover.pieces) >= 4
    assert cover.verify_partition()
    assert cover.verify_at((0, 0))
    assert cover.verify_at((sp.I, 1))
    assert cover.verify_at((-sp.I, 1))
    assert cover.verify_at((1, 1))
    # Deleting a nonempty stratum must invalidate the coverage proof.
    from alggeom.cover_proofs import check_parameter_partition

    assert not check_parameter_partition(cover.pieces[1:], cover.parameters)


def test_membership_remainder_tampering():
    a, x, y = sp.symbols("a x y")
    original = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,)).certificate()
    first = original.pieces[0]
    witness = first.basis_membership[0]
    # The equality still holds after changing the dividend and remainder
    # together; an ideal-membership certificate must prohibit this.
    altered = replace(witness, remainder=sp.Integer(1), dividend=witness.dividend + 1)
    forged = replace(first, basis_membership=(altered, *first.basis_membership[1:]))
    assert not forged.verify()


def test_invalid_parameter_binding_tampering():
    a, x, y = sp.symbols("a x y")
    original = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,)).certificate()
    exceptional = next(piece for piece in original.pieces if piece.parameter_bindings)
    # The identities hold at a=0, but not on the claimed punctured line.
    false_stratum = replace(exceptional, condition=sp.Ne(a, 0))
    assert not false_stratum.verify()


def test_missing_denominator_guard_tampering():
    a, x, y = sp.symbols("a x y")
    original = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,)).certificate()
    generic = next(piece for piece in original.pieces if not piece.parameter_bindings)
    assert generic.nonzero_guards
    unguarded = replace(generic, nonzero_guards=())
    assert not unguarded.verify()


def test_nonlinear_partition_algebraic_intersections():
    from types import SimpleNamespace

    from alggeom.cover_proofs import check_parameter_partition

    a, b = sp.symbols("a b")
    pieces = (
        SimpleNamespace(condition=sp.Ne(a**2 + b**2, 0)),
        SimpleNamespace(condition=sp.And(sp.Eq(a - sp.I * b, 0), sp.Ne(b, 0))),
        SimpleNamespace(condition=sp.And(sp.Eq(a + sp.I * b, 0), sp.Ne(b, 0))),
        SimpleNamespace(condition=sp.And(sp.Eq(a, 0), sp.Eq(b, 0))),
    )
    assert check_parameter_partition(pieces, (a, b))
    assert not check_parameter_partition(pieces[:-1], (a, b))
    # Removing a guard causes two branches to overlap at the origin.
    overlapping = (*pieces[:1], SimpleNamespace(condition=sp.Eq(a - sp.I * b, 0)), *pieces[2:])
    assert not check_parameter_partition(overlapping, (a, b))


def test_nonlinear_full_cover_certificate():
    """Replay every stratum, including both complex branches and their intersection."""
    a, b, x, y = sp.symbols("a b x y")
    cover = comprehensive_groebner_cover(((a**2 + b**2) * x, a * y), (x, y), (a, b))
    proof = cover.certificate()
    assert len(proof.pieces) == 7
    assert proof.verify()
    assert any(piece.parameter_bindings == ((a, 0), (b, 0)) for piece in proof.pieces)
    assert all(piece.verify() for piece in proof.pieces)


def test_s_polynomial_witness_tampering():
    a, x, y = sp.symbols("a x y")
    proof = comprehensive_groebner_cover((a * x, x + y), (x, y), (a,)).certificate()
    piece = next(part for part in proof.pieces if part.critical_pairs)
    pair = piece.critical_pairs[0]
    altered = replace(pair, quotients=(pair.quotients[0] + 1, *pair.quotients[1:]))
    assert not replace(piece, critical_pairs=(altered, *piece.critical_pairs[1:])).verify()
