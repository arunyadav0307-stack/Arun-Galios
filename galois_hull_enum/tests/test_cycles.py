"""Unit tests for core/cycles.py -- the frozen f^{#} map and Algorithm E1."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.field import FiniteField  # noqa: E402
from core import poly as P  # noqa: E402
from core.factor import factor_xn_minus_lambda, is_irreducible  # noqa: E402
from core.cycles import (  # noqa: E402
    cycle_data, e1_cycle_data, hash_map, hash_map_monic,
)


def _all_monic(F, m):
    if m == 0:
        yield [1]
        return
    for code in range(F.q ** m):
        coef = []
        c = code
        for _ in range(m):
            coef.append(c % F.q)
            c //= F.q
        yield coef + [1]


# --------------------------------------------------------------------------
# the f^{#} map itself
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e", [(2, 1), (2, 4), (3, 2), (3, 3), (5, 2)])
@pytest.mark.parametrize("k", [0, 1])
def test_hash_map_matches_the_frozen_definition(p, e, k):
    """f^{#}(x) = sum_i f_0^{-p^j} f_i^{p^j} x^{m-i},  j = e-k.

    Checked against a literal, independent transcription of that formula.
    """
    if k > e - 1:
        pytest.skip("k out of range")
    F = FiniteField(p, e)
    j = e - k
    count = 0
    for f in _all_monic(F, 3):
        if f[0] == 0:                      # f^{#} requires f(0) != 0
            continue
        m = P.deg(f)
        c = F.inv(F.frobenius(f[0], j))
        expect = P.trim([F.mul(c, F.frobenius(f[i], j)) for i in range(m, -1, -1)])
        got = hash_map(F, f, j)
        assert got == expect, f"f^{{#}} mismatch at {P.format_poly(F, f)}"
        count += 1
        if count >= 60:
            break
    assert count > 0


@pytest.mark.parametrize("p,e", [(2, 3), (3, 2), (3, 4), (5, 2)])
@pytest.mark.parametrize("k", [0, 1, 2])
def test_hash_map_preserves_degree_and_is_monic(p, e, k):
    if k > e - 1:
        pytest.skip("k out of range")
    F = FiniteField(p, e)
    j = e - k
    for f in _all_monic(F, 3):
        if f[0] == 0:
            continue
        g = hash_map(F, f, j)
        assert P.deg(g) == P.deg(f), "f^{#} must preserve degree"
        assert P.plead(F, g) == 1, "f^{#} must be monic (leading coeff is f_0^{-p^j} f_0^{p^j} = 1)"
        assert g[0] != 0, "f^{#}(0) must be nonzero so the map can be iterated"


@pytest.mark.parametrize("p,e", [(2, 4), (3, 3), (5, 2)])
def test_hash_map_commutes_with_frobenius(p, e):
    """#(sigma(f)) = sigma(#(f)) -- both use powers of the same sigma."""
    F = FiniteField(p, e)
    for k in range(e):
        j = e - k
        for f in _all_monic(F, 2):
            if f[0] == 0:
                continue
            lhs = hash_map(F, P.frobenius_coeffs(F, f, 1), j)
            rhs = P.frobenius_coeffs(F, hash_map(F, f, j), 1)
            assert lhs == rhs


@pytest.mark.parametrize("p,e", [(2, 4), (3, 2), (3, 3), (3, 4), (5, 2)])
def test_k_zero_gives_the_classical_reciprocal_hence_cycles_of_length_at_most_two(p, e):
    """For k = 0 we have j = e, sigma^e = id, so
        f^{#} = f_0^{-1} x^{m} f(1/x) = f^{*},
    the classical reciprocal, and (f^{*})^{*} = f. Hence for k = 0 every cycle
    has length <= 2.

    THIS IS THE STRUCTURAL REASON the prior art (PA-1, whose maps f -> f^{*} and
    f -> f^{dagger} are involutions) cannot reach cycle length >= 3
    (BLUEPRINT section 8.1, Class C1). Verified by computation, not asserted.
    """
    F = FiniteField(p, e)
    # the reciprocal really is an involution
    for f in _all_monic(F, 2):
        if f[0] == 0:
            continue
        star1 = hash_map(F, f, e)
        star2 = hash_map(F, star1, e)
        assert star2 == f, "for k=0 the #-map must be an involution"

    # and on actual instance data, no cycle longer than 2 occurs for k=0
    for n in (5, 7, 9, 15, 21):
        for lam in (F.one, F.gen):
            cd = e1_cycle_data(F, n, lam, 0)
            if not cd.permutation_ok:
                continue
            assert max(cd.cycle_lengths) <= 2, (
                f"k=0 produced a cycle of length {max(cd.cycle_lengths)} > 2"
            )


def test_hash_map_rejects_bad_input():
    F = FiniteField(3, 2)
    with pytest.raises(ValueError):
        hash_map(F, [], 1)              # zero polynomial
    with pytest.raises(ValueError):
        hash_map(F, [0, 1], 1)          # f(0) == 0


# --------------------------------------------------------------------------
# cycle_data / Algorithm E1 structural identities
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e,n,k", [
    (3, 2, 3, 1), (5, 2, 7, 1), (3, 3, 15, 2), (3, 4, 21, 2),
    (3, 4, 5, 3), (3, 4, 17, 3), (3, 4, 25, 3),
    (2, 4, 15, 1), (2, 4, 15, 2), (2, 4, 15, 3),
    (5, 2, 3, 1), (7, 2, 5, 1), (2, 8, 7, 3),
])
def test_e1_structural_identities(p, e, n, k):
    """The identities E1 must satisfy (BLUEPRINT 5.1 / 3.1)."""
    if k > e - 1:
        pytest.skip("k out of range")
    F = FiniteField(p, e)
    for lam in [F.one, F.gen, F.gen_pow((F.q - 1) // 2 if (F.q - 1) % 2 == 0 else 1)]:
        cd = e1_cycle_data(F, n, lam, k)
        # n = n' P with gcd(n', p) = 1
        assert cd.n_prime * (p ** cd.nu) == n
        assert cd.n_prime % p != 0
        assert cd.P == p ** cd.nu
        # mu^P = lambda
        assert F.powi(cd.mu, cd.P) == lam
        # every factor irreducible and monic
        for f, m in cd.factors:
            assert is_irreducible(F, f)
            assert P.plead(F, f) == 1
            assert m == cd.P
        if not cd.permutation_ok:
            continue
        # sum of cycle lengths == number of distinct irreducible factors
        assert cd.sum_a == cd.distinct_irreducible_factors, (
            "sum_c a(c) must equal the number of distinct irreducible factors"
        )
        # cycles partition the factors
        all_idx = sorted(i for c in cd.cycles for i in c["orbit"])
        assert all_idx == list(range(cd.distinct_irreducible_factors))
        # every cycle has a well-defined degree, shared by all its members
        for c in cd.cycles:
            assert c["a"] >= 1
            assert c["d"] == P.deg(c["representative"])
        # predicted code count is (P+1)^{sum a}
        assert cd.predicted_code_count == (cd.P + 1) ** cd.sum_a


def test_cycle_data_detects_non_permutation():
    """When the standing hypothesis (H) is violated the #-map need not permute
    the irreducible factors; cycle_data must report permutation_ok = False
    rather than silently producing a bogus cycle decomposition."""
    F = FiniteField(3, 4)
    n_prime, nu, Pw, mu, fac = factor_xn_minus_lambda(F, 5, F.gen)
    distinct = [f for f, _ in fac]
    # j = 3 (k = 1); with this lambda (H) will generally fail
    cycles, ok = cycle_data(F, distinct, 3)
    if ok:
        # if it did permute, the structural identities must still hold
        assert sum(c["a"] for c in cycles) == len(distinct)
    else:
        assert cycles == []


# --------------------------------------------------------------------------
# Agreement with the [PROBE-DONE] reference in BLUEPRINT.md Annex A.3
#
# These expected values were produced by the throwaway prototype during
# blueprinting. They are [PROBE-DONE], NOT paper results. They are used here
# as an INDEPENDENTLY RECORDED cross-check of the new implementation -- the
# point being that a from-scratch re-implementation lands on the same cycle
# data. They are not the P1 acceptance gate (that is the CP1 quartet plus the
# regression sweep).
# --------------------------------------------------------------------------
A3_REFERENCE = [
    # (p, e, k, n, lambda_spec, expected_B, expected_shapes, expected_code_count)
    (3, 2, 1, 3,  "prime2",  1, [(1, 1)],          4),
    (5, 2, 1, 7,  "g^8",     3, [(1, 1), (1, 3), (1, 3)], 8),
    (3, 3, 2, 15, "prime2",  2, [(1, 1), (1, 4)],  16),
    (3, 4, 2, 21, "g^8",     2, [(1, 1), (2, 3)],  64),
    (3, 4, 3, 5,  "1",       2, [(1, 1), (4, 1)],  32),
    (3, 4, 3, 17, "1",       2, [(1, 1), (4, 4)],  32),
    (3, 4, 3, 25, "1",       3, [(1, 1), (4, 1), (4, 5)], 512),
]


@pytest.mark.parametrize("p,e,k,n,spec,expB,exp_shapes,exp_count", A3_REFERENCE)
def test_agrees_with_blueprint_A3_probe_reference(p, e, k, n, spec, expB, exp_shapes, exp_count):
    F = FiniteField(p, e)
    if spec == "prime2":
        lam = F.from_int(2)
    elif spec == "1":
        lam = F.one
    elif spec == "g^8":
        lam = F.gen_pow(8)
    else:
        raise ValueError(spec)
    cd = e1_cycle_data(F, n, lam, k)
    assert cd.permutation_ok
    assert cd.B == expB, f"B: got {cd.B}, A.3 says {expB}"
    assert cd.cycle_shapes == exp_shapes, (
        f"cycle shapes: got {cd.cycle_shapes}, A.3 says {exp_shapes}"
    )
    assert cd.predicted_code_count == exp_count, (
        f"code count: got {cd.predicted_code_count}, A.3 says {exp_count}"
    )


def test_a_ge_three_cycles_are_reachable():
    """Class C1/C2 concern cycle length a >= 3, which PA-1 cannot express.
    Verify computationally that the implementation produces them."""
    F = FiniteField(3, 4)
    found = []
    for n, lam in ((5, F.one), (17, F.one), (25, F.one)):
        cd = e1_cycle_data(F, n, lam, 3)
        found.extend(cd.cycle_lengths)
    assert max(found) >= 3, f"no cycle of length >= 3 found; got {found}"
