"""Unit tests for the P3 (CP3) enumerator and oracles.

Covers Algorithm E2 (enumeration/enumerator.py) and Algorithm E3
(oracle/brute.py): the substitution z <- y^d, the global product, the
independent #-map, the flat (cycle-free) oracle, the frozen E3 oracle, the
factorisation certificate, and the two-oracle agreement that caught the
stale-cache defect.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import zz
from core.cycles import e1_cycle_data, hash_map
from core.factor import factor_xn_minus_lambda
from core.field import FiniteField
from enumeration.enumerator import (cycle_factor, cycle_factor_bruteforce,
                                    enumerator_distribution,
                                    enumerator_poly, run_e2,
                                    structural_checks, substitute_z_power_y)
from oracle.brute import (certify_factorisation, hash_map_independent,
                          hull_dimension_poly, oracle_distribution_cycle,
                          oracle_distribution_flat, oracle_stats)


def _element_of_order(F, r):
    """A deterministic element of multiplicative order exactly r in F_q^*."""
    if r == 1:
        return F.one
    for a in F.nonzero_elements():
        if F.multiplicative_order(a) == r:
            return a
    raise LookupError(f"no element of order {r} in F_{F.q}")


# --------------------------------------------------------------------------
# z <- y^d
# --------------------------------------------------------------------------
def test_substitute_z_power_y():
    # 2 + 3z at d=1 -> 2 + 3y
    assert substitute_z_power_y((2, 3), 1) == (2, 3)
    # 2 + 3z at d=3 -> 2 + 3y^3
    assert substitute_z_power_y((2, 3), 3) == (2, 0, 0, 3)
    # 1 + z + z^2 at d=2 -> 1 + y^2 + y^4
    assert substitute_z_power_y((1, 1, 1), 2) == (1, 0, 1, 0, 1)
    with pytest.raises(ValueError):
        substitute_z_power_y((1, 1), -1)


def test_cycle_factor_agrees_with_bruteforce_only_route():
    for (P, a, d) in [(1, 1, 1), (1, 2, 3), (3, 1, 1), (2, 3, 2), (1, 4, 5)]:
        assert cycle_factor(a, P, d) == cycle_factor_bruteforce(a, P, d)


def test_cycle_factor_rejects_bad_input():
    with pytest.raises(ValueError):
        cycle_factor(0, 1, 1)
    with pytest.raises(ValueError):
        cycle_factor(1, 0, 1)
    with pytest.raises(ValueError):
        cycle_factor(1, 1, 0)


# --------------------------------------------------------------------------
# the independent #-map
# --------------------------------------------------------------------------
# k must lie in [0, e-1], so it is fixed per field alongside (p, e)
@pytest.mark.parametrize("p,e,n,lam_int,k", [
    (3, 2, 3, 2, 1), (3, 2, 5, 1, 0), (3, 2, 9, 1, 1),
    (3, 4, 5, 1, 3), (3, 4, 17, 1, 2), (3, 4, 25, 1, 3),
    (2, 4, 5, 1, 3), (2, 4, 15, 1, 2),
    (5, 2, 7, 1, 1), (5, 2, 13, 1, 0),
])
def test_hash_map_independent_matches_core_cycles(p, e, n, lam_int, k):
    """The oracle's own #-map must reproduce core.cycles.hash_map -- they are
    separate code paths, so this is a real cross-check, not a tautology."""
    F = FiniteField(p, e)
    lam = F.from_int(lam_int)
    try:
        _np, _nu, _P, _mu, facs = factor_xn_minus_lambda(F, n, lam)
    except Exception:                                       # noqa: BLE001
        pytest.skip("instance not admissible")
    j = F.e - k
    for (f, _m) in facs:
        assert hash_map_independent(F, f, j) == hash_map(F, f, j)


def test_hash_map_independent_is_an_involution_at_k_zero():
    """For k = 0, j = e and sigma^e = id, so (f^#)^# = f."""
    F = FiniteField(3, 2)
    for (f, _m) in factor_xn_minus_lambda(F, 3, F.from_int(2))[4]:
        once = hash_map_independent(F, f, F.e)          # j = e - k = e
        twice = hash_map_independent(F, once, F.e)
        assert zz.ztrim(twice) == zz.ztrim(f)


# --------------------------------------------------------------------------
# the factorisation certificate
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e,n,lam_int", [
    (3, 2, 3, 2), (3, 2, 5, 1), (3, 3, 15, 1), (3, 4, 5, 1),
    (2, 4, 5, 1), (5, 2, 7, 1),
])
def test_certify_factorisation_accepts_genuine_factorisations(p, e, n, lam_int):
    F = FiniteField(p, e)
    lam = F.from_int(lam_int)
    _np, _nu, P, _mu, facs = factor_xn_minus_lambda(F, n, lam)
    cert = certify_factorisation(F, n, lam, facs)
    assert cert["product_ok"] and cert["irreducible_ok"]
    assert cert["distinct_ok"] and cert["uniform_multiplicity"]
    assert cert["P"] == P
    assert cert["n_factors"] == len(facs)


def test_certify_factorisation_rejects_a_wrong_factorisation():
    F = FiniteField(3, 2)
    _np, _nu, _P, _mu, facs = factor_xn_minus_lambda(F, 3, F.from_int(2))
    # perturb one factor's multiplicity
    bad = [(f, m + 1) for (f, m) in facs]
    with pytest.raises(AssertionError):
        certify_factorisation(F, 3, F.from_int(2), bad)


# --------------------------------------------------------------------------
# E2: structural checks
# --------------------------------------------------------------------------
def test_enumerator_structural_checks_on_a_known_instance():
    """F_81, n=25, k=3, lambda=1: B=3, |C|=512, LCD=8=2^3 (P2 / A.3)."""
    F = FiniteField(3, 4)
    cd = e1_cycle_data(F, 25, F.one, 3)
    chk = structural_checks(cd.cycle_shapes, cd.P,
                            cd.predicted_code_count, cd.B)
    assert chk["n_codes_E2_N1"] == 512 == chk["n_codes_predicted"]
    assert chk["lcd_count_E2"] == 8 == chk["lcd_count_2^B"]
    assert chk["all_coeffs_nonneg_int"]
    assert chk["n_codes_match"] and chk["lcd_match"]
    # sum_c a(c) must equal the number of distinct irreducible factors
    assert chk["sum_a"] == len(factor_xn_minus_lambda(F, 25, F.one)[4]) == 9


def test_enumerator_poly_of_the_empty_cycle_set():
    """No factors at all: N(y) = 1, one code (the zero code)."""
    assert enumerator_poly([], 1) == (1,)
    assert enumerator_distribution([], 1) == {0: 1}


def test_run_e2_summary():
    F = FiniteField(3, 4)
    cd = e1_cycle_data(F, 5, F.one, 3)
    res = run_e2(cd.cycle_shapes, cd.P, cd.predicted_code_count, cd.B)
    assert res.n_codes == 32
    assert res.dist == {0: 4, 1: 24, 2: 4}
    assert "OK" in res.summary()


# --------------------------------------------------------------------------
# the two oracles agree with each other and with E2
# --------------------------------------------------------------------------
PILOT = [
    (3, 2, 3, 2, 1),      # X1   F_9  P=3
    (5, 2, 7, 1, 1),      # X2   F_25 P=1
    (3, 3, 15, 1, 2),     # X3   F_27 P=3
    (3, 4, 21, 1, 2),     # X4   F_81 P=3  a=2
    (3, 4, 5, 1, 3),      # X5   F_81 P=1  a=4
    (3, 4, 17, 1, 3),     # X6   F_81 P=1  a=4
    (2, 4, 5, 1, 3),      # X8a  F_16 P=1  a=4
    (3, 4, 25, 1, 3),     # X7   F_81 P=1  512 codes
]


@pytest.mark.parametrize("p,e,n,lam_int,k", PILOT)
def test_two_oracles_agree_with_e2(p, e, n, lam_int, k):
    F = FiniteField(p, e)
    lam = F.from_int(lam_int)
    cd = e1_cycle_data(F, n, lam, k)
    assert cd.permutation_ok
    e2 = enumerator_distribution(cd.cycle_shapes, cd.P)
    oc = oracle_distribution_cycle(F, cd, n, lam, k)
    _np, _nu, _P, _mu, facs = factor_xn_minus_lambda(F, n, lam)
    of, st = oracle_distribution_flat(F, n, lam, k, facs)
    assert e2 == oc == of
    assert st["n_codes"] == cd.predicted_code_count
    assert of.get(0, 0) == 2 ** cd.B


def test_the_two_oracles_enumerate_the_same_code_set():
    """The stale-cache regression test. oracle_flat indexes codes by a FLAT
    exponent vector while oracle_cycle indexes them by the cycle product; both
    must visit exactly the same set of generator polynomials, so their
    histograms must coincide -- including on instances with a >= 3 cycles,
    which is where the removed cache used to go wrong."""
    F = FiniteField(3, 4)
    n, k, lam = 25, 3, F.one
    cd = e1_cycle_data(F, n, lam, k)
    _np, _nu, _P, _mu, facs = factor_xn_minus_lambda(F, n, lam)
    oc = oracle_distribution_cycle(F, cd, n, lam, k)
    of, _st = oracle_distribution_flat(F, n, lam, k, facs)
    assert oc == of


def test_oracle_flat_code_count_is_the_divisor_count():
    """|C| must equal (P+1)^(number of distinct irreducible factors) -- the
    flat oracle's own count, independent of E1's sum_a."""
    F = FiniteField(3, 4)
    n, k, lam = 25, 3, F.one
    _np, _nu, P, _mu, facs = factor_xn_minus_lambda(F, n, lam)
    _hist, st = oracle_distribution_flat(F, n, lam, k, facs)
    assert st["n_codes"] == (P + 1) ** len(facs) == 512


def test_oracle_stats():
    assert oracle_stats({0: 4, 1: 24, 2: 4}) == {
        "n_codes": 32, "sum_dims": 32, "mean_exact": "32/32",
        "min_dim": 0, "max_dim": 2, "support": [0, 1, 2], "lcd_count": 4,
    }


# --------------------------------------------------------------------------
# the frozen feasibility bound
# --------------------------------------------------------------------------
def test_instances_over_the_frozen_E3_bound_are_not_oracle_feasible():
    """X10 has |C| = 279,936 > 10^5. The frozen Blueprint's validation-scope
    amendment excludes it from oracle checking, so the driver must flag it
    rather than run it. Here we only assert the bound arithmetic."""
    F = FiniteField(5, 2)
    # order 6 is X10a's r; note from_int(2) has order 4 here and fails (H)
    lam = _element_of_order(F, 6)
    cd = e1_cycle_data(F, 65, lam, 1)
    assert cd.permutation_ok
    assert cd.predicted_code_count == 279936 > 10 ** 5


# --------------------------------------------------------------------------
# hull_dimension_poly guards
# --------------------------------------------------------------------------
def test_hull_dimension_poly_rejects_a_non_divisor():
    F = FiniteField(3, 2)
    with pytest.raises(ValueError):
        hull_dimension_poly(F, 3, F.from_int(2), 1, [1, 0, 1])   # x^2+1 ∤ x^3-2


def test_hull_dimension_poly_on_the_whole_space():
    """g = 1 gives C = F_q^n, whose k-Galois dual is {0}, so
    hull = C cap C^{perp_k} = {0} and the dimension is 0 -- NOT n."""
    F = FiniteField(3, 2)
    assert hull_dimension_poly(F, 3, F.from_int(2), 1, [1]) == 0


def test_hull_dimension_poly_on_the_repetition_code():
    """g = x^n - lambda gives the repetition code, hull dim 0 (it is LCD)."""
    F = FiniteField(3, 2)
    g = [F.neg(F.from_int(2))] + [0] * 2 + [1]
    assert hull_dimension_poly(F, 3, F.from_int(2), 1, g) == 0
