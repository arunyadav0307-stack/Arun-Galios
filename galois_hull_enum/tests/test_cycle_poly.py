"""Unit tests for core/cycle_poly.py and core/hull_dim.py.

Covers the frozen objects of BLUEPRINT.md sections 3.2/3.3:
  Definition 2 (local weight), Definition 3 (cycle polynomial),
  Definition 5 (restricted), T3 (transfer matrix + corollaries),
  T4 (P = 1 closed form), T5 (LCD count), T1 (hull-dimension identity).
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import zz                                    # noqa: E402
from core import poly as P                             # noqa: E402
from core.cycle_poly import (                          # noqa: E402
    closed_form_P1, cycle_poly_bruteforce,
    cycle_poly_bruteforce_restricted, cycle_poly_trace,
    cycle_poly_trace_restricted, degree_bound, det_I_minus_tT,
    dict_to_poly, local_weight, restricted_transfer_matrix,
    t_degree, transfer_matrix, verify_recurrence,
)
from core.field import FiniteField                      # noqa: E402
from core.hull_dim import (                            # noqa: E402
    hull_dimension_poly, parity_check, polylcm,
    verify_T1_on_instance,
)


def _words(P, a):
    import itertools
    return list(itertools.product(range(P + 1), repeat=a))


# --------------------------------------------------------------------------
# Definition 2 -- the local weight
# --------------------------------------------------------------------------
@pytest.mark.parametrize("P", [1, 2, 3, 4])
@pytest.mark.parametrize("a", [1, 2, 3, 4, 5])
def test_local_weight_is_a_cyclic_nearest_neighbour_sum(P, a):
    """w_a(u) = sum_m min{u_{m-1}, P-u_m} -- check against a term-by-term
    evaluation, and that it decomposes as sum_m phi(u_{m-1}, u_m)."""
    for u in _words(P, a):
        w = local_weight(u, P)
        by_terms = sum(min(u[(m - 1) % a], P - u[m]) for m in range(a))
        assert w == by_terms
        # the interaction form: phi(x,y) = min{x, P-y}
        as_phi = sum(min(u[(m - 1) % a], P - u[m]) for m in range(a))
        assert w == as_phi


@pytest.mark.parametrize("P", [1, 2, 3, 5])
def test_local_weight_range(P):
    """w_a^{(P)}(u) in [0, floor(aP/2)] (Definition 2)."""
    for a in range(1, 6):
        vals = [local_weight(u, P) for u in _words(P, a)]
        assert min(vals) == 0
        assert max(vals) == (a * P) // 2


def test_local_weight_extremal_words_are_zero_and_all_P():
    """T5's content at the level of a single word: exactly 0^a and P^a have
    weight 0."""
    for P in (1, 2, 3):
        for a in range(1, 6):
            zero = [u for u in _words(P, a) if local_weight(u, P) == 0]
            assert zero == [tuple([0] * a), tuple([P] * a)]


def test_local_weight_matches_frozen_examples():
    """Hand-evaluated instances of Definition 2. Each is worked out term by
    term in the comment so the expected value is auditable."""
    # P = 1, u = (0,0,0): every term min{0, 1} = 0  -> 0
    assert local_weight((0, 0, 0), 1) == 0
    # P = 1, u = (1,1,1): every term min{1, 0} = 0  -> 0
    assert local_weight((1, 1, 1), 1) == 0
    # P = 1, u = (0,1,0): m=0 min{u_2=0, 1-0}=0 ; m=1 min{u_0=0, 1-1}=0 ;
    #                     m=2 min{u_1=1, 1-0}=1  -> 1
    assert local_weight((0, 1, 0), 1) == 1
    # P = 2, u = (1,1): m=0 min{u_1=1, 2-1}=1 ; m=1 min{u_0=1, 2-1}=1 -> 2
    assert local_weight((1, 1), 2) == 2
    # P = 2, u = (2,0): m=0 min{u_1=0, 2-2}=0 ; m=1 min{u_0=2, 2-0}=2 -> 2
    assert local_weight((2, 0), 2) == 2


# --------------------------------------------------------------------------
# Definition 3 / T3 -- cycle polynomial by brute force == trace
# --------------------------------------------------------------------------
@pytest.mark.parametrize("P,a", [
    (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 8), (1, 12),
    (2, 1), (2, 2), (2, 3), (2, 4), (2, 5),
    (3, 1), (3, 2), (3, 3), (3, 4),
    (4, 1), (4, 2), (4, 3), (5, 1), (5, 2), (5, 3),
])
def test_trace_equals_bruteforce_word_enumeration(P, a):
    """T3: W_a^{(P)} = tr(T_P^a), against Definition 3's direct enumeration."""
    bf = dict_to_poly(cycle_poly_bruteforce(P, a))
    tr = cycle_poly_trace(P, a)
    assert bf == tr, f"P={P} a={a}: bf={zz.zformat(bf)} trace={zz.zformat(tr)}"


@pytest.mark.parametrize("P", [0, 1, 2, 3, 4, 5])
def test_transfer_matrix_shape_and_entries(P):
    """T_P(z) = (z^{min{x,P-y}})_{0<=x,y<=P}: (P+1)x(P+1), index set 0..P."""
    T = transfer_matrix(P)
    n = P + 1
    assert len(T) == n and all(len(r) == n for r in T)
    for x in range(n):
        for y in range(n):
            e = min(x, P - y)
            expect = (1,) if e == 0 else (0,) * e + (1,)
            assert T[x][y] == expect
    # diagonal entries are z^{min{x,P-x}}
    for x in range(n):
        e = min(x, P - x)
        assert T[x][x] == ((1,) if e == 0 else (0,) * e + (1,))


@pytest.mark.parametrize("P,a", [(1, 1), (1, 3), (2, 2), (2, 3), (3, 1), (3, 2), (4, 2)])
def test_W_value_at_one_is_total_word_count(P, a):
    """W_a^{(P)}(1) = (P+1)^a -- the total number of local states."""
    assert zz.zval(cycle_poly_trace(P, a), 1) == (P + 1) ** a


@pytest.mark.parametrize("P,a", [
    (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
    (2, 1), (2, 2), (2, 3), (2, 4), (2, 5),
    (3, 1), (3, 2), (3, 3), (3, 4), (4, 1), (4, 2), (4, 3), (5, 1), (5, 2),
])
def test_degree_is_floor_aP_over_2(P, a):
    poly = dict_to_poly(cycle_poly_bruteforce(P, a))
    assert zz.zdeg(poly) == degree_bound(P, a) == (a * P) // 2


@pytest.mark.parametrize("P,a", [(1, 1), (1, 2), (2, 1), (2, 3), (3, 2), (4, 1), (5, 2)])
def test_constant_term_is_two(P, a):
    """T5: c_{a,0}^{(P)} = 2."""
    poly = dict_to_poly(cycle_poly_bruteforce(P, a))
    assert zz.zcoeff(poly, 0) == 2


@pytest.mark.parametrize("P", [0, 1, 2, 3, 4, 5])
def test_det_I_minus_tT_has_t_degree_P_plus_1(P):
    """T3 corollary: deg_t det(I - t T_P(z)) = P+1."""
    assert t_degree(P) == P + 1


@pytest.mark.parametrize("P", [0, 1, 2, 3, 4, 5])
def test_det_at_z_equals_one_is_one_minus_P_plus_1_t(P):
    """At z = 1, T_P(1) = J_{P+1}, so det(I - t T_P(1)) = 1 - (P+1) t."""
    at1 = [zz.zval(c, 1) for c in det_I_minus_tT(P)]
    while at1 and at1[-1] == 0:
        at1.pop()
    assert at1 == [1, -(P + 1)]


@pytest.mark.parametrize("P", [0, 1, 2, 3, 4, 5])
def test_order_P_plus_1_recurrence(P):
    """T3 corollary: W_a^{(P)} obeys an order-(P+1) recurrence over Z[z]."""
    ok, coeffs, order = verify_recurrence(P, 2 * (P + 1) + 6)
    assert ok and order == P + 1


# --------------------------------------------------------------------------
# T4 -- the P = 1 closed form
# --------------------------------------------------------------------------
@pytest.mark.parametrize("a", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
def test_closed_form_P1_matches_trace(a):
    cf = closed_form_P1(a)
    tr = cycle_poly_trace(1, a)
    assert cf == tr


def test_closed_form_P1_verbatim_coefficient_lists():
    """The lists quoted verbatim in BLUEPRINT.md section 2.4."""
    quoted = {1: [2], 2: [2, 2], 3: [2, 6], 4: [2, 12, 2],
              5: [2, 20, 10], 6: [2, 30, 30, 2]}
    for a, exp in quoted.items():
        assert list(cycle_poly_trace(1, a)) == exp


def test_T1_matrix_for_P1():
    """T4's derivation: T_1(z) = [[1,1],[z,1]], eigenvalues 1 +/- sqrt z."""
    T = transfer_matrix(1)
    assert T == [[(1,), (1,)], [(0, 1), (1,)]]
    # characteristic polynomial lambda^2 - 2 lambda + (1 - z)
    chi = zz.zmat_charpoly_via_traces(T)
    assert chi[0] == (1,)
    assert chi[1] == (-2,)
    assert chi[2] == zz.zsub((1,), (0, 1))       # 1 - z


# --------------------------------------------------------------------------
# Definition 5 -- the restricted cycle polynomial
# --------------------------------------------------------------------------
@pytest.mark.parametrize("P,a", [
    (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (3, 3), (3, 4),
])
def test_restricted_trace_matches_bruteforce_no_adjacent_equal(P, a):
    bf = dict_to_poly(cycle_poly_bruteforce_restricted(P, a))
    tr = cycle_poly_trace_restricted(P, a)
    assert bf == tr, f"P={P} a={a}: bf={zz.zformat(bf)} trace={zz.zformat(tr)}"


@pytest.mark.parametrize("P", [0, 1, 2, 3, 4])
def test_restricted_T_has_zero_diagonal(P):
    Tt = restricted_transfer_matrix(P)
    for i in range(P + 1):
        assert Tt[i][i] == ()
        for j in range(P + 1):
            if i != j:
                assert Tt[i][j] == transfer_matrix(P)[i][j]


@pytest.mark.parametrize("P", [0, 1, 2, 3, 4])
def test_restricted_a_equals_one_special_case(P):
    """Definition 5: tilde W_1 := W_1, NOT tr(tilde T_P) (which is 0)."""
    assert cycle_poly_trace_restricted(P, 1) == cycle_poly_trace(P, 1)
    assert zz.mat_trace(restricted_transfer_matrix(P)) == ()


def test_restricted_P1_odd_a_is_zero():
    """With only two symbols, a cycle with no cyclically adjacent equal
    exponents is impossible for odd a -- so tilde W_a^{(1)} = 0 for odd a."""
    for a in (3, 5, 7):
        assert cycle_poly_trace_restricted(1, a) == ()
        assert dict_to_poly(cycle_poly_bruteforce_restricted(1, a)) == ()


# --------------------------------------------------------------------------
# T1 -- the frozen hull-dimension identity
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e,k,n", [
    (3, 2, 1, 3), (5, 2, 1, 7), (3, 3, 2, 15), (3, 4, 2, 21),
    (3, 4, 3, 5), (3, 4, 3, 17), (3, 4, 3, 25),
    (2, 4, 1, 15), (2, 4, 2, 15), (2, 4, 3, 15),
    (2, 8, 3, 7), (7, 2, 1, 5),
])
def test_T1_hull_dimension_identity(p, e, k, n):
    """T1: dim = sum_c ord_j(q) w_a^{(P)}(u), checked against the polynomial
    route n - deg lcm(g, h^{#}) on several word assignments."""
    F = FiniteField(p, e)
    for lam in {F.one, F.from_int(2) if p > 2 else F.one, F.gen}:
        res = verify_T1_on_instance(F, n, lam, k)
        for r in res:
            if "skipped" in r:
                continue                      # hypothesis (H) fails
            assert r["agree"], (
                f"T1 mismatch at {r['instance']} words={r['words']}: "
                f"lhs={r['lhs_polynomial_route']} rhs={r['rhs_local_weight_route']}"
            )


def test_parity_check_and_lcm_are_consistent():
    """h = (x^n - lambda)/g and lcm(g, h) = g*h/gcd(g,h).

    NOTE on this instance: over F_3, 2 = -1, so x^3 - 2 = x^3 + 1 = (x+1)^3
    by the Freshman's dream. Hence g = x+1 and h = (x+1)^2 are NOT coprime,
    gcd = x+1, and lcm = (x+1)^2 = h -- not g*h. The test asserts exactly
    that, so it exercises the non-coprime branch of polylcm.
    """
    F = FiniteField(3, 2)
    n, lam = 3, F.from_int(2)
    xn = [F.neg(lam)] + [0] * (n - 1) + [1]
    g = P.poly_from_int_coeffs(F, [1, 1])          # x + 1
    h = parity_check(F, n, lam, g)
    assert P.deg(h) == 2
    assert P.trim(P.pmul(F, g, h)) == P.trim(xn)
    # non-coprime: gcd(g, h) = x + 1
    assert P.deg(P.pgcd(F, g, h)) == 1
    # lcm = g*h/gcd = (x+1)^2
    lcm = polylcm(F, g, h)
    assert P.trim(lcm) == P.trim(h)
    assert P.trim(lcm) == P.trim(P.pmul(F, g, g))
    # and a coprime case, for contrast: take a genuine irreducible factor of
    # x^5 - 1 over F_9. (x+1 does NOT divide x^7-1, since (-1)^7-1 = 1 != 0
    # in F_3 -- only even n makes x+1 a factor of x^n-1.)
    from core.factor import factor_xn_minus_lambda
    F9 = FiniteField(3, 2)
    _np, _nu, _Pw, _mu, fac = factor_xn_minus_lambda(F9, 5, F9.one)
    g5 = fac[0][0]
    h5 = parity_check(F9, 5, F9.one, g5)
    assert P.deg(P.pgcd(F9, g5, h5)) == 0           # coprime
    assert P.trim(polylcm(F9, g5, h5)) == P.trim(P.pmul(F9, g5, h5))


def test_hull_dimension_of_the_zero_code_is_zero():
    """g = 1 (the whole space) has hull dimension 0."""
    F = FiniteField(3, 2)
    assert hull_dimension_poly(F, 3, F.from_int(2), 1, [1]) == 0


# --------------------------------------------------------------------------
# PA-1 boundary
# --------------------------------------------------------------------------
@pytest.mark.parametrize("P", [1, 2, 3, 4, 5, 6, 7])
def test_W1_matches_PA1_per_cycle_object(P):
    """PA-1's per-cycle object for a self-reciprocal factor:
    sum_b ||{b, P-b}|| z^b."""
    coeffs = [len({b, P - b}) for b in range(P // 2 + 1)]
    assert list(cycle_poly_trace(P, 1)) == coeffs


@pytest.mark.parametrize("P", [1, 2, 3, 4, 5, 6, 7])
def test_W2_matches_PA1_per_cycle_object(P):
    """PA-1's per-cycle object for a reciprocal pair:
    sum_{a=0}^{P} 2^{1-floor(a/P)} (a+1) z^a."""
    coeffs = [(2 ** (1 - (a // P))) * (a + 1) for a in range(P + 1)]
    assert list(cycle_poly_trace(P, 2)) == coeffs


@pytest.mark.parametrize("P", [1, 2, 3, 4, 5, 6])
def test_frozen_pointwise_identity(P):
    """BLUEPRINT section 4.1 step 1: P - max{u_m, P-u_{m-1}} = min{P-u_m, u_{m-1}}."""
    for um in range(P + 1):
        for um1 in range(P + 1):
            assert P - max(um, P - um1) == min(P - um, um1)


def test_P1_reduction_gives_PA1_Cor12_shape():
    """At P = 1: W_1 = 2, W_2 = 2(1+z), so
    N(y) = 2^B prod_{c: a(c)=2} (1 + y^{d(c)}), which is PA-1 Cor 12's
    2^{s+t} |h(l)| with |h| from Eq. (26)."""
    assert cycle_poly_trace(1, 1) == (2,)
    assert cycle_poly_trace(1, 2) == (2, 2)


@pytest.mark.parametrize("p,e,n", [(3, 2, 5), (3, 2, 7), (5, 2, 3), (3, 3, 5), (7, 2, 5)])
def test_k_zero_forces_cycle_length_at_most_two(p, e, n):
    """PA-1's scope: for k = 0 the #-map is the classical reciprocal, an
    involution, so every cycle has length <= 2."""
    from core.cycles import e1_cycle_data
    F = FiniteField(p, e)
    for lam in (F.one, F.from_int(2) if p > 2 else F.one):
        cd = e1_cycle_data(F, n, lam, 0)
        if cd.permutation_ok:
            assert max(cd.cycle_lengths) <= 2


def test_k_three_produces_cycle_length_ge_three():
    """Outside PA-1's scope (k not in {0, e/2}) cycle lengths >= 3 occur."""
    from core.cycles import e1_cycle_data
    F = FiniteField(3, 4)
    lens = set()
    for n, lam in ((5, F.one), (17, F.one), (25, F.one)):
        cd = e1_cycle_data(F, n, lam, 3)
        lens.update(cd.cycle_lengths)
    assert max(lens) >= 3


# --------------------------------------------------------------------------
# edge cases
# --------------------------------------------------------------------------
def test_P_zero_is_degenerate_but_consistent():
    """P = 0 means a single local state; W_a^{(0)} = 1 for all a."""
    for a in (1, 2, 3, 5):
        assert cycle_poly_trace(0, a) == (1,)
        assert dict_to_poly(cycle_poly_bruteforce(0, a)) == (1,)


def test_a_one_is_the_self_reciprocal_case():
    """a = 1: W_1^{(P)} = sum_x z^{min{x,P-x}}."""
    for P in (1, 2, 3, 4, 5):
        expect = {}
        for x in range(P + 1):
            e = min(x, P - x)
            expect[e] = expect.get(e, 0) + 1
        assert dict_to_poly(cycle_poly_bruteforce(P, 1)) == dict_to_poly(expect)


def test_a_two_is_the_reciprocal_pair_case():
    """a = 2: w_2(v,w) = min{v+w, 2P-v-w} (derived from Definition 2)."""
    for P in (1, 2, 3):
        for v in range(P + 1):
            for w in range(P + 1):
                assert local_weight((v, w), P) == min(v + w, 2 * P - v - w)


def test_bruteforce_rejects_invalid_a():
    with pytest.raises(ValueError):
        cycle_poly_trace(2, 0)
    with pytest.raises(ValueError):
        cycle_poly_trace_restricted(2, 0)
