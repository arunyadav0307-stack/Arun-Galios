"""Unit tests for core/factor.py.

Includes DEDICATED REGRESSION TESTS for the two defects recorded in
BLUEPRINT.md A.7.1:

  DEFECT (a)  even-characteristic absolute trace must include the i=0 term.
  DEFECT (b)  mu = sigma^{e - (nu mod e)}(lambda); sigma^{e-nu} is undefined
              (negative exponent) when nu > e, e.g. q=3, n=9.
"""

import os
import random
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.field import FiniteField  # noqa: E402
from core import poly as P  # noqa: E402
from core.factor import (  # noqa: E402
    absolute_trace_polynomial, absolute_trace_polynomial_direct, ddf, edf,
    factor_poly, factor_squarefree, factor_xn_minus_lambda, is_irreducible,
    squarefree_decomposition, v_p,
)


# --------------------------------------------------------------------------
# irreducibility test, cross-checked against a brute-force enumeration
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e,maxdeg", [(2, 1, 5), (2, 2, 3), (3, 1, 5), (3, 2, 3), (5, 1, 4)])
def test_is_irreducible_against_brute_force(p, e, maxdeg):
    """Independent oracle: enumerate ALL monic polynomials up to maxdeg and
    decide irreducibility by trial division. Compare with is_irreducible."""
    F = FiniteField(p, e)

    def brute_irreducible(f):
        m = P.deg(f)
        if m <= 0:
            return False
        if m == 1:
            return True
        # trial-divide by every monic polynomial of degree 1 .. m//2
        def divide(a, b):
            return P.pmod(F, a, b) == []
        for d in range(1, m // 2 + 1):
            for cand in _all_monic(F, d):
                if divide(f, cand):
                    return False
        return True

    checked = 0
    for m in range(1, maxdeg + 1):
        for f in _all_monic(F, m):
            assert is_irreducible(F, f) == brute_irreducible(f), (
                f"irreducibility disagreement at {P.format_poly(F, f)} over F_{F.q}"
            )
            checked += 1
    assert checked > 0


def _all_monic(F, m):
    """All monic polynomials of degree m over F (deterministic order)."""
    if m == 0:
        yield [1]
        return
    total = F.q ** m
    for code in range(total):
        coef = []
        c = code
        for _ in range(m):
            coef.append(c % F.q)
            c //= F.q
        yield coef + [1]


# --------------------------------------------------------------------------
# factorisation: product identity + irreducibility of every factor
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e", [(2, 1), (2, 3), (2, 8), (3, 1), (3, 2), (3, 4), (5, 2), (7, 2)])
def test_factor_poly_is_correct(p, e):
    F = FiniteField(p, e)
    rng = random.Random(8080 + p * 10 + e)
    for _ in range(40):
        # build a polynomial as a product of random factors, then check the
        # factoriser recovers a factorisation with the same degree profile
        target_deg = rng.randrange(1, 8)
        f = [1]
        pieces = []
        built = 0
        while built < target_deg:
            d = rng.randrange(1, min(4, target_deg - built) + 1)
            g = P.pmonic(F, [rng.randrange(F.q) for _ in range(d)] + [1])
            if P.deg(g) != d:
                continue
            f = P.pmul(F, f, g)
            pieces.append(g)
            built += d
        f = P.pmonic(F, f)
        fac = factor_poly(F, f, seed=rng.randrange(10 ** 6))
        # product identity
        prod = [1]
        for irr, m in fac:
            for _ in range(m):
                prod = P.pmul(F, prod, irr)
        assert P.trim(prod) == f, "factorisation product != original"
        # every factor irreducible and monic
        for irr, m in fac:
            assert m >= 1
            assert P.plead(F, irr) == 1
            assert is_irreducible(F, irr)
        # degree additivity
        assert sum(P.deg(i) * m for i, m in fac) == P.deg(f)


@pytest.mark.parametrize("p,e", [(2, 4), (3, 2), (5, 2), (2, 8)])
def test_squarefree_decomposition(p, e):
    F = FiniteField(p, e)
    rng = random.Random(555)
    for _ in range(25):
        a = P.pmonic(F, [rng.randrange(F.q) for _ in range(rng.randrange(1, 5))] + [1])
        b = P.pmonic(F, [rng.randrange(F.q) for _ in range(rng.randrange(1, 5))] + [1])
        # f = a * b^2  -> squarefree decomposition must find a part of mult 1
        f = P.pmul(F, a, P.pmul(F, b, b))
        dec = squarefree_decomposition(F, f)
        prod = [1]
        for g, m in dec:
            assert m >= 1
            for _ in range(m):
                prod = P.pmul(F, prod, g)
            # each g must be squarefree: gcd(g, g') == 1
            assert P.deg(P.pgcd(F, g, P.pderiv(F, g))) == 0
        assert P.trim(prod) == f


# --------------------------------------------------------------------------
# DEFECT (a) -- even characteristic absolute trace must include i = 0
# --------------------------------------------------------------------------
def test_defect_a_absolute_trace_includes_i_equals_zero():
    """The two implementations (transitive and literal direct sum) must agree,
    and BOTH must include the i=0 term. Omitting it is exactly the defect that
    made even-characteristic splitting fail."""
    F = FiniteField(2, 4)          # F_16
    # f = product of two distinct degree-2 irreducibles over F_16
    rng = random.Random(20240)
    facs = []
    for f in _all_monic(F, 2):
        if is_irreducible(F, f):
            facs.append(f)
        if len(facs) >= 2:
            break
    f = P.pmul(F, facs[0], facs[1])
    d = 2
    h = [1, 1, 0, F.gen]           # an arbitrary residue class

    t_fast = absolute_trace_polynomial(F, h, f, d)
    t_direct = absolute_trace_polynomial_direct(F, h, f, F.e * d)

    # (1) the fast transitive implementation must equal the LITERAL formula
    #     sum_{i=0}^{ed-1} h^{2^i}, which is what the Blueprint specifies.
    assert t_fast == t_direct, "transitive trace != literal direct sum"

    # (2) the BUGGY version (i starting at 1, i.e. omitting h^{2^0} = h) must
    #     give a DIFFERENT answer. This is the regression test for DEFECT (a).
    buggy = P.trim([])
    cur = P.pmod(F, h, f)
    for _ in range(1, F.e * d):
        cur = P.pmod(F, P.pmul(F, cur, cur), f)
        buggy = P.padd(F, buggy, cur)
    assert buggy != t_fast, "the i=0 term does not distinguish this instance"

    # (3) Tr(h) must land in the prime field F_2, i.e. satisfy y^2 = y mod f.
    sq = P.pmod(F, P.pmul(F, t_fast, t_fast), f)
    assert sq == t_fast, "Tr(h) must lie in F_2, i.e. satisfy y^2 = y"

    # (4) FUNCTIONAL requirement: across random h, the CORRECT trace must
    #     actually split f non-trivially at least sometimes -- otherwise EDF
    #     could never terminate. (The buggy trace splits never, because it is
    #     not F_2-valued and gcd(f, buggy) stays trivial.)
    def buggy_trace(hh):
        """The DEFECT (a) version: sum_{i=1}^{ed-1} h^{2^i}, i.e. i starts at 1.
        Note the accumulator starts at ZERO and the first addition is h^2."""
        acc = P.trim([])
        cur = P.pmod(F, hh, f)
        for _ in range(1, F.e * d):
            cur = P.pmod(F, P.pmul(F, cur, cur), f)
            acc = P.padd(F, acc, cur)
        return acc

    rng2 = random.Random(999)
    splits_correct = 0
    for _ in range(24):
        hh = P.trim([rng2.randrange(F.q) for _ in range(P.deg(f))])
        if P.deg(hh) <= 0:
            continue
        tc = absolute_trace_polynomial(F, hh, f, d)
        gc = P.pgcd(F, f, tc)
        if 0 < P.deg(gc) < P.deg(f):
            splits_correct += 1
    assert splits_correct > 0, "correct trace never splits -- EDF cannot work"

    # (5) The DEFINING property the EDF relies on is that Tr(h) is F_2-valued
    #     (y^2 = y). The i=0-omitting version is NOT F_2-valued: it equals
    #     Tr(h) + h, whose square is Tr(h) + h^2 != Tr(h) + h in general. That
    #     is the real signature of DEFECT (a) -- it is not that the buggy
    #     trace can never split by accident (it can), but that it is not a
    #     trace at all, so splitting is unreliable.
    not_f2_valued = 0
    for _ in range(24):
        hh = P.trim([rng2.randrange(F.q) for _ in range(P.deg(f))])
        if P.deg(hh) <= 0:
            continue
        bg = buggy_trace(hh)
        if P.pmod(F, P.pmul(F, bg, bg), f) != bg:
            not_f2_valued += 1
    assert not_f2_valued > 0, (
        "the i=0-omitting 'trace' unexpectedly IS F_2-valued; "
        "this test instance does not exhibit DEFECT (a)"
    )


def test_defect_a_even_characteristic_factorisation_completes():
    """Full factorisation in even characteristic must terminate and be
    correct -- the end-to-end consequence of the defect-(a) repair."""
    for e in (1, 2, 3, 4, 5, 6, 8):
        F = FiniteField(2, e)
        rng = random.Random(4711)
        for _ in range(6):
            f = P.pmonic(F, [rng.randrange(F.q) for _ in range(rng.randrange(2, 9))] + [1])
            fac = factor_poly(F, f, seed=rng.randrange(10 ** 6))
            prod = [1]
            for irr, m in fac:
                assert is_irreducible(F, irr)
                for _ in range(m):
                    prod = P.pmul(F, prod, irr)
            assert P.trim(prod) == f


# --------------------------------------------------------------------------
# DEFECT (b) -- nu > e
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e,n", [
    (3, 1, 9),    # nu=2 > e=1   <-- the Blueprint's named crash case
    (3, 1, 27),   # nu=3 > e=1
    (5, 1, 25),   # nu=2 > e=1
    (3, 2, 27),   # nu=3 > e=2
    (2, 1, 4),    # nu=2 > e=1 (even characteristic as well)
    (2, 1, 8),    # nu=3 > e=1
])
def test_defect_b_nu_gt_e(p, e, n):
    """x^n - lambda must factor correctly when nu = v_p(n) exceeds e.
    The Blueprint's named crash case is q = 3, n = 9 (nu=2 > e=1)."""
    F = FiniteField(p, e)
    nu = v_p(n, p)
    assert nu > e, f"test setup: nu={nu} should exceed e={e}"
    for lam in F.nonzero_elements():
        n_prime, nu2, Pw, mu, fac = factor_xn_minus_lambda(F, n, lam)
        assert nu2 == nu
        assert Pw == p ** nu
        assert F.powi(mu, Pw) == lam, "mu^{p^nu} must equal lambda"
        prod = [1]
        for f, m in fac:
            assert is_irreducible(F, f)
            for _ in range(m):
                prod = P.pmul(F, prod, f)
        target = [F.neg(lam)] + [0] * (n - 1) + [1]
        assert P.trim(prod) == P.trim(target)


def test_defect_b_mu_is_the_unique_p_power_nu_root():
    F = FiniteField(3, 1)
    # q=3, e=1, n=9 -> nu=2, P=9. In F_3, x -> x^9 = x, so mu = lam.
    for lam in (1, 2):
        _np, nu, Pw, mu, _fac = factor_xn_minus_lambda(F, 9, lam)
        assert nu == 2 and Pw == 9
        assert F.powi(mu, 9) == lam
        assert mu == F.frobenius(lam, (1 - (2 % 1)) % 1)


# --------------------------------------------------------------------------
# x^n - lambda structural invariants (used by the P1 regression)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e", [(2, 1), (2, 4), (2, 8), (3, 1), (3, 4), (5, 2)])
@pytest.mark.parametrize("n", [3, 5, 7, 9, 15, 21, 25, 35])
def test_xn_minus_lambda_structure(p, e, n):
    F = FiniteField(p, e)
    rng = random.Random(n * 100 + p)
    for lam in [F.one, F.gen, F.gen_pow((F.q - 1) // 2)]:
        n_prime, nu, Pw, mu, fac = factor_xn_minus_lambda(F, n, lam)
        # n = n' p^nu and gcd(n', p) = 1
        assert n_prime * (p ** nu) == n
        assert n_prime % p != 0
        # mu^P = lambda
        assert F.powi(mu, Pw) == lam
        # every factor irreducible, monic, and with multiplicity exactly P
        for f, m in fac:
            assert m == Pw
            assert P.plead(F, f) == 1
            assert is_irreducible(F, f)
        # product identity against x^n - lambda
        prod = [1]
        for f, m in fac:
            for _ in range(m):
                prod = P.pmul(F, prod, f)
        assert P.trim(prod) == P.trim([F.neg(lam)] + [0] * (n - 1) + [1])
        # degree additivity
        assert sum(P.deg(f) * m for f, m in fac) == n


def test_factor_xn_minus_lambda_rejects_bad_input():
    F = FiniteField(3, 2)
    with pytest.raises(ValueError):
        factor_xn_minus_lambda(F, 0, F.one)
    with pytest.raises(ValueError):
        factor_xn_minus_lambda(F, 5, 0)


# --------------------------------------------------------------------------
# DDF / EDF contracts
# --------------------------------------------------------------------------
@pytest.mark.parametrize("p,e", [(2, 4), (3, 2), (5, 2)])
def test_ddf_groups_by_degree(p, e):
    F = FiniteField(p, e)
    rng = random.Random(1357)
    for _ in range(15):
        # product of distinct irreducibles of assorted degrees
        f = [1]
        degs = []
        for _ in range(rng.randrange(1, 4)):
            d = rng.randrange(1, 4)
            for cand in _all_monic(F, d):
                if is_irreducible(F, cand) and P.pmod(F, f, cand) != []:
                    f = P.pmul(F, f, cand)
                    degs.append(d)
                    break
        if not degs:
            continue
        parts = ddf(F, f)
        # every returned piece must be a product of equal-degree irreducibles
        for g, d in parts:
            assert P.deg(g) % d == 0
            for irr in edf(F, g, d, random.Random(1)):
                assert P.deg(irr) == d
                assert is_irreducible(F, irr)
