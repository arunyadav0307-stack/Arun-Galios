"""Unit tests for core/poly.py -- dense polynomial arithmetic over F_q."""

import os
import random
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.field import FiniteField  # noqa: E402
from core import poly as P  # noqa: E402

FIELDS = [(2, 1), (2, 3), (2, 8), (3, 1), (3, 2), (3, 4), (5, 2), (7, 2)]


@pytest.mark.parametrize("p,e", FIELDS)
def test_arithmetic_identities(p, e):
    F = FiniteField(p, e)
    rng = random.Random(9001 + p)
    for _ in range(150):
        da, db = rng.randrange(0, 7), rng.randrange(0, 7)
        a = P.trim([rng.randrange(F.q) for _ in range(da)])
        b = P.trim([rng.randrange(F.q) for _ in range(db)])
        c = P.trim([rng.randrange(F.q) for _ in range(5)])
        # additive group
        assert P.padd(F, a, b) == P.padd(F, b, a)
        assert P.padd(F, P.padd(F, a, b), c) == P.padd(F, a, P.padd(F, b, c))
        assert P.psub(F, a, b) == P.padd(F, a, P.pneg(F, b))
        assert P.deg(P.psub(F, a, a)) == -1
        # multiplicative
        assert P.pmul(F, a, b) == P.pmul(F, b, a)
        assert P.pmul(F, P.pmul(F, a, b), c) == P.pmul(F, a, P.pmul(F, b, c))
        # distributivity
        assert P.pmul(F, a, P.padd(F, b, c)) == P.padd(F, P.pmul(F, a, b), P.pmul(F, a, c))
        # degree of a product
        if P.deg(a) >= 0 and P.deg(b) >= 0:
            assert P.deg(P.pmul(F, a, b)) == P.deg(a) + P.deg(b)
        # zero handling
        assert P.pmul(F, a, []) == []
        assert P.padd(F, a, []) == a


@pytest.mark.parametrize("p,e", FIELDS)
def test_divmod_identity(p, e):
    """a = q*b + r with deg r < deg b."""
    F = FiniteField(p, e)
    rng = random.Random(31337)
    for _ in range(150):
        da = rng.randrange(1, 9)
        db = rng.randrange(1, 5)
        a = P.trim([rng.randrange(F.q) for _ in range(da)])
        b = P.trim([rng.randrange(F.q) for _ in range(db)])
        if P.deg(b) < 0:
            continue
        q, r = P.pdivmod(F, a, b)
        assert P.deg(r) < P.deg(b) or P.deg(r) == -1
        assert P.padd(F, P.pmul(F, q, b), r) == a


@pytest.mark.parametrize("p,e", FIELDS)
def test_gcd_and_egcd(p, e):
    F = FiniteField(p, e)
    rng = random.Random(5150)
    for _ in range(80):
        a = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 7))])
        b = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 7))])
        if P.deg(a) < 0 or P.deg(b) < 0:
            continue
        g = P.pgcd(F, a, b)
        assert P.deg(g) >= 0
        assert P.plead(F, g) == 1          # monic
        assert P.pmod(F, a, g) == []
        assert P.pmod(F, b, g) == []
        gg, s, t = P.pegcd(F, a, b)
        assert gg == g
        assert P.padd(F, P.pmul(F, s, a), P.pmul(F, t, b)) == g


@pytest.mark.parametrize("p,e", FIELDS)
def test_powmod(p, e):
    F = FiniteField(p, e)
    rng = random.Random(1010)
    for _ in range(60):
        m = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(2, 6))])
        if P.plead(F, m) == 0:
            continue
        m = P.pmonic(F, m)
        base = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 5))])
        # a^(b+c) = a^b a^c  mod m
        b, c = rng.randrange(0, 40), rng.randrange(0, 40)
        lhs = P.ppowmod(F, base, b + c, m)
        rhs = P.pmod(F, P.pmul(F, P.ppowmod(F, base, b, m), P.ppowmod(F, base, c, m)), m)
        assert lhs == rhs
        assert P.ppowmod(F, base, 0, m) == [1]


@pytest.mark.parametrize("p,e", FIELDS)
def test_derivative_rules(p, e):
    F = FiniteField(p, e)
    rng = random.Random(2020)
    for _ in range(60):
        a = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 6))])
        b = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 6))])
        # (a+b)' = a' + b'
        assert P.pderiv(F, P.padd(F, a, b)) == P.padd(F, P.pderiv(F, a), P.pderiv(F, b))
        # (ab)' = a'b + ab'
        lhs = P.pderiv(F, P.pmul(F, a, b))
        rhs = P.padd(F, P.pmul(F, P.pderiv(F, a), b), P.pmul(F, a, P.pderiv(F, b)))
        assert lhs == rhs
    # derivative of x^n is n x^{n-1} reduced mod p
    for n in range(1, 6):
        xn = [0] * n + [1]
        d = P.pderiv(F, xn)
        expect_n = n % p
        if expect_n == 0 or n - 1 < 0:
            assert d == []
        else:
            assert d == P.pscal(F, F.from_int(expect_n), [0] * (n - 1) + [1])


@pytest.mark.parametrize("p,e", FIELDS)
def test_eval_is_a_ring_hom(p, e):
    F = FiniteField(p, e)
    rng = random.Random(3030)
    for _ in range(80):
        a = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 6))])
        b = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 6))])
        x = rng.randrange(F.q)
        assert P.peval(F, P.padd(F, a, b), x) == F.add(P.peval(F, a, x), P.peval(F, b, x))
        assert P.peval(F, P.pmul(F, a, b), x) == F.mul(P.peval(F, a, x), P.peval(F, b, x))


def test_frobenius_on_coefficients_and_the_quotient_ring():
    """frobenius_basis / apply_frobenius compute the q-power map on
    F_q[x]/(f). Cross-check against literal ppowmod."""
    for (p, e) in [(2, 3), (3, 2), (5, 2), (2, 8)]:
        F = FiniteField(p, e)
        rng = random.Random(6060)
        # use an irreducible f so that F_q[x]/(f) is a field
        from core.factor import squarefree_decomposition
        for trial in range(3):
            f = P.pmonic(F, [rng.randrange(F.q) for _ in range(4)] + [1])
            if P.deg(f) != 4:
                continue
            # only need f monic nonzero for the identity h^q = (h mod f)^q
            basis = P.frobenius_basis(F, f)
            assert len(basis) == P.deg(f)
            assert basis[0] == [1]
            for _ in range(6):
                h = P.trim([rng.randrange(F.q) for _ in range(rng.randrange(1, 5))])
                fast = P.apply_frobenius(F, h, basis, f)
                slow = P.ppowmod(F, h, F.q, f)
                assert fast == slow, "Frobenius basis disagrees with literal powering"


def test_frobenius_coeffs_is_sigma_on_each_coefficient():
    F = FiniteField(3, 4)
    a = [1, 2, 0, F.gen]
    assert P.frobenius_coeffs(F, a, 0) == P.trim(a)
    assert P.frobenius_coeffs(F, a, 1) == P.trim([F.frobenius(c, 1) for c in a])
    assert P.frobenius_coeffs(F, a, 4) == P.trim(a)     # sigma^e = id


def test_poly_from_int_coeffs_and_format():
    F = FiniteField(3, 3)
    # P1 Example 6: x^4 + 2x^3 + x^2 + 2x + 1
    f = P.poly_from_int_coeffs(F, [1, 2, 1, 2, 1])
    assert f == [1, 2, 1, 2, 1]
    assert P.format_poly(F, f) == "x^4+2x^3+x^2+2x+1"
    assert P.format_poly(F, []) == "0"
    assert P.format_poly(F, [1]) == "1"
    assert P.format_poly(F, [0, 1]) == "x"
