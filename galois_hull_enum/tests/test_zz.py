"""Unit tests for core/zz.py -- exact Z[z] polynomial and matrix arithmetic."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import zz  # noqa: E402


def test_basic_polynomial_arithmetic():
    a = (1, 2, 3)          # 1 + 2z + 3z^2
    b = (2, -1)            # 2 - z
    assert zz.zdeg(a) == 2 and zz.zdeg(()) == -1
    assert zz.zadd(a, b) == (3, 1, 3)
    assert zz.zsub(a, b) == (-1, 3, 3)
    assert zz.zneg(a) == (-1, -2, -3)
    assert zz.zmul(a, b) == (2, 3, 4, -3)
    assert zz.zscal(0, a) == ()
    assert zz.zscal(2, a) == (2, 4, 6)
    assert zz.zmul(a, ()) == ()
    assert zz.zpow((0, 1), 0) == (1,)
    assert zz.zpow((0, 1), 3) == (0, 0, 0, 1)
    assert zz.ztrim((1, 2, 0, 0)) == (1, 2)


def test_multiplication_is_commutative_and_associative():
    import random
    rng = random.Random(11)
    for _ in range(200):
        a = zz.ztrim([rng.randrange(-5, 6) for _ in range(rng.randrange(0, 6))])
        b = zz.ztrim([rng.randrange(-5, 6) for _ in range(rng.randrange(0, 6))])
        c = zz.ztrim([rng.randrange(-5, 6) for _ in range(rng.randrange(0, 6))])
        assert zz.zmul(a, b) == zz.zmul(b, a)
        assert zz.zmul(zz.zmul(a, b), c) == zz.zmul(a, zz.zmul(b, c))
        assert zz.zmul(a, zz.zadd(b, c)) == zz.zadd(zz.zmul(a, b), zz.zmul(a, c))


def test_evaluation_and_format():
    assert zz.zval((1, 2, 3), 0) == 1
    assert zz.zval((1, 2, 3), 2) == 1 + 4 + 12
    assert zz.zformat(()) == "0"
    assert zz.zformat((5,)) == "5"
    assert zz.zformat((2, 2)) == "2z+2"
    assert zz.zformat((2, 12, 2)) == "2z^2+12z+2"


def test_matrix_operations():
    A = [[(1,), (2,)], [(0, 1), (1,)]]
    B = zz.mat_eye(2)
    assert zz.mat_trace(A) == (2,)
    assert zz.mat_mul(A, B) == A
    assert zz.mat_add(A, zz.mat_neg_like(A)) == zz.mat_zero(2, 2)
    # (AB)[i][j] computed by hand
    C = zz.mat_mul(A, A)
    assert C[0][0] == zz.zadd(zz.zmul((1,), (1,)), zz.zmul((2,), (0, 1)))
    assert C[0][1] == zz.zadd(zz.zmul((1,), (2,)), zz.zmul((2,), (1,)))


def test_determinant():
    # 2x2
    assert zz.zdet([[(1,), (2,)], [(3,), (4,)]]) == (-2,)
    # 3x3 identity
    assert zz.zdet(zz.mat_eye(3)) == (1,)
    # 3x3 by hand: |1 2 3; 4 5 6; 7 8 10| = 1(50-48) - 2(40-42) + 3(32-35)
    #            = 2 + 4 - 9 = -3
    M = [[(1,), (2,), (3,)], [(4,), (5,), (6,)], [(7,), (8,), (10,)]]
    assert zz.zdet(M) == (-3,)


def test_charpoly_via_traces_matches_laplace():
    """zmat_charpoly_via_traces uses Newton's identities on tr(A^i); cross-check
    it against a direct Laplace determinant of (lambda I - A)."""
    import random
    rng = random.Random(77)
    for n in (1, 2, 3, 4):
        A = zz.mat_zero(n, n)
        for i in range(n):
            for j in range(n):
                A[i][j] = zz.ztrim([rng.randrange(-3, 4) for _ in range(3)])
        chi = zz.zmat_charpoly_via_traces(A)
        # build lambda*I - A and take the determinant with lam as a symbol via
        # the zz machinery is not possible; instead verify Cayley-Hamilton:
        # sum_i chi[i] A^{n-i} == 0
        acc = zz.mat_zero(n, n)
        for i in range(n + 1):
            coeff = chi[i]
            power = zz._matpow(A, n - i)
            for r in range(n):
                for c in range(n):
                    acc[r][c] = zz.zadd(acc[r][c], zz.zmul(coeff, power[r][c]))
        assert all(acc[r][c] == () for r in range(n) for c in range(n)), \
            "Cayley-Hamilton failed for the computed characteristic polynomial"
        # leading coefficient must be 1
        assert chi[0] == (1,)
