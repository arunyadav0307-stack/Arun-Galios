r"""Cycle-local enumeration: local weights, cycle polynomials, transfer matrix.

This module implements, verbatim and independently, the frozen objects of
`BLUEPRINT.md` sections 3.2 and 3.3:

  Definition 2  w_a^{(P)}(u) = sum_{m=0}^{a-1} min{u_{m-1}, P - u_m},  u in [0,P]^a
  Definition 3  W_a^{(P)}(z) = sum_{u in [0,P]^a} z^{w_a^{(P)}(u)}
  Definition 5  \tilde W_a^{(P)} = tr(\tilde T_P^a),  \tilde T_P = T_P o (J - I)
  T3            W_a^{(P)}(z) = tr(T_P(z)^a),  T_P(z) = (z^{min{x,P-y}})_{0<=x,y<=P}

INDEPENDENT IMPLEMENTATION. Nothing is imported from or adapted from the
throwaway blueprinting prototype, which lives outside this repository.
Exact arithmetic throughout: Z[z] polynomials are tuples of Python ints
(see core/zz.py), so no floating point ever enters a computation.
"""

from __future__ import annotations

import itertools
from typing import Dict, List, Sequence, Tuple

from . import zz
from .zz import Matrix, Z

__all__ = [
    "local_weight", "cycle_poly_bruteforce", "cycle_poly_bruteforce_restricted",
    "transfer_matrix", "restricted_transfer_matrix",
    "cycle_poly_trace", "cycle_poly_trace_restricted",
    "all_ones_minus_identity", "ones", "dict_to_poly",
    "bivariate_gf_trace", "det_I_minus_tT", "t_degree",
    "verify_recurrence", "closed_form_P1", "degree_bound", "max_weight_word",
    "hull_dim_from_words",
]


# --------------------------------------------------------------------------
# Definition 2 -- the local weight
# --------------------------------------------------------------------------
def local_weight(u: Sequence[int], P: int) -> int:
    """w_a^{(P)}(u) = sum_{m=0}^{a-1} min{u_{m-1}, P - u_m}, indices mod a.

    This is a NEAREST-NEIGHBOUR CYCLIC INTERACTION: w_a = sum_m phi(u_{m-1},u_m)
    with phi(x,y) = min{x, P-y}. That structure is what makes the transfer
    matrix work.
    """
    a = len(u)
    if a == 0:
        return 0
    total = 0
    for m in range(a):
        total += min(u[(m - 1) % a], P - u[m])
    return total


def hull_dim_from_words(words_by_cycle: Sequence[Tuple[Sequence[int], int]], P: int) -> int:
    """T1's RHS: sum over cycles of ord_j(q) * w_a^{(P)}(u).

    `words_by_cycle` is a sequence of (u, d) pairs, one per cycle.
    """
    total = 0
    for u, d in words_by_cycle:
        total += d * local_weight(u, P)
    return total


# --------------------------------------------------------------------------
# Definition 3 -- the cycle polynomial, by direct enumeration of local states
# --------------------------------------------------------------------------
def cycle_poly_bruteforce(P: int, a: int) -> Dict[int, int]:
    """W_a^{(P)}(z) as {exponent: coefficient}, by enumerating ALL words
    u in [0,P]^a and applying Definition 2 directly.

    This is the STRUCTURALLY INDEPENDENT route: it never mentions a transfer
    matrix. It is the oracle for T3.
    """
    counts: Dict[int, int] = {}
    for u in itertools.product(range(P + 1), repeat=a):
        w = local_weight(u, P)
        counts[w] = counts.get(w, 0) + 1
    return dict(sorted(counts.items()))


def cycle_poly_bruteforce_restricted(P: int, a: int) -> Dict[int, int]:
    """The restricted analogue: same enumeration, but only over words with NO
    cyclically adjacent equal exponents (u_m != u_{m-1} for all m).

    This is P3's sample space C(n,q,lambda); Definition 5 encodes it by zeroing
    the diagonal of T_P.
    """
    counts: Dict[int, int] = {}
    for u in itertools.product(range(P + 1), repeat=a):
        if any(u[m] == u[(m - 1) % a] for m in range(a)):
            continue
        w = local_weight(u, P)
        counts[w] = counts.get(w, 0) + 1
    return dict(sorted(counts.items()))


def _dict_to_poly(counts: Dict[int, int]) -> Z:
    """{exponent: coeff} -> Z[z] tuple."""
    if not counts:
        return ()
    top = max(counts)
    return zz.ztrim([counts.get(i, 0) for i in range(top + 1)])


def dict_to_poly(counts: Dict[int, int]) -> Z:
    return _dict_to_poly(counts)


# --------------------------------------------------------------------------
# T3 -- the transfer matrix
# --------------------------------------------------------------------------
def transfer_matrix(P: int) -> Matrix:
    """T_P(z) = ( z^{min{x, P-y}} )_{0 <= x,y <= P}  over Z[z].

    Index set: the exponents 0..P of a single irreducible factor, i.e. the
    possible values of u_m in Definition 1. ROW index x is the "previous"
    state u_{m-1}, COLUMN index y is the "current" state u_m, and the entry
    z^{min{x, P-y}} is exactly phi(x, y) = min{x, P-y} from Definition 2.
    """
    n = P + 1
    M = zz.mat_zero(n, n)
    for x in range(n):
        for y in range(n):
            e = min(x, P - y)
            M[x][y] = (1,) if e == 0 else (0,) * e + (1,)
    return M


def ones(n: int) -> Matrix:
    """The all-ones matrix J_n."""
    return [[(1,) for _ in range(n)] for _ in range(n)]


def all_ones_minus_identity(n: int) -> Matrix:
    """J_n - I_n."""
    M = ones(n)
    for i in range(n):
        M[i][i] = zz.zsub(M[i][i], (1,))
    return M


def restricted_transfer_matrix(P: int) -> Matrix:
    """\tilde T_P = T_P o (J - I): T_P with its diagonal zeroed (Definition 5)."""
    T = transfer_matrix(P)
    JmI = all_ones_minus_identity(P + 1)
    return zz.mat_hadamard(T, JmI)


def cycle_poly_trace(P: int, a: int) -> Z:
    """W_a^{(P)}(z) = tr(T_P(z)^a) -- T3, computed by exact matrix powering."""
    if a < 1:
        raise ValueError("a must be >= 1")
    T = transfer_matrix(P)
    Pw = zz._matpow(T, a)
    return zz.mat_trace(Pw)


def cycle_poly_trace_restricted(P: int, a: int) -> Z:
    """\tilde W_a^{(P)}(z).

    Definition 5: for a >= 2 it is tr(\tilde T_P^a); for a = 1 it is DEFINED
    to be W_1^{(P)}(z) (not tr(\tilde T_P), which is 0 because the diagonal
    is zeroed). This special case is honoured exactly.
    """
    if a < 1:
        raise ValueError("a must be >= 1")
    if a == 1:
        return cycle_poly_trace(P, 1)
    Tt = restricted_transfer_matrix(P)
    return zz.mat_trace(zz._matpow(Tt, a))


# --------------------------------------------------------------------------
# T3 corollaries -- bivariate GF, determinant degree, recurrence
# --------------------------------------------------------------------------
def bivariate_gf_trace(P: int, N: int) -> Z:
    """sum_{a=1}^{N} W_a^{(P)}(z) * t^a, truncated at a = N, expressed as a
    polynomial in t with Z[z] coefficients.

    Here `t` is represented by a second variable tracked manually: the result
    is a Z[z]-polynomial per power of t, returned as a list indexed by the
    power of t. (T3 states the full sum is tr(tT(I-tT)^{-1}).)
    """
    out: List[Z] = [()]                       # out[k] = coeff of t^k
    for a in range(1, N + 1):
        Wa = cycle_poly_trace(P, a)
        if len(out) < a + 1:
            out.extend([()] * (a + 1 - len(out)))
        out[a] = zz.zadd(out[a], Wa)
    return out


def det_I_minus_tT(P: int) -> List[Z]:
    """det(I - t T_P(z)) as a polynomial in t, Z[z] coefficients.

    Returns a list indexed by the power of t (0..P+1). T3 says the t-degree
    is exactly P+1.
    """
    n = P + 1
    T = transfer_matrix(P)
    I = zz.mat_eye(n)
    # det(I - tT) = sum_{k=0}^{n} (-1)^k t^k * (sum of k x k principal minors)
    # computed by expansion over subsets -- exact, and small enough here.
    coeffs: List[Z] = [()] * (n + 1)
    for k in range(n + 1):
        acc: Z = ()
        for rows in itertools.combinations(range(n), k):
            sub = [[T[i][j] for j in rows] for i in rows]
            acc = zz.zadd(acc, zz.zdet(sub))
        coeffs[k] = acc if k % 2 == 0 else zz.zneg(acc)
    # trim leading zero coefficients
    while coeffs and not coeffs[-1]:
        coeffs.pop()
    return coeffs


def t_degree(P: int) -> int:
    """t-degree of det(I - t T_P(z)) (T3 says P+1)."""
    c = det_I_minus_tT(P)
    return len(c) - 1


def verify_recurrence(P: int, a_max: int) -> Tuple[bool, List[Z], int]:
    """Verify that W_a^{(P)} satisfies a linear recurrence in a of order P+1
    with coefficients in Z[z] (T3, third bullet).

    Method: the sequence W_1, ..., W_{2(P+1)} determines the unique order-(P+1)
    recurrence; check it then predicts W_{P+2}, ..., W_{a_max}.
    Returns (ok, recurrence_coeffs, order).
    """
    order = P + 1
    seq = [cycle_poly_trace(P, a) for a in range(1, 2 * order + 2)]
    if len(seq) < 2 * order:
        return False, [], order
    # solve for c_0..c_{order-1} with  W_{a} = sum_{i=0}^{order-1} c_i W_{a-1-i}
    # using the `order` equations starting at a = order+1 .. 2*order.
    # Solve the linear system over Q(z) by Gaussian elimination on Z[z] with
    # rational reconstruction avoided: use exact fraction-free elimination via
    # sympy-free approach -> solve the integer-linear system by trying the
    # unique solution from Cayley-Hamilton instead.
    #
    # Cleaner and fully rigorous: the characteristic polynomial of T_P(z) gives
    # the recurrence by Cayley-Hamilton. chi_{T_P}(lambda) = lambda^{n} + c_1
    # lambda^{n-1} + ... + c_n, and T_P^n + c_1 T_P^{n-1} + ... + c_n I = 0.
    # Taking traces: W_{a+n} = -sum_{i=1}^{n} c_i W_{a+n-i}.
    chi = zz.zmat_charpoly_via_traces(transfer_matrix(P))
    # chi = [lambda^n coeff, ..., constant]; c_i = chi[i]
    # recurrence: W_{a+n} = -sum_{i=1..n} c_i W_{a+n-i}
    n = P + 1
    ok = True
    for a in range(1, a_max - n + 1):
        acc: Z = ()
        for i in range(1, n + 1):
            term = zz.zmul(chi[i], cycle_poly_trace(P, a + n - i))
            acc = zz.zadd(acc, term)
        acc = zz.zadd(acc, cycle_poly_trace(P, a + n))
        if acc != ():
            ok = False
            break
    # recurrence coefficients, as [c_1, ..., c_n] with W_{a+n} = -sum c_i W_{a+n-i}
    coeffs = [zz.zneg(chi[i]) for i in range(1, n + 1)]
    return ok, coeffs, n


# --------------------------------------------------------------------------
# T4 -- the P = 1 closed form
# --------------------------------------------------------------------------
def closed_form_P1(a: int) -> Z:
    """W_a^{(1)}(z) = (1+sqrt z)^a + (1-sqrt z)^a = 2 sum_{i<=a/2} C(a,2i) z^i.

    Computed from the BINOMIAL expansion, i.e. from the closed form, NOT from
    the transfer matrix -- so agreement with cycle_poly_trace(P=1, a) is a
    genuine check of T4.
    """
    import math
    out = [0] * (a // 2 + 1)
    for i in range(a // 2 + 1):
        out[i] = 2 * math.comb(a, 2 * i)
    return zz.ztrim(out)


# --------------------------------------------------------------------------
# T3/T4 degree statements
# --------------------------------------------------------------------------
def degree_bound(P: int, a: int) -> int:
    """floor(a P / 2) -- the claimed degree of W_a^{(P)}."""
    return (a * P) // 2


def max_weight_word(P: int, a: int) -> Tuple[int, ...]:
    """A word attaining the maximum local weight floor(aP/2), found by brute
    force. Needed for the degree claim: the degree of W equals the maximum
    weight attained, so a maximiser must be exhibited.
    """
    best = None
    best_w = -1
    for u in itertools.product(range(P + 1), repeat=a):
        w = local_weight(u, P)
        if w > best_w:
            best_w = w
            best = u
    return best, best_w
