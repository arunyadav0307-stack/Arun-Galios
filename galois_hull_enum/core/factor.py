"""Polynomial factorisation over F_{p^e}: distinct-degree + equal-degree.

INDEPENDENT IMPLEMENTATION. Written from scratch for P1; nothing is imported,
copied or adapted from the throwaway blueprinting prototype, which lives
outside this repository.

Two defects recorded in BLUEPRINT.md A.7.1 are handled explicitly, and each
has a dedicated regression test (see tests/test_factor.py):

  DEFECT (a) -- even characteristic.  The equal-degree splitting for q even
  uses the absolute trace  T(h) = sum_{i=0}^{ed-1} h^{2^i}  from F_{2^{ed}}
  down to F_2.  The Blueprint's root cause was that the i=0 term h^{2^0} = h
  was OMITTED, so gcd(T(h), f) was taken with the wrong polynomial and no
  non-trivial split occurred.  Both the fast (trace-transitivity) and the
  literal (direct-sum) implementations below INCLUDE i = 0.

  DEFECT (b) -- nu > e.  In x^n - lambda with n = n' p^nu we need the unique
  mu with mu^{p^nu} = lambda.  Since sigma = (x -> x^p) has order e,
  sigma^nu = sigma^{nu mod e} and therefore
        mu = sigma^{e - (nu mod e)}(lambda),
  NOT sigma^{e-nu}(lambda), which is undefined/fractional when nu > e.
  The crash case named in the Blueprint is q = 3, n = 9 (nu = 2 > e = 1).
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Sequence, Tuple

from .poly import (
    apply_frobenius, deg, frobenius_basis, is_zero, padd, pderiv, pdivmod,
    pgcd, pmod, pmul, pmonic, pneg, ppowmod, psub, trim,
)
from .field import prime_factors

__all__ = [
    "is_irreducible", "squarefree_decomposition", "ddf", "edf",
    "factor_squarefree", "factor_poly", "factor_xn_minus_lambda",
    "absolute_trace_polynomial", "absolute_trace_polynomial_direct",
    "v_p",
]


def v_p(n: int, p: int) -> int:
    """p-adic valuation of n (n > 0)."""
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


# --------------------------------------------------------------------------
# irreducibility  (the INDEPENDENT test used by the P1 regression, A.7.1)
# --------------------------------------------------------------------------
def is_irreducible(F, f: Sequence[int]) -> bool:
    """Standard test: f monic of degree m > 0 is irreducible iff
    x^{q^m} = x (mod f) and gcd(x^{q^{m/ell}} - x, f) = 1 for all prime ell | m.
    """
    f = pmonic(F, f)
    m = deg(f)
    if m <= 0:
        return False
    if m == 1:
        return True
    x = [0, 1]
    if ppowmod(F, x, F.q ** m, f) != x:
        return False
    for ell in prime_factors(m):
        g = ppowmod(F, x, F.q ** (m // ell), f)
        d = pgcd(F, psub(F, g, x), f)
        if deg(d) > 0:
            return False
    return True


# --------------------------------------------------------------------------
# squarefree decomposition (Yun), with the derivative-zero case
# --------------------------------------------------------------------------
def _pth_root(F, f: Sequence[int]) -> List[int]:
    """The unique g with g(x)^p = f(x), for f whose formal derivative is 0.

    If f' = 0 then f(x) = h(x^p), and h(x^p) = (g(x))^p where g has
    coefficients g_i = sigma^{-1}(f_{p i}) = sigma^{e-1}(f_{p i}).
    """
    p = F.p
    f = trim(f)
    m = deg(f)
    if m < 0:
        return []
    if m % p != 0:
        raise ValueError("p-th root requested but degree is not divisible by p")
    gcoef = [0] * (m // p + 1)
    for i in range(m // p + 1):
        c = f[p * i] if p * i < len(f) else 0
        if c:
            gcoef[i] = F.frobenius(c, F.e - 1)
    g = trim(gcoef)
    # mandatory self-check: g^p must reconstruct f
    chk = [0] * (p * max(0, len(g) - 1) + 1)
    for i, c in enumerate(g):
        chk[p * i] = F.frobenius(c, 1)
    if trim(chk) != f:
        raise ValueError("p-th root extraction failed self-check")
    return g


def squarefree_decomposition(F, f: Sequence[int]) -> List[Tuple[List[int], int]]:
    """f = prod_i g_i^i with the g_i squarefree and pairwise coprime.

    MUSSE'S ALGORITHM, not Yun's.

    This matters: Yun's algorithm is only valid when the characteristic does
    not divide any multiplicity. Over F_2 with f = (x+1)^3 it terminates
    after one iteration and reports multiplicity 1 instead of 3, because the
    integer 3 that drives the multiplicity counter is reduced mod p.
    Musser's algorithm tracks multiplicities by repeated gcd against the
    gcd(f, f') accumulator and recurses on p-th roots, which is valid in
    small characteristic. (Discovered by the unit test
    test_factor_poly_is_correct over F_2; see docs/P1_REPORT.md.)
    """
    f = pmonic(F, f)
    acc: Dict[int, List[int]] = {}
    if deg(f) <= 0:
        return []
    _musser(F, f, 1, acc)
    out = [(pmonic(F, g), i) for i, g in sorted(acc.items())]
    return out


def _musser(F, f: Sequence[int], mult_scale: int, acc: Dict[int, List[int]]) -> None:
    """Musser's squarefree decomposition; accumulates multiplicity -> product."""
    p = F.p
    f = pmonic(F, f)
    if deg(f) <= 0:
        return
    fp = pderiv(F, f)

    if is_zero(fp):
        # f = g^p ; recurse on g with all multiplicities scaled by p
        _musser(F, _pth_root(F, f), mult_scale * p, acc)
        return

    c = pgcd(F, f, fp)
    w = pdivmod(F, f, c)[0]
    i = 1
    while deg(w) > 0:
        y = pgcd(F, w, c)
        z = pdivmod(F, w, y)[0]
        if deg(z) > 0:
            m = i * mult_scale
            acc[m] = pmul(F, acc[m], pmonic(F, z)) if m in acc else pmonic(F, z)
        i += 1
        w = y
        c = pdivmod(F, c, y)[0]

    if deg(c) > 0:
        # the remaining part is a p-th power: recurse on its p-th root
        _musser(F, _pth_root(F, c), mult_scale * p, acc)


# --------------------------------------------------------------------------
# distinct-degree factorisation
# --------------------------------------------------------------------------
def ddf(F, f: Sequence[int]) -> List[Tuple[List[int], int]]:
    """[(g_d, d), ...] where g_d is the product of all irreducible factors of
    the SQUAREFREE monic f having degree exactly d."""
    fstar = pmonic(F, f)
    out: List[Tuple[List[int], int]] = []
    if deg(fstar) <= 0:
        return out

    basis = frobenius_basis(F, fstar)
    h = [0, 1]          # h = x^{q^d} mod fstar
    d = 0
    x = [0, 1]
    while deg(fstar) > 0:
        d += 1
        if 2 * d > deg(fstar):
            # nothing of degree >= d can remain except fstar itself
            out.append((pmonic(F, fstar), deg(fstar)))
            fstar = []
            break
        h = apply_frobenius(F, h, basis, fstar)
        g = pgcd(F, fstar, psub(F, h, x))
        if deg(g) > 0:
            out.append((pmonic(F, g), d))
            fstar = pdivmod(F, fstar, g)[0]
            if deg(fstar) <= 0:
                break
            h = pmod(F, h, fstar)
            basis = frobenius_basis(F, fstar)
    return out


# --------------------------------------------------------------------------
# equal-degree factorisation
# --------------------------------------------------------------------------
def _edf_split(F, f: Sequence[int], d: int, rng: random.Random) -> List[List[int]]:
    """Split f, a product of r >= 2 monic irreducibles all of degree d."""
    n = deg(f)
    r = n // d
    if r <= 1:
        return [pmonic(F, f)]

    factors: List[List[int]] = [pmonic(F, f)]
    guard = 0
    while len(factors) < r:
        guard += 1
        if guard > 200 * (r + 1):
            raise RuntimeError("EDF failed to split (giving up)")

        # random polynomial of degree < n
        h = [rng.randrange(F.q) for _ in range(n)]
        h = trim(h)
        if deg(h) <= 0:
            continue

        if F.p == 2:
            t = absolute_trace_polynomial(F, h, f, d)
            if deg(t) <= 0:
                continue
            g = t
        else:
            e = (F.q ** d - 1) // 2
            g1 = ppowmod(F, h, e, f)
            g = psub(F, g1, [1])
            if deg(g) <= 0:
                continue

        new_factors: List[List[int]] = []
        changed = False
        for u in factors:
            if deg(u) <= d:
                new_factors.append(u)
                continue
            c = pgcd(F, u, g)
            if 0 < deg(c) < deg(u):
                v = pdivmod(F, u, c)[0]
                new_factors.append(pmonic(F, c))
                new_factors.append(pmonic(F, v))
                changed = True
            else:
                new_factors.append(u)
        if changed:
            factors = new_factors
    return factors


def edf(F, f: Sequence[int], d: int, rng: random.Random) -> List[List[int]]:
    """Factor a product of equal-degree irreducibles into irreducibles."""
    f = pmonic(F, f)
    if deg(f) == 0:
        return []
    if deg(f) == d:
        return [f]
    return _edf_split(F, f, d, rng)


def factor_squarefree(F, f: Sequence[int], rng: random.Random) -> List[List[int]]:
    """Factor a SQUAREFREE monic polynomial into monic irreducibles."""
    out: List[List[int]] = []
    for g, d in ddf(F, f):
        out.extend(edf(F, g, d, rng))
    return out


def factor_poly(F, f: Sequence[int], seed: int = 0) -> List[Tuple[List[int], int]]:
    """Full factorisation: [(monic irreducible, multiplicity), ...]."""
    rng = random.Random(seed)
    f = pmonic(F, f)
    out: List[Tuple[List[int], int]] = []
    if deg(f) <= 0:
        return out
    for g, mult in squarefree_decomposition(F, f):
        for irr in factor_squarefree(F, g, rng):
            out.append((irr, mult))
    out.sort(key=lambda t: (deg(t[0]), t[0]))
    return out


# --------------------------------------------------------------------------
# absolute trace (even characteristic) -- DEFECT (a)
# --------------------------------------------------------------------------
def absolute_trace_polynomial(F, h: Sequence[int], f: Sequence[int], d: int) -> List[int]:
    """Tr_{F_{q^d}/F_2}(h) in F_q[x]/(f), for q = 2^e.  USED IN PRODUCTION.

    Uses transitivity of the trace:
        Tr_{F_{2^{ed}}/F_2} = Tr_{F_{2^e}/F_2} o Tr_{F_{2^{ed}}/F_{2^e}},
    i.e.  S = sum_{i=0}^{d-1} h^{q^i}   (the i = 0 term is h itself), then
          T = sum_{i=0}^{e-1} S^{2^i}   (the i = 0 term is S itself).
    Both sums start at i = 0.  This is the repair for DEFECT (a).
    """
    if F.p != 2:
        raise ValueError("absolute trace to F_2 requires even characteristic")
    if deg(f) <= 0:
        return []

    # inner trace F_{q^d} -> F_q : sum_{i=0}^{d-1} h^{q^i}
    basis = frobenius_basis(F, f)
    s = pmod(F, h, f)
    acc = list(s)
    cur = list(s)
    for _ in range(1, d):
        cur = apply_frobenius(F, cur, basis, f)
        acc = padd(F, acc, cur)

    # outer trace F_q -> F_2 : sum_{i=0}^{e-1} acc^{2^i}
    total = list(acc)
    sq = list(acc)
    for _ in range(1, F.e):
        sq = pmod(F, pmul(F, sq, sq), f)
        total = padd(F, total, sq)
    return trim(total)


def absolute_trace_polynomial_direct(F, h: Sequence[int], f: Sequence[int], m: int) -> List[int]:
    """LITERAL formula  sum_{i=0}^{m-1} h^{2^i} mod f,  m = e*d.

    Provided for cross-checking the fast transitive implementation above.
    The i = 0 term IS included (this is the DEFECT (a) repair).
    """
    if deg(f) <= 0:
        return []
    total = pmod(F, h, f)
    cur = pmod(F, h, f)
    for _ in range(1, m):
        cur = pmod(F, pmul(F, cur, cur), f)
        total = padd(F, total, cur)
    return trim(total)


# --------------------------------------------------------------------------
# x^n - lambda   (Algorithm E1, steps 1-3) -- includes DEFECT (b) fix
# --------------------------------------------------------------------------
def factor_xn_minus_lambda(F, n: int, lam: int):
    """Factor x^n - lambda over F_q.

    n = n' p^nu with gcd(n', p) = 1;  P = p^nu.
    x^n - lambda = (x^{n'} - mu)^{P} where mu^{P} = lambda.

    Returns (n_prime, nu, P, mu, factors) where `factors` is the list of
    [(monic irreducible, multiplicity)] of x^n - lambda. Every irreducible
    factor of the squarefree part x^{n'} - mu occurs with multiplicity P.

    DEFECT (b): mu = sigma^{e - (nu mod e)}(lambda). Using e - nu would be a
    negative (fractional) exponent when nu > e, which is precisely the crash
    recorded in BLUEPRINT.md A.7.1 for q = 3, n = 9.
    """
    if n <= 0:
        raise ValueError("n must be positive")
    if lam == 0:
        raise ValueError("lambda must be nonzero (constacyclic)")

    p, e = F.p, F.e
    nu = v_p(n, p)
    P = p ** nu
    n_prime = n // P

    # mu = sigma^{e - (nu mod e)}(lam);  sigma has order e, so sigma^nu =
    # sigma^{nu mod e} and mu^{p^nu} = sigma^{nu mod e}(mu) = lam.
    shift = (e - (nu % e)) % e
    mu = F.frobenius(lam, shift)

    # mandatory self-check: mu^{P} must equal lambda
    if F.powi(mu, P) != lam:
        raise AssertionError(
            f"mu^(p^nu) != lambda  (q={F.q}, n={n}, nu={nu}, P={P}); "
            f"this is the DEFECT (b) failure mode"
        )

    if n_prime < 1:
        raise ValueError("n_prime must be >= 1")

    # x^{n'} - mu : constant term -mu, then zeros, then the leading 1.
    g = [F.neg(mu)] + [0] * (n_prime - 1) + [1]

    # x^{n'} - mu is squarefree: its derivative is n' x^{n'-1} and gcd(n',p)=1

    sqf = pgcd(F, g, pderiv(F, g))
    if deg(sqf) != 0:
        raise AssertionError(
            "x^n_prime - mu is not squarefree "
            f"(q={F.q}, n={n}, n_prime={n_prime}, nu={nu})"
        )

    if n_prime == 1:
        dist = [pmonic(F, g)]
    else:
        dist = factor_squarefree(F, g, random.Random(0))

    factors = [(f, P) for f in dist]
    return n_prime, nu, P, mu, factors
