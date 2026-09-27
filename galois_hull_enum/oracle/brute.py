r"""Algorithm E3 -- the independent brute-force validation oracle.

BLUEPRINT.md section 5.3, verbatim:

    Iterate u in prod_c [0,P]^{a(c)}; build g = prod f_m^{u_m},
    h = prod f_m^{P-u_m}; h^{#} from the polynomial h;
    dim = n - deg lcm(g, h^{#}); histogram.
    Feasible while |C| <= 10^5 (pure Python) or <= 10^7 (NumPy-accelerated).

WHAT MAKES THIS AN ORACLE RATHER THAN A SECOND ENUMERATOR
    E3 computes the hull dimension of ONE code at a time from actual
    polynomials: it builds g, divides x^n - lambda by it, applies the #-map to
    the quotient, and takes a gcd. It never forms W_a^{(P)}, never multiplies
    per-cycle factors, and never sees the transfer matrix T_P. The cycle
    decomposition is used ONLY as an indexing device -- which of the (P+1)^t
    exponent vectors to visit -- never in the dimension formula.

TWO ORACLES ARE PROVIDED, ON PURPOSE
    oracle_distribution_cycle  -- E3 exactly as the Blueprint writes it:
                                  indexed by the cycle product prod_c [0,P]^a.
    oracle_distribution_flat   -- the same mathematics indexed by a FLAT
                                  exponent vector over the distinct
                                  irreducible factors, with NO use of the
                                  cycle decomposition at all.

    If E2, oracle_cycle and oracle_flat all agree, then the agreement cannot
    be an artefact of the cycle indexing, because oracle_flat shares no cycle
    code with E2. That is the answer to the "hidden implementation
    dependency" requirement; see docs/P3_REPORT.md section B.4.

AN INDEPENDENT #-MAP
    hash_map_independent is written directly from the frozen definition
    (BLUEPRINT section 3.1) and does NOT import core.cycles.hash_map, so it is
    a separate code path. certify_factorisation independently re-derives that
    the factorisation of x^n - lambda really is one.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from core import poly as P
from core.field import FiniteField
from core.factor import factor_xn_minus_lambda, is_irreducible

__all__ = [
    "hash_map_independent", "hull_dimension_poly",
    "oracle_distribution_cycle", "oracle_distribution_flat",
    "certify_factorisation", "oracle_stats",
]


# --------------------------------------------------------------------------
# an INDEPENDENT #-map (frozen definition, fresh implementation)
# --------------------------------------------------------------------------
def hash_map_independent(F: FiniteField, f: Sequence[int], j: int) -> List[int]:
    r"""f^{#}(x) = sum_{i=0}^{m} f_0^{-p^j} f_i^{p^j} x^{m-i}.

    Written straight from BLUEPRINT section 3.1, with sigma^j applied one
    coefficient at a time. core.cycles.hash_map is NOT imported.
    """
    f = P.trim(list(f))
    m = P.deg(f)
    if m < 0:
        raise ValueError("f^# is not defined for the zero polynomial")
    f0 = f[0]
    if f0 == 0:
        raise ValueError("f^# requires f(0) != 0")

    def sigma_pow(a: int) -> int:
        out = a
        for _ in range(j):
            out = F.frobenius(out, 1)
        return out

    scale = F.inv(sigma_pow(f0))
    out = [F.mul(scale, sigma_pow(f[i])) for i in range(m + 1)]
    out.reverse()                       # x^{m-i}
    return P.trim(out)


def hull_dimension_poly(F: FiniteField, n: int, lam: int, k: int,
                        g: Sequence[int]) -> int:
    """dim hull_k(C) = n - deg lcm(g, h^{#}) with h = (x^n - lambda)/g.

    This is P1 Remark 2 / BLUEPRINT section 2.0.4. It is the ONLY place the
    oracle's dimension comes from; the transfer matrix is nowhere involved.
    """
    j = F.e - k
    xn_minus_lam = [F.neg(lam)] + [0] * (n - 1) + [1]
    _q, r = P.pdivmod(F, xn_minus_lam, g)
    if r:
        raise ValueError("g does not divide x^n - lambda")
    h = _q
    h_sharp = hash_map_independent(F, h, j)
    g_monic = P.pmonic(F, g)
    gcd = P.pgcd(F, g_monic, h_sharp)
    if not gcd:
        raise AssertionError("gcd(g, h^#) is zero")
    lcm = P.pmonic(F, P.pmul(F, g_monic, P.pdivmod(F, h_sharp, gcd)[0]))
    if P.deg(lcm) > n:
        raise AssertionError(
            f"deg lcm(g,h^#) = {P.deg(lcm)} > n = {n}; the frozen identity "
            f"requires lcm | x^n - lambda, so hypothesis (H) is violated")
    return n - P.deg(lcm)


# --------------------------------------------------------------------------
# independent certification of the factorisation both sides consume
# --------------------------------------------------------------------------
def certify_factorisation(F: FiniteField, n: int, lam: int,
                          factors: Sequence[Tuple[Sequence[int], int]]
                          ) -> Dict[str, object]:
    """Re-derive, from scratch, that `factors` really is the factorisation of
    x^n - lambda.

    Checks, independently of core/factor.py's internals:
      1. prod f_i^{mult_i} == x^n - lambda            (as polynomials)
      2. each f_i is monic, irreducible, and f_i(0) != 0
      3. the factors are pairwise distinct
      4. every multiplicity is the same P
    Raises AssertionError on any failure.
    """
    xn = [F.neg(lam)] + [0] * (n - 1) + [1]
    prod = [1]
    mults = set()
    for (f, m) in factors:
        fm = P.pmonic(F, f)
        if P.deg(fm) < 1:
            raise AssertionError("factor is constant")
        if not is_irreducible(F, fm):
            raise AssertionError("factor is not irreducible")
        if fm[0] == 0:
            raise AssertionError("factor has zero constant term")
        for _ in range(m):
            prod = P.pmul(F, prod, fm)
        mults.add(m)
    if P.trim(prod) != P.trim(xn):
        raise AssertionError("product of factors != x^n - lambda")
    keys = [tuple(P.pmonic(F, f)) for (f, _m) in factors]
    if len(set(keys)) != len(keys):
        raise AssertionError("factors are not pairwise distinct")
    if len(mults) != 1:
        raise AssertionError(f"multiplicities not uniform: {mults}")
    return {"product_ok": True, "irreducible_ok": True, "distinct_ok": True,
            "uniform_multiplicity": True, "P": mults.pop(),
            "n_factors": len(factors)}


# --------------------------------------------------------------------------
# odometer over words of mixed lengths (shared by both oracles)
# --------------------------------------------------------------------------
def _incr_word(w: List[int], Pv: int) -> Optional[List[int]]:
    """Increment a mixed-radix digit string (each digit in [0, Pv]).

    Returns the new word, or None on overflow.
    """
    out = list(w)
    for i in range(len(out) - 1, -1, -1):
        if out[i] < Pv:
            out[i] += 1
            return out
        out[i] = 0
    return None


# --------------------------------------------------------------------------
# Oracle A -- E3 exactly as the Blueprint writes it (cycle-indexed)
# --------------------------------------------------------------------------
def oracle_distribution_cycle(F: FiniteField, cd, n: int, lam: int, k: int
                              ) -> Dict[int, int]:
    """E3 indexed by the cycle product: u in prod_c [0,P]^{a(c)}.

    `cd` is a core.cycles.CycleData. Every code is materialised as an actual
    generator polynomial and its hull dimension computed by the polynomial
    route.
    """
    Pv = cd.P
    lengths = [c["a"] for c in cd.cycles]
    hist: Dict[int, int] = {}
    if not lengths:
        # no factors: only the zero code, g = 1
        hist[hull_dimension_poly(F, n, lam, k, [1])] = 1
        return hist

    words = [[0] * a for a in lengths]
    while True:
        g = [1]
        for cyc, u in zip(cd.cycles, words):
            for idx, exp in zip(cyc["orbit"], u):
                if exp:
                    f = cd.factors[idx][0]
                    for _ in range(exp):
                        g = P.pmul(F, g, f)
        d = hull_dimension_poly(F, n, lam, k, g)
        hist[d] = hist.get(d, 0) + 1

        i = len(words) - 1
        nxt = _incr_word(words[i], Pv)
        while nxt is None and i > 0:
            words[i] = [0] * lengths[i]
            i -= 1
            nxt = _incr_word(words[i], Pv)
        if nxt is None:
            return hist
        words[i] = nxt


# --------------------------------------------------------------------------
# Oracle B -- flat enumeration, NO cycle decomposition anywhere
# --------------------------------------------------------------------------
def oracle_distribution_flat(F: FiniteField, n: int, lam: int, k: int,
                             factors: Sequence[Tuple[Sequence[int], int]],
                             with_certification: bool = True
                             ) -> Tuple[Dict[int, int], Dict[str, object]]:
    """E3 indexed by a FLAT exponent vector over the distinct irreducible
    factors of x^n - lambda.

    This is the hidden-dependency-free oracle. It calls
    core.factor.factor_xn_minus_lambda ONCE to learn the factors (the CP1
    kernel, independently certified here and in P1) and then enumerates every
    monic divisor of x^n - lambda. It never calls core.cycles, never forms a
    cycle, and never touches T_P.

    Returns (histogram, stats).
    """
    if with_certification:
        certify_factorisation(F, n, lam, factors)

    Pv = factors[0][1]
    flist = [P.pmonic(F, f) for (f, _m) in factors]
    t = len(flist)
    if t == 0:
        return {0: 1}, {"n_factors": 0, "n_codes": 1, "P": Pv}

    # precompute f_i^e for e = 0..P
    powers: List[List[List[int]]] = []
    for f in flist:
        row = [[1]]
        for _ in range(Pv):
            row.append(P.pmul(F, row[-1], f))
        powers.append(row)

    # g is rebuilt from the precomputed powers on EVERY step. An earlier
    # version cached suffix[i] = prod_{l>=i} f_l^{e_l} and updated it
    # incrementally; that cache went STALE on an odometer carry (the wrapped
    # digit was reset to 0 without refreshing its suffix entry), which silently
    # produced wrong g's and a wrong histogram on exactly the instances with
    # a >= 3 cycles. The cache was removed rather than patched: correctness
    # matters more than a constant factor here, and the cost is dominated by
    # hull_dimension_poly's gcd anyway. See docs/P3_REPORT.md section G.
    exps = [0] * t
    hist: Dict[int, int] = {}
    n_codes = 0
    while True:
        g = [1]
        for i in range(t):
            g = P.pmul(F, g, powers[i][exps[i]])
        d = hull_dimension_poly(F, n, lam, k, g)
        hist[d] = hist.get(d, 0) + 1
        n_codes += 1

        i = t - 1
        while i >= 0:
            exps[i] += 1
            if exps[i] <= Pv:
                break
            exps[i] = 0
            i -= 1
        if i < 0:
            break

    stats = {
        "n_factors": t,
        "n_codes": n_codes,
        "n_codes_expected": (Pv + 1) ** t,
        "P": Pv,
        "n_codes_match_expected": n_codes == (Pv + 1) ** t,
    }
    return hist, stats


def oracle_stats(hist: Dict[int, int]) -> Dict[str, object]:
    """Summary statistics of a histogram (exact rationals as strings)."""
    total = sum(hist.values())
    if total == 0:
        return {"n_codes": 0}
    s = sum(d * c for d, c in hist.items())
    return {
        "n_codes": total,
        "sum_dims": s,
        "mean_exact": f"{s}/{total}",
        "min_dim": min(hist),
        "max_dim": max(hist),
        "support": sorted(hist),
        "lcd_count": hist.get(0, 0),
    }
