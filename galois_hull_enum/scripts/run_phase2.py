#!/usr/bin/env python3
r"""P2 -- cycle-polynomial / transfer-matrix foundation: full verification.

Executes ONLY P2-sized computations:
  * local weights and cycle polynomials (Definitions 2, 3)
  * the transfer matrix and T3's trace representation
  * T3's corollaries (bivariate GF, det degree, recurrence, degree, LCD)
  * T4's closed form at P = 1
  * Definition 5's restricted version
  * T1, the frozen hull-dimension identity, on real instances from P1
  * the PA-1 specialisation boundary
  * edge cases a = 1, a = 2, P = 1, a >= 3

No enumerator N over all codes (that is P3 / Algorithm E2); no moments
(P6); no limit law (P8).

Usage:  python3 scripts/run_phase2.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core import zz                                    # noqa: E402
from core import poly as P                             # noqa: E402
from core.cycle_poly import (                          # noqa: E402
    closed_form_P1, cycle_poly_bruteforce,
    cycle_poly_bruteforce_restricted, cycle_poly_trace,
    cycle_poly_trace_restricted, degree_bound, det_I_minus_tT,
    dict_to_poly, hull_dim_from_words, local_weight,
    restricted_transfer_matrix, t_degree, transfer_matrix,
    verify_recurrence,
)
from core.field import FiniteField                      # noqa: E402
from core.hull_dim import verify_T1_on_instance        # noqa: E402

RESULTS = ROOT / "results"
OUT_JSON = RESULTS / "phase2_summary.json"
OUT_CSV = RESULTS / "phase2_cycle_polys.csv"

summary = {
    "phase": "P2",
    "blueprint_commit": "076f9109142f7965242c365aa90f664e4acfb56b",
    "plan_commit": "96c78859cf2ec77763ef6bf1adc8e5f11260dcec",
    "p0_commit": "076f97770c802d9af205bea2224ed44eb4b0c043",
    "p1_commit": "7d467652718eb78b33120c37650ccdbd96e22e90",
}
checks = []          # (name, ok, detail)


def record(name, ok, detail=""):
    checks.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    return bool(ok)


def hr(title):
    print("\n" + "=" * 76)
    print(title)
    print("=" * 76)


# ==========================================================================
# 1. T3 -- the trace representation, brute force vs transfer matrix
# ==========================================================================
def check_trace_vs_bruteforce():
    hr("T3  W_a^{(P)}(z) = tr(T_P(z)^a)  --  brute force over local states vs trace")
    rows = []
    cases = []
    for a in (1, 2, 3, 4, 5, 6, 8, 12):
        cases.append((1, a))
    for P in (2, 3, 4, 5):
        for a in (1, 2, 3, 4):
            cases.append((P, a))
    n_ok = 0
    for P, a in cases:
        bf = dict_to_poly(cycle_poly_bruteforce(P, a))
        tr = cycle_poly_trace(P, a)
        ok = (bf == tr)
        n_ok += ok
        rows.append({"P": P, "a": a, "bruteforce": zz.zformat(bf),
                     "trace": zz.zformat(tr), "agree": ok})
        if not ok:
            print(f"    MISMATCH P={P} a={a}: bf={zz.zformat(bf)} trace={zz.zformat(tr)}")
    record(f"T3 trace == brute-force local-state enumeration ({len(cases)} (P,a) pairs)",
           n_ok == len(cases), f"{n_ok}/{len(cases)}")
    return rows


# ==========================================================================
# 2. T3 corollaries
# ==========================================================================
def check_t3_corollaries():
    hr("T3 corollaries")
    # (a) t-degree of det(I - t T_P) is P+1
    degs = {P: t_degree(P) for P in range(0, 6)}
    record("deg_t det(I - t T_P(z)) = P+1  for P = 0..5",
           all(degs[P] == P + 1 for P in degs),
           str(degs))
    # (b) at z = 1, T_P(1) = J_{P+1} so det(I - t T_P(1)) = 1 - (P+1) t.
    # NB: evaluate then TRIM trailing zeros -- the top t-coefficient of
    # det(I - t T_P(z)) is a power of (1-z), which vanishes at z = 1.
    for P in range(0, 6):
        coeffs = det_I_minus_tT(P)
        at1 = [zz.zval(c, 1) for c in coeffs]
        while at1 and at1[-1] == 0:
            at1.pop()
        if at1 != [1, -(P + 1)]:
            record(f"det(I - t T_{P}(1)) = 1 - {P+1} t", False, f"got {at1}")
            return
    record("det(I - t T_P(1)) = 1 - (P+1) t   (T_P(1) = J_{P+1})", True,
           "P = 0..5")
    # (c) order-(P+1) recurrence via Cayley-Hamilton
    recs = {}
    allok = True
    for P in range(0, 6):
        ok, coeffs, order = verify_recurrence(P, 2 * (P + 1) + 6)
        recs[P] = {"order": order, "verified": ok}
        allok = allok and ok
    record("W_a^{(P)} satisfies an order-(P+1) recurrence in a over Z[z] "
           "(Cayley-Hamilton on T_P)", allok,
           "P = 0..5")
    # (d) deg W = floor(aP/2)
    degcases = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
                (2, 1), (2, 2), (2, 3), (2, 4), (2, 5),
                (3, 1), (3, 2), (3, 3), (3, 4),
                (4, 1), (4, 2), (4, 3), (5, 1), (5, 2), (5, 3)]
    nok = 0
    for P, a in degcases:
        d = zz.zdeg(dict_to_poly(cycle_poly_bruteforce(P, a)))
        nok += (d == degree_bound(P, a))
    record("deg W_a^{(P)} = floor(aP/2)", nok == len(degcases),
           f"{nok}/{len(degcases)}")
    # (e) W(1) = (P+1)^a
    nok = 0
    for P, a in [(1, 1), (1, 5), (1, 12), (2, 3), (3, 4), (5, 2)]:
        nok += (zz.zval(cycle_poly_trace(P, a), 1) == (P + 1) ** a)
    record("W_a^{(P)}(1) = (P+1)^a  (total word count)", nok == 6, f"{nok}/6")
    return degs


# ==========================================================================
# 3. T5 -- LCD count c_{a,0} = 2
# ==========================================================================
def check_t5():
    hr("T5  c_{a,0}^{(P)} = 2, extremal words 0^a and P^a")
    cases = [(1, 1), (1, 2), (1, 3), (1, 5), (2, 1), (2, 2), (2, 3),
             (3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (5, 1)]
    nok = 0
    for P, a in cases:
        poly = dict_to_poly(cycle_poly_bruteforce(P, a))
        nok += (zz.zcoeff(poly, 0) == 2)
        # and confirm the two extremal words are exactly the weight-0 words
        zero_words = [u for u in _words(P, a) if local_weight(u, P) == 0]
        if zero_words != [tuple([0] * a), tuple([P] * a)]:
            nok -= 1
    record("c_{a,0}^{(P)} = 2 with extremal words 0^a and P^a", nok == len(cases),
           f"{nok}/{len(cases)}")


def _words(P, a):
    import itertools
    return list(itertools.product(range(P + 1), repeat=a))


# ==========================================================================
# 4. T4 -- the P = 1 closed form
# ==========================================================================
def check_t4():
    hr("T4  W_a^{(1)}(z) = (1+sqrt z)^a + (1-sqrt z)^a = 2 sum_i C(a,2i) z^i")
    nok = 0
    a_max = 14
    for a in range(1, a_max + 1):
        cf = closed_form_P1(a)
        tr = cycle_poly_trace(1, a)
        nok += (cf == tr)
    record(f"closed form == trace for a = 1..{a_max} at P = 1", nok == a_max,
           f"{nok}/{a_max}")
    # explicit coefficient lists quoted in the frozen Blueprint section 2.4
    quoted = {1: [2], 2: [2, 2], 3: [2, 6], 4: [2, 12, 2], 5: [2, 20, 10],
              6: [2, 30, 30, 2]}
    nok = 0
    for a, exp in quoted.items():
        nok += (list(cycle_poly_trace(1, a)) == exp)
    record("coefficient lists match the frozen Blueprint section 2.4 verbatim "
           "(a = 1..6)", nok == len(quoted), f"{nok}/{len(quoted)}")


# ==========================================================================
# 5. Definition 5 -- the restricted cycle polynomial
# ==========================================================================
def check_restricted():
    hr("Definition 5  tilde W_a^{(P)} = tr(tilde T_P^a), tilde T_P = T_P o (J-I)")
    cases = [(1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (3, 3), (3, 4)]
    nok = 0
    for P, a in cases:
        bf = dict_to_poly(cycle_poly_bruteforce_restricted(P, a))
        tr = cycle_poly_trace_restricted(P, a)
        nok += (bf == tr)
        if bf != tr:
            print(f"    MISMATCH P={P} a={a}")
    record("restricted trace == brute force over words with no cyclically "
           f"adjacent equal exponents ({len(cases)} pairs)", nok == len(cases),
           f"{nok}/{len(cases)}")
    # the a = 1 special case
    nok = 0
    for P in range(0, 5):
        nok += (cycle_poly_trace_restricted(P, 1) == cycle_poly_trace(P, 1))
        nok += (zz.mat_trace(restricted_transfer_matrix(P)) == ())
    record("tilde W_1^{(P)} := W_1^{(P)} (Definition 5 special case; "
           "tr(tilde T_P) = 0)", nok == 10, f"{nok}/10")


# ==========================================================================
# 6. T1 -- the frozen hull-dimension identity
# ==========================================================================
def check_T1():
    hr("T1  dim hull_k(C) = sum_c ord_j(q) w_a^{(P)}(u^{(c)})   "
        "(polynomial route vs local-weight route)")
    cases = [
        (3, 2, 1, 3), (5, 2, 1, 7), (3, 3, 2, 15), (3, 4, 2, 21),
        (3, 4, 3, 5), (3, 4, 3, 17), (3, 4, 3, 25),
        (2, 4, 1, 15), (2, 4, 2, 15), (2, 4, 3, 15),
        (2, 8, 3, 7), (7, 2, 1, 5),
    ]
    total = 0
    agreed = 0
    skipped = 0
    mismatches = []
    for (p, e, k, n) in cases:
        F = FiniteField(p, e)
        for lam in {F.one, F.from_int(2) if p > 2 else F.one, F.gen}:
            res = verify_T1_on_instance(F, n, lam, k)
            for r in res:
                if "skipped" in r:
                    skipped += 1
                    continue
                total += 1
                if r["agree"]:
                    agreed += 1
                else:
                    mismatches.append(r)
    record(f"T1 verified on {total} (instance, word) checks",
           agreed == total and mismatches == [],
           f"{agreed}/{total} agree, {skipped} skipped (hypothesis (H) fails)")
    for m in mismatches[:5]:
        print("     ", m)
    return {"total": total, "agreed": agreed, "skipped": skipped}


# ==========================================================================
# 7. The a >= 3 example, independently constructed
# ==========================================================================
def check_a_ge_3():
    hr("a >= 3 -- independently constructed example, verified by brute-force "
       "enumeration of local states")
    out = {}

    # --- (i) purely combinatorial constructions -------------------------
    combos = [(1, 3), (1, 4), (2, 3), (2, 4), (3, 3), (3, 4)]
    nok = 0
    for P, a in combos:
        bf = dict_to_poly(cycle_poly_bruteforce(P, a))
        tr = cycle_poly_trace(P, a)
        nok += (bf == tr)
    record(f"a >= 3: trace == brute-force word enumeration ({len(combos)} "
           f"(P,a) pairs)", nok == len(combos), f"{nok}/{len(combos)}")

    # --- (ii) a REAL 4-cycle from the P1 implementation -----------------
    # F_81, n = 5, k = 3, lambda = 1 gives cycle shapes [(1,1),(4,1)]
    from core.cycles import e1_cycle_data
    F = FiniteField(3, 4)
    cd = e1_cycle_data(F, 5, F.one, 3)
    shapes = cd.cycle_shapes
    has_a4 = any(a >= 3 for a, d in shapes)
    record("a real instance with a cycle of length >= 3 was constructed from "
           "the P1 core", has_a4, f"F_81 n=5 k=3 lambda=1, shapes={shapes}")

    # the 4-cycle's polynomial, by brute force over its 2^4 = 16 local states
    P_, a4 = 1, 4
    bf = dict_to_poly(cycle_poly_bruteforce(P_, a4))
    tr = cycle_poly_trace(P_, a4)
    cf = closed_form_P1(a4)
    record(f"W_{a4}^{{(P={P_})}}(z) by local-state brute force == trace == closed form",
           bf == tr == cf,
           f"W = {zz.zformat(tr)}")

    # enumerate the 16 local states explicitly and show the histogram
    states = []
    for u in _words(P_, a4):
        states.append({"u": list(u), "w": local_weight(u, P_)})
    hist = {}
    for s in states:
        hist[s["w"]] = hist.get(s["w"], 0) + 1
    out["local_states_P1_a4"] = states
    out["histogram_P1_a4"] = hist

    # --- (iii) T1 on the real 4-cycle instance --------------------------
    res = verify_T1_on_instance(F, 5, F.one, 3)
    n_ok = sum(1 for r in res if r.get("agree"))
    record("T1 verified on the real 4-cycle instance "
           "(F_81, n=5, k=3, lambda=1)", n_ok == len(res),
           f"{n_ok}/{len(res)} words")
    out["T1_on_4cycle_instance"] = res

    # --- (iv) the weight of the 4-cycle word, by hand, for one word -----
    u = (1, 0, 1, 0)
    w = local_weight(u, 1)
    by_hand = sum(min(u[(m - 1) % 4], 1 - u[m]) for m in range(4))
    record("local weight computed by the cyclic formula matches the "
           "term-by-term evaluation on one a=4 word", w == by_hand,
           f"u={u}, w={w}")

    # --- (v) a >= 3 with P > 1 (repeated roots) -------------------------
    bf = dict_to_poly(cycle_poly_bruteforce(2, 3))
    tr = cycle_poly_trace(2, 3)
    record("a = 3 with P = 2 > 1 (repeated roots): brute force == trace",
           bf == tr, f"W = {zz.zformat(tr)}")
    out["W_P2_a3"] = zz.zformat(tr)
    return out


# ==========================================================================
# 8. The PA-1 specialisation boundary
# ==========================================================================
def check_PA1_boundary():
    hr("PA-1 specialisation and the novelty boundary")
    out = {}

    # --- (a) per-cycle objects, P = 1..7 (frozen section 2.0.1 table) ----
    ok1 = True
    ok2 = True
    for P in range(1, 8):
        W1 = cycle_poly_trace(P, 1)
        # PA-1: sum_b ||{b, P-b}|| z^b
        coeffs = [0] * (P // 2 + 1)
        for b in range(P // 2 + 1):
            coeffs[b] = len({b, P - b})
        if list(W1) != coeffs:
            ok1 = False
            print(f"    P={P}: W_1={zz.zformat(W1)} vs PA-1 {coeffs}")
        # PA-1: sum_a 2^{1-floor(a/P)} (a+1) z^a, over a = 0..P.
        # The range is 0..P because deg W_2^{(P)} = floor(2P/2) = P (T3).
        W2 = cycle_poly_trace(P, 2)
        pa = [(2 ** (1 - (a // P))) * (a + 1) for a in range(0, P + 1)]
        if list(W2) != pa:
            ok2 = False
            print(f"    P={P}: W_2={zz.zformat(W2)} vs PA-1 {pa}")
    record("W_1^{(P)} == PA-1's per-cycle object sum_b ||{b,P-b}|| z^b, "
           "P = 1..7", ok1)
    record("W_2^{(P)} == PA-1's per-cycle object "
           "sum_a 2^{1-floor(a/P)}(a+1) z^a, P = 1..7", ok2)

    # --- (b) the pointwise identity behind T1 ---------------------------
    # BLUEPRINT section 4.1, step 1, verbatim:
    #   P - max{u_m, P - u_{m-1}} = min{P - u_m, u_{m-1}}
    # which, since min is commutative, is exactly phi(u_{m-1}, u_m) =
    # min{u_{m-1}, P - u_m}, the frozen local weight of Definition 2.
    # NOTE: the pairing matters. max{u_m, P-u_{m-1}} must be paired with
    # min{P-u_m, u_{m-1}} -- NOT with min{u_m, P-u_{m-1}}, which is false.
    okp = True
    for P in range(1, 7):
        for um in range(0, P + 1):
            for um1 in range(0, P + 1):
                if P - max(um, P - um1) != min(P - um, um1):
                    okp = False
                if min(P - um, um1) != min(um1, P - um):
                    okp = False
    record("pointwise identity P - max{u_m, P-u_{m-1}} = min{P-u_m, u_{m-1}} "
           "= phi(u_{m-1},u_m)  (BLUEPRINT section 4.1, step 1)", okp)

    # --- (c) P = 1 reduction: what it collapses to -----------------------
    # W_1^{(1)} = 2 ; W_2^{(1)} = 2(1+z) ; hence
    #   N(y) = 2^B * prod_{c: a(c)=2} (1 + y^{d(c)})
    # which is exactly PA-1 Cor 12:  # = 2^{s+t} |h(l)|, |h| the coefficient of
    # X^l in Eq. (26), whose chi=0 factors collapse to 1 at P=1 (since
    # floor(p^nu/2) = 0) and whose chi=1 factors are (1 + X^{ord_j(q)}).
    w1 = cycle_poly_trace(1, 1)
    w2 = cycle_poly_trace(1, 2)
    record("P = 1: W_1^{(1)} = 2 and W_2^{(1)} = 2(1+z), so "
           "N(y) = 2^B prod_{a(c)=2} (1 + y^{d(c)})",
           w1 == (2,) and w2 == (2, 2),
           f"W_1={zz.zformat(w1)}, W_2={zz.zformat(w2)}")

    # structural check on a real PA-1-scope instance: k = 0 => a <= 2
    from core.cycles import e1_cycle_data
    nok = 0
    tested = 0
    for (p, e, n) in [(3, 2, 5), (3, 2, 7), (5, 2, 3), (3, 3, 5), (7, 2, 5)]:
        F = FiniteField(p, e)
        for lam in (F.one, F.from_int(2) if p > 2 else F.one):
            cd = e1_cycle_data(F, n, lam, 0)
            if not cd.permutation_ok:
                continue
            tested += 1
            nok += (max(cd.cycle_lengths) <= 2)
    record("for k = 0 every cycle has length <= 2 (PA-1's scope: "
           "f -> f* is an involution)", nok == tested and tested > 0,
           f"{nok}/{tested} instances")

    # and a >= 3 needs k not in {0, e/2}
    F = FiniteField(3, 4)
    lens = set()
    for n, lam in ((5, F.one), (17, F.one), (25, F.one)):
        cd = e1_cycle_data(F, n, lam, 3)
        lens.update(cd.cycle_lengths)
    record("for k = 3 (k not in {0, e/2}) cycle lengths >= 3 occur",
           max(lens) >= 3, f"lengths seen = {sorted(lens)}")
    out["cycle_lengths_k3_F81"] = sorted(lens)
    return out


# ==========================================================================
# 9. Edge cases
# ==========================================================================
def check_edge_cases():
    hr("Edge cases: a = 1, a = 2, P = 1, a >= 3")
    nok = 0
    tot = 0
    for (P, a) in [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (3, 1), (3, 2),
                   (3, 3), (4, 1), (4, 2), (4, 3), (5, 1)]:
        tot += 1
        bf = dict_to_poly(cycle_poly_bruteforce(P, a))
        tr = cycle_poly_trace(P, a)
        ok = (bf == tr and zz.zval(tr, 1) == (P + 1) ** a
              and zz.zcoeff(tr, 0) == 2
              and zz.zdeg(tr) == degree_bound(P, a))
        nok += ok
        if not ok:
            print(f"    MISMATCH P={P} a={a}")
    record(f"all edge cases satisfy trace == brute force, W(1)=(P+1)^a, "
           f"c_0=2, deg=floor(aP/2)", nok == tot, f"{nok}/{tot}")

    # P = 0 is degenerate (a single state) -- record rather than assert
    p0 = cycle_poly_trace(0, 3)
    record("P = 0 (degenerate, single state) is handled without error",
           p0 == (1,), f"W_3^{(0)} = {zz.zformat(p0)}")


# ==========================================================================
# 10. Independent SYMBOLIC verification with sympy
# ==========================================================================
def check_sympy_symbolic():
    hr("Independent symbolic verification (sympy, over Q[z] and Q[z,t])")
    try:
        import sympy as sp
    except ImportError:
        record("sympy available", False, "sympy not installed")
        return {}

    z = sp.Symbol("z")
    t = sp.Symbol("t")
    out = {}

    # (a) T3 symbolically: tr(T_P(z)^a) == our exact polynomial
    nok = 0
    cases = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
             (2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2), (3, 3),
             (4, 1), (4, 2), (5, 1)]
    for P, a in cases:
        n = P + 1
        T = sp.Matrix(n, n, lambda i, j: z ** min(i, P - j))
        theirs = sp.expand(sp.trace(T ** a))
        ours = sum(zz.zcoeff(cycle_poly_trace(P, a), i) * z ** i
                   for i in range(zz.zdeg(cycle_poly_trace(P, a)) + 1))
        nok += (sp.expand(theirs - ours) == 0)
    record(f"sympy confirms tr(T_P(z)^a) for {len(cases)} (P,a) pairs "
           "over Q[z]", nok == len(cases), f"{nok}/{len(cases)}")

    # (b) T4 symbolically: (1+sqrt z)^a + (1-sqrt z)^a.
    # Reduce in Q(w)[z] with w^2 = z: only even powers of w occur, so w^{2i}
    # becomes z^i and the result is a genuine Z[z] polynomial.
    w = sp.Symbol("w")
    nok = 0
    for a in range(1, 11):
        expr = sp.expand((1 + w) ** a + (1 - w) ** a)
        terms = sp.Poly(expr, w).terms()
        theirs = sum(c * z ** (e // 2) for (e,), c in terms)
        ours = sum(zz.zcoeff(closed_form_P1(a), i) * z ** i
                   for i in range(zz.zdeg(closed_form_P1(a)) + 1))
        nok += (sp.expand(theirs - ours) == 0)
    record("sympy confirms W_a^{(1)} = (1+sqrt z)^a + (1-sqrt z)^a, a = 1..10",
           nok == 10, f"{nok}/10")

    # (c) the bivariate GF: tr(tT(I-tT)^{-1}) = sum_{a>=1} W_a t^a
    nok = 0
    for P in (1, 2, 3):
        n = P + 1
        T = sp.Matrix(n, n, lambda i, j: z ** min(i, P - j))
        G = sp.expand(sp.trace(t * T * (sp.eye(n) - t * T).inv()))
        ser = sp.series(G, t, 0, 8).removeO()
        ser = sp.expand(ser)
        for a in range(1, 8):
            ours = sum(zz.zcoeff(cycle_poly_trace(P, a), i) * z ** i
                       for i in range(zz.zdeg(cycle_poly_trace(P, a)) + 1))
            coeff = sp.expand(ser.coeff(t, a))
            if sp.expand(coeff - ours) != 0:
                nok -= 1
                break
        else:
            nok += 1
    record("sympy confirms G_P(z,t) = tr(tT_P(I-tT_P)^{-1}) = "
           "sum_{a=1..7} W_a t^a for P = 1,2,3", nok == 3, f"{nok}/3")

    # (d) determinant degree, symbolically
    nok = 0
    for P in range(1, 6):
        n = P + 1
        T = sp.Matrix(n, n, lambda i, j: z ** min(i, P - j))
        d = sp.Poly(sp.expand((sp.eye(n) - t * T).det()), t)
        nok += (d.degree() == P + 1)
    record("sympy confirms deg_t det(I - t T_P(z)) = P+1 for P = 1..5",
           nok == 5, f"{nok}/5")

    # (e) the restricted matrix
    nok = 0
    for P, a in [(1, 3), (1, 4), (2, 3), (2, 4), (3, 3)]:
        n = P + 1
        T = sp.Matrix(n, n, lambda i, j: z ** min(i, P - j))
        JmI = sp.ones(n, n) - sp.eye(n)
        Tt = sp.matrix_multiply_elementwise(T, JmI)
        theirs = sp.expand(sp.trace(Tt ** a)) if a >= 2 else sp.expand(sp.trace(T ** a))
        ours = sum(zz.zcoeff(cycle_poly_trace_restricted(P, a), i) * z ** i
                   for i in range(zz.zdeg(cycle_poly_trace_restricted(P, a)) + 1))
        nok += (sp.expand(theirs - ours) == 0)
    record("sympy confirms tr(tilde T_P^a) for a >= 2 (Definition 5)",
           nok == 5, f"{nok}/5")

    # (f) characteristic polynomial -> recurrence, symbolically
    nok = 0
    for P in (1, 2, 3):
        n = P + 1
        T = sp.Matrix(n, n, lambda i, j: z ** min(i, P - j))
        chi = sp.Poly(sp.expand(T.charpoly().as_expr()), sp.Symbol("lambda"))
        # Cayley-Hamilton gives the recurrence; check it predicts W_{a+n}
        coeffs = [chi.coeff_monomial(sp.Symbol("lambda") ** (n - i))
                  for i in range(0, n + 1)]
        good = True
        for a in range(1, 6):
            acc = 0
            for i in range(0, n + 1):
                Wi = cycle_poly_trace(P, a + n - i)
                acc += coeffs[i] * sum(zz.zcoeff(Wi, j) * z ** j
                                       for j in range(zz.zdeg(Wi) + 1))
            if sp.expand(acc) != 0:
                good = False
                break
        nok += good
    record("sympy confirms the Cayley-Hamilton recurrence for P = 1,2,3",
           nok == 3, f"{nok}/3")
    return out


# ==========================================================================
# main
# ==========================================================================
def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    trace_rows = check_trace_vs_bruteforce()
    check_t3_corollaries()
    check_t5()
    check_t4()
    check_restricted()
    t1 = check_T1()
    a3 = check_a_ge_3()
    pa1 = check_PA1_boundary()
    check_edge_cases()
    check_sympy_symbolic()

    elapsed = time.perf_counter() - t0

    # machine-readable CSV of every (P, a) polynomial computed
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        import csv
        w = csv.writer(fh)
        w.writerow(["P", "a", "poly_bruteforce", "poly_trace", "agree",
                    "deg", "deg_bound", "constant_term", "value_at_1"])
        for P in range(1, 6):
            for a in range(1, 13):
                if (P + 1) ** a > 3_000_000:
                    continue
                bf = dict_to_poly(cycle_poly_bruteforce(P, a))
                tr = cycle_poly_trace(P, a)
                w.writerow([P, a, zz.zformat(bf), zz.zformat(tr), bf == tr,
                            zz.zdeg(tr), degree_bound(P, a),
                            zz.zcoeff(tr, 0), zz.zval(tr, 1)])

    n_pass = sum(1 for _n, ok, _d in checks if ok)
    n_fail = sum(1 for _n, ok, _d in checks if not ok)

    hr("P2 SUMMARY")
    print(f"  checks : {n_pass} passed, {n_fail} failed   ({elapsed:.1f} s)")
    for name, ok, detail in checks:
        if not ok:
            print(f"    FAILED: {name}  {detail}")

    summary.update({
        "elapsed_seconds": round(elapsed, 2),
        "checks_total": n_pass + n_fail,
        "checks_passed": n_pass,
        "checks_failed": n_fail,
        "checks": [{"name": n, "pass": ok, "detail": d} for n, ok, d in checks],
        "T1": t1,
        "a_ge_3": a3,
        "PA1_boundary": pa1,
        "verdict": "PASS" if n_fail == 0 else "FAIL",
    })
    OUT_JSON.write_text(json.dumps(summary, indent=2, default=str),
                        encoding="utf-8")
    print(f"\n  summary -> {OUT_JSON.relative_to(ROOT)}")
    print(f"  detail  -> {OUT_CSV.relative_to(ROOT)}")
    print(f"\n  P2 VERDICT: {summary['verdict']}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
