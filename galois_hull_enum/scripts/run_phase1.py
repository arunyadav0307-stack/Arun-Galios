#!/usr/bin/env python3
"""P1 -- run the factorisation regression suite and the CP1 quartet gate.

Executes ONLY P1-sized computations:
  * factorisation of x^n - lambda over F_q  (Algorithm E1, steps 1-3)
  * the two A.7.1 checks:  product identity, and irreducibility of every
    reported factor
  * structural invariants: n = n' p^nu, P = p^nu, mu^P = lambda,
    gcd(n', p) = 1, every factor monic with multiplicity exactly P

No cycle polynomials, no enumerators, no hull dimensions, no moments
(those are P2/P3/P6).

Long-run discipline (EXECUTION_PLAN.md): pilot -> estimate -> checkpoint ->
full run.

Usage:
    python3 scripts/run_phase1.py --pilot 60     # pilot only, prints estimate
    python3 scripts/run_phase1.py --limit 500    # first 500 instances
    python3 scripts/run_phase1.py                # the full suite
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core.field import FiniteField  # noqa: E402
from core import poly as P  # noqa: E402
from core.factor import factor_xn_minus_lambda, is_irreducible, v_p  # noqa: E402

RESULTS = ROOT / "results"
MANIFEST = RESULTS / "regression_manifest.csv"
OUT_CSV = RESULTS / "phase1_factorisation.csv"
OUT_JSON = RESULTS / "phase1_summary.json"
CHECKPOINT_EVERY = 500


# --------------------------------------------------------------------------
# the CP1 quartet gate (BLUEPRINT.md A.1 / section 6.2 CP1)
# --------------------------------------------------------------------------
def cp1_quartet():
    """Reproduce the four published factorisations.

    CP1-1 (F_9,  n=3,  lambda=2)  -> (x+1)^3                          EXACT
    CP1-2 (F_27, n=15, lambda=2)  -> (x+1)^3 (x^4+2x^3+x^2+2x+1)^3    EXACT
    CP1-3 (F_25, n=7,  lambda=g^8)-> degrees [1,3,3]                  SHAPE
    CP1-4 (F_81, n=21, lambda=g^8)-> degrees [1,3,3]                  SHAPE

    CP1-1/CP1-2 are compared STRING-EXACTLY against the printed source
    polynomials (all their coefficients lie in the prime field, so the
    comparison carries no primitive-element convention). CP1-3/CP1-4 are
    compared by SHAPE, which is what BLUEPRINT.md A.1 records ("identical
    shape") and is Frobenius-invariant.
    """
    rows = []

    F = FiniteField(3, 2)
    _np, _nu, _Pw, _mu, fac = factor_xn_minus_lambda(F, 3, F.from_int(2))
    got = [(P.format_poly(F, f), m) for f, m in fac]
    rows.append(("CP1-1", "P1 Ex.4", "F_9", 3, "2", got, [("x+1", 3)], "EXACT",
                 got == [("x+1", 3)]))

    F = FiniteField(3, 3)
    _np, _nu, _Pw, _mu, fac = factor_xn_minus_lambda(F, 15, F.from_int(2))
    got = [(P.format_poly(F, f), m) for f, m in fac]
    exp = [("x+1", 3), ("x^4+2x^3+x^2+2x+1", 3)]
    rows.append(("CP1-2", "P1 Ex.6", "F_27", 15, "2", got, exp, "EXACT",
                 got == exp))

    F = FiniteField(5, 2)
    lam = F.gen_pow(8)
    _np, _nu, _Pw, _mu, fac = factor_xn_minus_lambda(F, 7, lam)
    got = sorted(P.deg(f) for f, m in fac)
    rows.append(("CP1-3", "P1 Ex.5", "F_25", 7, "g^8 (ord 3)", got, [1, 3, 3],
                 "SHAPE", got == [1, 3, 3]))

    F = FiniteField(3, 4)
    lam = F.gen_pow(8)
    _np, _nu, _Pw, _mu, fac = factor_xn_minus_lambda(F, 21, lam)
    got = sorted(P.deg(f) for f, m in fac)
    rows.append(("CP1-4", "P3 Ex.3.6", "F_81", 21, "g^8 (ord 10)", got,
                 [1, 3, 3], "SHAPE", got == [1, 3, 3]))

    return rows


# --------------------------------------------------------------------------
# one regression instance
# --------------------------------------------------------------------------
def run_instance(F, q, p, e, n, lam):
    """Returns a result dict; `ok` is the conjunction of every check."""
    t0 = time.perf_counter()
    rec = {
        "q": q, "p": p, "e": e, "n": n, "lambda": lam,
        "nu": None, "n_prime": None, "P": None, "mu": None,
        "n_factors": None, "degrees": None, "multiplicities": None,
        "ok_product": False, "ok_irreducible": False, "ok_structure": False,
        "ok_mu": False, "seconds": None, "error": "",
    }
    try:
        nu = v_p(n, p)
        n_prime = n // (p ** nu)
        Pw = p ** nu
        n_p, nu2, Pw2, mu, fac = factor_xn_minus_lambda(F, n, lam)

        # structural invariants
        rec["nu"], rec["n_prime"], rec["P"], rec["mu"] = nu2, n_p, Pw2, mu
        struct = (nu2 == nu and n_p == n_prime and Pw2 == Pw
                  and n_prime * (p ** nu) == n and n_prime % p != 0)
        rec["ok_structure"] = bool(struct)

        # mu^P = lambda   (DEFECT (b) check)
        rec["ok_mu"] = bool(F.powi(mu, Pw2) == lam)

        # product identity:  prod f_i^{m_i} == x^n - lambda
        prod = [1]
        for f, m in fac:
            for _ in range(m):
                prod = P.pmul(F, prod, f)
        target = [F.neg(lam)] + [0] * (n - 1) + [1]
        rec["ok_product"] = bool(P.trim(prod) == P.trim(target))

        # every factor irreducible, monic, with multiplicity exactly P
        ok_irr = True
        for f, m in fac:
            if m != Pw2 or P.plead(F, f) != 1 or not is_irreducible(F, f):
                ok_irr = False
                break
        rec["ok_irreducible"] = bool(ok_irr)

        rec["n_factors"] = len(fac)
        rec["degrees"] = "|".join(str(P.deg(f)) for f, m in fac)
        rec["multiplicities"] = "|".join(str(m) for f, m in fac)
    except Exception as exc:                     # noqa: BLE001
        rec["error"] = f"{type(exc).__name__}: {exc}"
    rec["seconds"] = round(time.perf_counter() - t0, 6)
    rec["ok"] = bool(rec["ok_product"] and rec["ok_irreducible"]
                     and rec["ok_structure"] and rec["ok_mu"] and not rec["error"])
    return rec


def load_manifest():
    with MANIFEST.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=int, default=0,
                    help="run only the first N instances and print an estimate")
    ap.add_argument("--limit", type=int, default=0,
                    help="run only the first N instances (no estimate)")
    args = ap.parse_args()

    rows = load_manifest()
    total = len(rows)
    if args.pilot:
        rows = rows[:args.pilot]
    elif args.limit:
        rows = rows[:args.limit]

    print("=" * 74)
    print("P1 -- factorisation regression (Algorithm E1 steps 1-3)")
    print("=" * 74)
    print(f"manifest          : {MANIFEST.relative_to(ROOT)}")
    print(f"instances in file : {total}")
    print(f"instances to run  : {len(rows)}"
          f"{'  (PILOT)' if args.pilot else ('  (LIMITED)' if args.limit else '  (FULL)')}")

    RESULTS.mkdir(parents=True, exist_ok=True)

    fields_cache = {}
    out_rows = []
    t_start = time.perf_counter()
    for i, row in enumerate(rows, start=1):
        q = int(row["q"]); p = int(row["p"]); e = int(row["e"])
        n = int(row["n"]); exp = int(row["lambda_exponent"])
        if q not in fields_cache:
            fields_cache[q] = FiniteField(p, e)
        F = fields_cache[q]
        lam = F.gen_pow(exp)

        rec = run_instance(F, q, p, e, n, lam)
        rec["grid_id"] = row["grid_id"]
        rec["lambda_exponent"] = exp
        rec["char"] = row["char"]
        rec["structural_category"] = row["structural_category"]
        rec["defects"] = row["defects"]
        out_rows.append(rec)

        if i % CHECKPOINT_EVERY == 0:
            elapsed = time.perf_counter() - t_start
            nok = sum(1 for r in out_rows if r["ok"])
            rate = i / elapsed
            print(f"  [checkpoint] {i}/{len(rows)}  ok={nok}  "
                  f"elapsed={elapsed:7.1f}s  rate={rate:7.1f} inst/s  "
                  f"eta={(len(rows)-i)/rate if rate else 0:7.1f}s")
            _write_csv(out_rows)

    elapsed = time.perf_counter() - t_start
    _write_csv(out_rows)

    n_ok = sum(1 for r in out_rows if r["ok"])
    n_fail = len(out_rows) - n_ok
    print(f"\nran {len(out_rows)} instances in {elapsed:.1f}s "
          f"({len(out_rows)/elapsed:.1f} inst/s)")
    print(f"  PASS: {n_ok}    FAIL: {n_fail}")

    if args.pilot:
        rate = len(out_rows) / elapsed
        print(f"\n--- PILOT ESTIMATE ---")
        print(f"  pilot rate        : {rate:.1f} instances/s")
        print(f"  projected FULL    : {total / rate:.1f} s "
              f"for {total} instances")
        print("  (estimate committed to results/phase1_pilot_estimate.json;")
        print("   the full run may now proceed)")
        (RESULTS / "phase1_pilot_estimate.json").write_text(json.dumps({
            "pilot_instances": len(out_rows),
            "pilot_seconds": round(elapsed, 3),
            "rate_instances_per_second": round(rate, 3),
            "total_instances": total,
            "projected_full_seconds": round(total / rate, 2),
        }, indent=2), encoding="utf-8")
        return 0 if n_fail == 0 else 1

    # ---- breakdown ----
    def brk(key, val):
        sub = [r for r in out_rows if r[key] == val]
        return len(sub), sum(1 for r in sub if r["ok"])

    print("\n  breakdown by characteristic:")
    for v in ("even", "odd"):
        tot, ok = brk("char", v)
        print(f"    {v:5s}: {ok}/{tot} pass")
    print("  breakdown by defect coverage:")
    for d in ("defect_a_even_char", "defect_b_nu_gt_e"):
        sub = [r for r in out_rows if d in (r["defects"] or "")]
        print(f"    {d:22s}: {sum(1 for r in sub if r['ok'])}/{len(sub)} pass")
    rep = [r for r in out_rows if r["P"] > 1]
    print(f"    repeated roots (P>1)    : {sum(1 for r in rep if r['ok'])}/{len(rep)} pass")

    # ---- CP1 quartet ----
    print("\n" + "=" * 74)
    print("CP1 gate -- the four published factorisations")
    print("=" * 74)
    cp1 = cp1_quartet()
    for cid, src, fld, n, lam_s, got, exp, mode, ok in cp1:
        print(f"  {cid} {src:12s} {fld:5s} n={n:<3} lambda={lam_s:14s} "
              f"[{mode}] {'PASS' if ok else 'FAIL'}")
        print(f"      got      : {got}")
        print(f"      expected : {exp}")
    cp1_ok = all(r[-1] for r in cp1)

    summary = {
        "phase": "P1",
        "instances_run": len(out_rows),
        "instances_available": total,
        "passed": n_ok,
        "failed": n_fail,
        "elapsed_seconds": round(elapsed, 2),
        "rate_instances_per_second": round(len(out_rows) / elapsed, 2),
        "checks": {
            "product_identity": sum(1 for r in out_rows if r["ok_product"]),
            "all_factors_irreducible": sum(1 for r in out_rows if r["ok_irreducible"]),
            "structure_nu_nprime_P": sum(1 for r in out_rows if r["ok_structure"]),
            "mu_to_the_P": sum(1 for r in out_rows if r["ok_mu"]),
        },
        "by_characteristic": {v: {"total": brk("char", v)[0], "passed": brk("char", v)[1]}
                              for v in ("even", "odd")},
        "errors": [r["error"] for r in out_rows if r["error"]][:20],
        "cp1_quartet": {
            "all_pass": bool(cp1_ok),
            "detail": [{"id": r[0], "source": r[1], "field": r[2], "n": r[3],
                        "lambda": r[4], "got": str(r[5]), "expected": str(r[6]),
                        "mode": r[7], "pass": bool(r[8])} for r in cp1],
        },
    }
    OUT_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nsummary -> {OUT_JSON.relative_to(ROOT)}")
    print(f"detail  -> {OUT_CSV.relative_to(ROOT)}")

    verdict = (n_fail == 0) and cp1_ok
    print(f"\nP1 REGRESSION VERDICT: {'PASS' if verdict else 'FAIL'}")
    return 0 if verdict else 1


def _write_csv(rows):
    if not rows:
        return
    cols = ["grid_id", "q", "p", "e", "n", "lambda_exponent", "lambda", "char",
            "nu", "n_prime", "P", "mu", "n_factors", "degrees",
            "multiplicities", "ok_structure", "ok_mu", "ok_product",
            "ok_irreducible", "ok", "seconds", "error",
            "structural_category", "defects"]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
