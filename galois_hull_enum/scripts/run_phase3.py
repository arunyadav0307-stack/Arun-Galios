#!/usr/bin/env python3
r"""P3 (CP3) driver -- end-to-end enumerator vs independent brute force.

BLUEPRINT.md section 6.2, CP3:

    For every instance in section 5.5 with |C| <= 10^5: E2 output == E3
    histogram, and |C| = (P+1)^{sum a}, and #{dim = 0} = 2^B.
    Gate: >= 300 verified instances across the X14 sweep with zero mismatches.

Sub-commands (run in this order; the estimate is committed BEFORE the full run,
per the execution plan's pilot -> estimate -> checkpoint -> full-run rule):

    pilot      P3a  -- the small pilot, covering a=1, a=2, a>=3, P=1, P>1,
                       even and odd characteristic
    estimate   P3b  -- cost estimate derived from the pilot's measured
                       timings, plus checkpoint boundaries and the definition
                       of a complete run
    full       P3c  -- the frozen section 5.5 grid (X1-X10) plus the X14 sweep
    novelty    P3d  -- the sub-suite inaccessible to PA-1: genuine a >= 3
                       cycles, repeated-root P > 1, and general lambda

Every instance is run through THREE routes:

    E2             the enumerator (transfer matrix; core/cycle_poly factors)
    oracle_cycle   E3 as the Blueprint writes it, indexed by the cycle product
    oracle_flat    E3 indexed by a flat exponent vector, no cycle code at all

All three must agree on the COMPLETE distribution, coefficient by coefficient.
Exact integer arithmetic throughout; no floating point anywhere in the
mathematical path.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from p3lib import (GRID, MAX_CODES, RING_ROWS, X14_FIELDS, Instance,          # noqa: E402
                   element_of_order, grid_instances, peak_rss_mb, sweep_instances)

RESULTS = os.path.join(ROOT, "results")


# --------------------------------------------------------------------------
def compare(e2: Dict[int, int], oc: Dict[int, int], of: Dict[int, int]
            ) -> Dict[str, object]:
    """Coefficient-level comparison of the three distributions."""
    keys = sorted(set(e2) | set(oc) | set(of))
    rows = []
    mismatches = []
    for d in keys:
        a, b, c = e2.get(d, 0), oc.get(d, 0), of.get(d, 0)
        rows.append({"dim": d, "e2": a, "oracle_cycle": b, "oracle_flat": c})
        if not (a == b == c):
            mismatches.append({"dim": d, "e2": a, "oracle_cycle": b,
                               "oracle_flat": c,
                               "max_abs_diff": max(abs(a - b), abs(a - c),
                                                   abs(b - c))})
    total_codes = sum(r["e2"] for r in rows)
    return {
        "coefficients": rows,
        "n_coefficients_compared": len(rows),
        "n_codes_compared": total_codes,
        "n_mismatched_coefficients": len(mismatches),
        "mismatches": mismatches,
        "max_abs_coefficient_diff": max((m["max_abs_diff"] for m in mismatches),
                                        default=0),
        "complete_agreement": not mismatches,
    }


def structural(cd, e2: Dict[int, int], of_stats: Dict[str, object]) -> Dict[str, object]:
    """The CP3 structural predicates plus the code-count cross-checks."""
    from core import zz
    from enumeration.enumerator import enumerator_poly
    poly = enumerator_poly(cd.cycle_shapes, cd.P)
    sum_a = cd.sum_a
    n_factors = len(cd.factors)
    return {
        "sum_a": sum_a,
        "n_factors": n_factors,
        "sum_a_equals_n_factors": sum_a == n_factors,
        "n_codes_E2_N1": zz.zval(poly, 1),
        "n_codes_predicted_E1": cd.predicted_code_count,
        "n_codes_oracle_flat": of_stats.get("n_codes"),
        "n_codes_match": (zz.zval(poly, 1) == cd.predicted_code_count
                          == of_stats.get("n_codes")),
        "lcd_count_E2": poly[0] if poly else 0,
        "lcd_count_oracle": e2.get(0, 0),
        "lcd_count_2^B": 2 ** cd.B,
        "lcd_match": (poly[0] if poly else 0) == 2 ** cd.B,
        "all_coeffs_nonneg_int": all(isinstance(c, int) and c >= 0 for c in poly),
        "support": sorted(e2),
    }


def run_instance(inst: Instance, crosscheck: bool = True) -> Dict[str, object]:
    """Run one instance through all three routes and assemble the record."""
    e2, oc, of, timings = inst.run_all(crosscheck=crosscheck)
    cmp = compare(e2, oc, of)
    st = structural(inst.cd, e2, timings["oracle_flat_stats"])
    return {
        **inst.summary_dict(),
        "key": inst.key,
        "e2_distribution": {str(k): v for k, v in sorted(e2.items())},
        "oracle_cycle_distribution": {str(k): v for k, v in sorted(oc.items())},
        "oracle_flat_distribution": {str(k): v for k, v in sorted(of.items())},
        "comparison": cmp,
        "structural": st,
        "timings_seconds": {k: (round(v, 6) if isinstance(v, float) else v)
                            for k, v in timings.items()},
        "verified": bool(cmp["complete_agreement"] and st["n_codes_match"]
                         and st["lcd_match"] and st["all_coeffs_nonneg_int"]
                         and st["sum_a_equals_n_factors"]),
    }


# --------------------------------------------------------------------------
def write_outputs(records: List[Dict[str, object]], prefix: str) -> None:
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, f"{prefix}_instances.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "q", "p", "e", "n", "k", "lam", "r", "P", "nu",
                    "n_prime", "j", "B", "sum_a", "n_factors", "shapes",
                    "n_codes", "n_coefficients", "n_codes_compared",
                    "max_dim", "lcd", "verified", "t_e2", "t_oracle_cycle",
                    "t_oracle_flat"])
        for r in records:
            w.writerow([r["label"], r["q"], r["p"], r["e"], r["n"], r["k"],
                        r["lam"], r["r"], r["P"], r["nu"], r["n_prime"],
                        r["j"], r["B"], r["sum_a"], r["n_factors"],
                        "|".join(f"{a}x{d}" for a, d in r["shapes"]),
                        r["structural"]["n_codes_E2_N1"],
                        r["comparison"]["n_coefficients_compared"],
                        r["comparison"]["n_codes_compared"],
                        max((int(d) for d in r["e2_distribution"]), default=0),
                        r["structural"]["lcd_count_E2"], r["verified"],
                        r["timings_seconds"]["e2"],
                        r["timings_seconds"]["oracle_cycle"],
                        r["timings_seconds"]["oracle_flat"]])

    with open(os.path.join(RESULTS, f"{prefix}_coefficients.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "dim", "e2", "oracle_cycle", "oracle_flat",
                    "agree"])
        for r in records:
            for row in r["comparison"]["coefficients"]:
                w.writerow([r["label"], row["dim"], row["e2"],
                            row["oracle_cycle"], row["oracle_flat"],
                            row["e2"] == row["oracle_cycle"] == row["oracle_flat"]])


def aggregate(records: List[Dict[str, object]]) -> Dict[str, object]:
    verified = [r for r in records if r["verified"]]
    t_e2 = sum(r["timings_seconds"]["e2"] for r in records)
    t_oc = sum(r["timings_seconds"]["oracle_cycle"] for r in records)
    t_of = sum(r["timings_seconds"]["oracle_flat"] for r in records)
    return {
        "n_instances_run": len(records),
        "n_verified": len(verified),
        "n_failed": len(records) - len(verified),
        "n_coefficients_compared": sum(r["comparison"]["n_coefficients_compared"]
                                       for r in records),
        "n_codes_compared": sum(r["comparison"]["n_codes_compared"]
                                for r in records),
        "n_mismatched_coefficients": sum(
            r["comparison"]["n_mismatched_coefficients"] for r in records),
        "max_abs_coefficient_difference": max(
            (r["comparison"]["max_abs_coefficient_diff"] for r in records),
            default=0),
        "structural_failures": [r["label"] for r in records
                                if not r["structural"]["n_codes_match"]
                                or not r["structural"]["lcd_match"]
                                or not r["structural"]["all_coeffs_nonneg_int"]
                                or not r["structural"]["sum_a_equals_n_factors"]],
        "seconds": {"e2": round(t_e2, 3), "oracle_cycle": round(t_oc, 3),
                    "oracle_flat": round(t_of, 3),
                    "total": round(t_e2 + t_oc + t_of, 3)},
        "peak_rss_mb": round(peak_rss_mb(), 1),
    }


# --------------------------------------------------------------------------
def cmd_pilot(args) -> int:
    """P3a -- the small pilot."""
    print("=" * 78)
    print("P3a  PILOT -- enumerator vs two independent oracles")
    print("=" * 78)

    # The pilot is chosen to hit every regime the brief names, and is a strict
    # subset of the frozen section 5.5 grid (no invented parameters).
    rows = [r for r in GRID if r["id"] in
            ("X1", "X2", "X3", "X4", "X5", "X6", "X8a", "X8c", "X9")]
    insts = [i for i in grid_instances(rows) if i.prepare()]
    # X8 with n=15 needs lambda=1 (r=1): the Blueprint's |C| = 2^12 for that
    # cell forces r=1, since r=3 violates gcd(n',r)=1 there.
    F16 = insts[0].F.__class__(2, 4)
    insts.append(Instance(F16, 15, F16.one, 3, label="X8b ",
                          source="new (even q, l=4)"))
    insts[-1].prepare()

    print(f"\npilot instances: {len(insts)}")
    print(f"{'label':<7}{'q':<5}{'n':<4}{'k':<3}{'r':<4}{'P':<3}"
          f"{'shapes':<26}{'|C|':<8}{'E2=C':<6}{'E2=F':<6}{'t_orF':<8}")
    records = []
    for inst in insts:
        t0 = time.time()
        rec = run_instance(inst)
        rec["_wall"] = round(time.time() - t0, 4)
        records.append(rec)
        print(f"{rec['label']:<7}{rec['q']:<5}{rec['n']:<4}{rec['k']:<3}"
              f"{rec['r']:<4}{rec['P']:<3}"
              f"{'|'.join(f'{a}x{d}' for a, d in rec['shapes']):<26}"
              f"{rec['structural']['n_codes_E2_N1']:<8}"
              f"{'OK' if rec['e2_distribution'] == rec['oracle_cycle_distribution'] else 'BAD':<6}"
              f"{'OK' if rec['e2_distribution'] == rec['oracle_flat_distribution'] else 'BAD':<6}"
              f"{rec['timings_seconds']['oracle_flat']:<8}")

    agg = aggregate(records)
    print(f"\n--- pilot aggregate ---")
    print(f"instances run      : {agg['n_instances_run']}")
    print(f"verified           : {agg['n_verified']}")
    print(f"failed             : {agg['n_failed']}")
    print(f"coefficients cmp   : {agg['n_coefficients_compared']}")
    print(f"codes compared     : {agg['n_codes_compared']}")
    print(f"mismatched coeffs  : {agg['n_mismatched_coefficients']}")
    print(f"seconds            : {agg['seconds']}")
    print(f"peak RSS (MB)      : {agg['peak_rss_mb']}")

    # regime coverage check -- the brief's explicit pilot requirements
    regimes = {
        "a=1": any(1 in r["cycle_lengths"] for r in records),
        "a=2": any(2 in r["cycle_lengths"] for r in records),
        "a>=3": any(max(r["cycle_lengths"]) >= 3 for r in records),
        "P=1": any(r["P"] == 1 for r in records),
        "P>1": any(r["P"] > 1 for r in records),
        "odd characteristic": any(r["p"] % 2 == 1 for r in records),
        "even characteristic": any(r["p"] == 2 for r in records),
    }
    print("\n--- pilot regime coverage ---")
    for kname, v in regimes.items():
        print(f"  {kname:<24}: {'covered' if v else 'NOT COVERED'}")

    write_outputs(records, "phase3_pilot")
    with open(os.path.join(RESULTS, "phase3_pilot.json"), "w") as fh:
        json.dump({"phase": "P3a pilot", "records": records,
                   "aggregate": agg, "regimes": regimes,
                   "max_codes_bound": MAX_CODES}, fh, indent=1)
    print(f"\nwrote results/phase3_pilot.json, phase3_pilot_instances.csv, "
          f"phase3_pilot_coefficients.csv")
    print(f"\nP3a PILOT VERDICT: "
          f"{'PASS' if agg['n_failed'] == 0 and all(regimes.values()) else 'FAIL'}")
    return 0 if agg["n_failed"] == 0 and all(regimes.values()) else 1


# --------------------------------------------------------------------------
def cmd_estimate(args) -> int:
    """P3b -- cost estimate from the measured pilot."""
    print("=" * 78)
    print("P3b  COST ESTIMATE (derived from measured pilot timings)")
    print("=" * 78)
    path = os.path.join(RESULTS, "phase3_pilot.json")
    if not os.path.exists(path):
        print("P3a pilot not found; run `pilot` first.")
        return 1
    with open(path) as fh:
        pilot = json.load(fh)

    recs = pilot["records"]
    codes = [r["structural"]["n_codes_E2_N1"] for r in recs]
    t_of = [r["timings_seconds"]["oracle_flat"] for r in recs]
    t_oc = [r["timings_seconds"]["oracle_cycle"] for r in recs]
    t_e2 = [r["timings_seconds"]["e2"] for r in recs]
    ncoef = [r["comparison"]["n_coefficients_compared"] for r in recs]

    # per-code and per-coefficient rates, from the largest pilots (least noise)
    pairs = sorted(zip(codes, t_of))
    per_code = [t / c for c, t in pairs if c >= 32] or [t_of[-1] / max(codes[-1], 1)]
    rate_code = max(per_code)                       # conservative upper bound
    rate_coef = max((t / max(k, 1) for t, k in zip(t_of, ncoef)), default=0.0)
    rate_oc = max((t / max(c, 1) for t, c in zip(t_oc, codes)), default=0.0)
    rate_e2 = max(t_e2) if t_e2 else 0.0

    # how many instances does the sweep actually contain?
    print("\nenumerating the X14 sweep (this also measures sweep-construction cost)...")
    t0 = time.time()
    insts, rejected = sweep_instances()
    t_sweep = time.time() - t0
    codes_list = [i.cd.predicted_code_count for i in insts]

    total_codes = sum(codes_list)
    est_of = total_codes * rate_code
    est_oc = total_codes * rate_oc
    est_e2 = len(insts) * rate_e2
    est_total = est_of + est_oc + est_e2 + t_sweep

    # checkpoint boundaries: split the sweep into ~5 chunks of equal cost
    order = sorted(range(len(insts)),
                   key=lambda i: -insts[i].cd.predicted_code_count)
    chunk_cost = est_of / 5.0
    chunks: List[List[int]] = [[] for _ in range(5)]
    acc, ci = 0.0, 0
    for i in order:
        if ci < 4 and acc + insts[i].cd.predicted_code_count * rate_code > chunk_cost:
            ci += 1
            acc = 0.0
        chunks[ci].append(i)
        acc += insts[i].cd.predicted_code_count * rate_code

    estimate = {
        "measured_from_pilot": {
            "n_pilot_instances": len(recs),
            "pilot_code_counts": codes,
            "pilot_oracle_flat_seconds": t_of,
            "pilot_oracle_cycle_seconds": t_oc,
            "pilot_e2_seconds": t_e2,
            "per_code_seconds_conservative": rate_code,
            "per_coefficient_seconds_conservative": rate_coef,
            "per_code_oracle_cycle_seconds_conservative": rate_oc,
            "per_instance_e2_seconds_conservative": rate_e2,
        },
        "sweep": {
            "n_candidate_instances": len(insts),
            "n_rejected": len(rejected),
            "total_codes": total_codes,
            "max_codes_per_instance": max(codes_list) if codes_list else 0,
            "sweep_construction_seconds": round(t_sweep, 3),
        },
        "estimate_seconds": {
            "oracle_flat": round(est_of, 1),
            "oracle_cycle": round(est_oc, 1),
            "e2": round(est_e2, 1),
            "sweep_construction": round(t_sweep, 1),
            "total": round(est_total, 1),
            "total_minutes": round(est_total / 60.0, 1),
        },
        "checkpoint_boundaries": [
            {"checkpoint": f"CP3-{i + 1}", "n_instances": len(chunks[i]),
             "codes": sum(insts[j].cd.predicted_code_count for j in chunks[i]),
             "estimated_seconds": round(sum(insts[j].cd.predicted_code_count
                                            for j in chunks[i]) * rate_code, 1)}
            for i in range(5) if chunks[i]
        ],
        "complete_run_definition": (
            "A COMPLETE RUN is: every instance produced by the frozen X14 sweep "
            "filter (q in the frozen field list, n <= 60, gcd(n',r)=1, "
            "r | (1+p^(e-k)), |C| <= 10^5) is run through E2, oracle_cycle and "
            "oracle_flat; all three distributions agree on EVERY coefficient; "
            "|C| = (P+1)^sum_a = the flat oracle's code count; "
            "#{dim=0} = 2^B; all coefficients are nonnegative integers. "
            "Instances with |C| > 10^5 are NOT verified -- they are recorded as "
            "enumerator-only, exactly as the frozen Blueprint's validation-scope "
            "amendment requires. The run is complete only when the recorded "
            "n_verified equals the number of sweep instances and n_failed = 0."
        ),
        "verdict": "PROCEED" if est_total < 6 * 3600 else "TOO SLOW -- narrow the sweep",
    }

    print(f"\nsweep candidates          : {len(insts)}  (rejected {len(rejected)})")
    print(f"total codes to certify    : {total_codes}")
    print(f"sweep construction        : {t_sweep:.1f} s")
    print(f"conservative per-code rate: {rate_code * 1e6:.1f} us/code")
    print(f"\nESTIMATE (conservative, from measured pilot):")
    for kname, v in estimate["estimate_seconds"].items():
        print(f"  {kname:<22}: {v}")
    print(f"\ncheckpoint boundaries:")
    for cb in estimate["checkpoint_boundaries"]:
        print(f"  {cb['checkpoint']}: {cb['n_instances']} instances, "
              f"{cb['codes']} codes, ~{cb['estimated_seconds']} s")
    print(f"\nVERDICT: {estimate['verdict']}")

    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "phase3_estimate.json"), "w") as fh:
        json.dump(estimate, fh, indent=1)
    with open(os.path.join(RESULTS, "phase3_sweep_rejected.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["q", "n", "k", "lam", "r", "reason"])
        for r in rejected:
            w.writerow([r["q"], r["n"], r["k"], r["lam"], r["r"], r["reason"]])
    print(f"\nwrote results/phase3_estimate.json, phase3_sweep_rejected.csv")
    return 0 if estimate["verdict"] == "PROCEED" else 1


# --------------------------------------------------------------------------
def cmd_full(args) -> int:
    """P3c -- the frozen grid plus the X14 sweep."""
    print("=" * 78)
    print("P3c  FULL VALIDATION -- frozen section 5.5 grid + X14 sweep")
    print("=" * 78)

    records: List[Dict[str, object]] = []

    # ---- part 1: the frozen grid rows X1-X10 ----
    print("\n--- frozen grid X1-X10 ---")
    insts = [i for i in grid_instances(GRID) if i.prepare()]
    from core.field import FiniteField
    F16 = FiniteField(2, 4)
    insts.append(Instance(F16, 15, F16.one, 3, label="X8b ",
                          source="new (even q, l=4)"))
    insts[-1].prepare()

    # The frozen E3 feasibility bound is |C| <= 10^5. Rows above it (X10a, X10b)
    # are computed by E2 ONLY and flagged, exactly as the Blueprint's
    # validation-scope amendment requires. They are deliberately NOT pushed
    # through the oracle: that would be infeasible and would silently breach the
    # frozen feasibility bound.
    over = [i for i in insts if i.cd.predicted_code_count > MAX_CODES]
    insts = [i for i in insts if i.cd.predicted_code_count <= MAX_CODES]

    for inst in insts:
        rec = run_instance(inst)
        records.append(rec)
        _print_row(rec)

    flagged = []
    for inst in over:
        from enumeration.enumerator import enumerator_poly, structural_checks
        poly = enumerator_poly(inst.cd.cycle_shapes, inst.cd.P)
        flagged.append({
            **inst.summary_dict(),
            "status": "ENUMERATOR ONLY -- |C| exceeds the frozen E3 bound "
                      f"{MAX_CODES}; NOT verified against an oracle",
            "e2_distribution": {str(i): c for i, c in enumerate(poly) if c},
            "structural": structural_checks(inst.cd.cycle_shapes, inst.cd.P,
                                            inst.cd.predicted_code_count,
                                            inst.cd.B),
        })
        print(f"{inst.label:<7}{inst.F.q:<5}{inst.n:<4}{inst.k:<3}"
              f"{inst.F.multiplicative_order(inst.lam):<4}{inst.cd.P:<3}"
              f"{'|'.join(f'{a}x{d}' for a, d in inst.cd.cycle_shapes):<26}"
              f"{inst.cd.predicted_code_count:<8}ENUMERATOR ONLY")

    # ---- part 2: the X14 sweep (checkpointed, resumable) ----
    print(f"\n--- X14 sweep ---")
    ckpt = os.path.join(RESULTS, "phase3_sweep_checkpoint.jsonl")
    done_keys = set()
    if args.resume and os.path.exists(ckpt):
        with open(ckpt) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    done_keys.add(json.loads(line)["key"])
        print(f"resuming: {len(done_keys)} instances already checkpointed")

    t0 = time.time()
    sweep, rejected = sweep_instances()
    print(f"sweep produced {len(sweep)} certifiable instances "
          f"({len(rejected)} rejected) in {time.time() - t0:.1f} s")

    t_start = time.time()
    ck = open(ckpt, "a")
    n_new = 0
    for idx, inst in enumerate(sweep, 1):
        if inst.key in done_keys:
            continue
        rec = run_instance(inst)
        records.append(rec)
        ck.write(json.dumps(rec) + "\n")
        ck.flush()
        os.fsync(ck.fileno())
        n_new += 1
        if n_new % 25 == 0:
            el = time.time() - t_start
            rate = el / max(n_new, 1)
            print(f"  [{n_new} new / {len(sweep) - len(done_keys)} to do] "
                  f"verified so far {sum(1 for r in records if r['verified'])}  "
                  f"elapsed {el:.0f}s  ({rate:.2f}s/instance, "
                  f"eta {rate * (len(sweep) - len(done_keys) - n_new) / 60:.0f} min)  "
                  f"peak RSS {peak_rss_mb():.0f} MB", flush=True)
    ck.close()
    if args.resume and done_keys:
        # fold the previously checkpointed records back in
        with open(ckpt) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    r = json.loads(line)
                    rk = r.get("key")
                    if rk is not None and rk not in {x.get("key") for x in records}:
                        records.append(r)

    agg = aggregate(records)
    _print_aggregate(agg)

    # novelty sub-suite statistics (P3d reports on the same records)
    nov = novelty_stats(records)

    write_outputs(records, "phase3_full")
    with open(os.path.join(RESULTS, "phase3_summary.json"), "w") as fh:
        json.dump({
            "phase": "P3c full validation",
            "records": records,
            "aggregate": agg,
            "novelty_subset": nov,
            "enumerator_only_flagged": flagged,
            "n_sweep_rejected": len(rejected),
            "checkpoint_file": "results/phase3_sweep_checkpoint.jsonl",
            "ring_rows_deferred_to_CP6": RING_ROWS,
            "max_codes_bound": MAX_CODES,
        }, fh, indent=1)
    with open(os.path.join(RESULTS, "phase3_full_rejected.csv"), "w",
              newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["q", "n", "k", "lam", "r", "reason"])
        for r in rejected:
            w.writerow([r["q"], r["n"], r["k"], r["lam"], r["r"], r["reason"]])

    print(f"\nwrote results/phase3_summary.json, phase3_full_instances.csv, "
          f"phase3_full_coefficients.csv, phase3_full_rejected.csv")
    ok = agg["n_failed"] == 0 and agg["n_verified"] >= 300
    print(f"\nP3c FULL-VALIDATION VERDICT: "
          f"{'PASS' if ok else 'FAIL'} "
          f"({agg['n_verified']} verified, {agg['n_failed']} failed)")
    return 0 if ok else 1


# --------------------------------------------------------------------------
def novelty_stats(records: List[Dict[str, object]]) -> Dict[str, object]:
    """The P3d sub-suite: the region PA-1 cannot reach."""
    def sel(pred):
        return [r for r in records if pred(r)]

    a_ge3 = sel(lambda r: max(r["cycle_lengths"]) >= 3)
    p_gt1 = sel(lambda r: r["P"] > 1)
    both = sel(lambda r: max(r["cycle_lengths"]) >= 3 and r["P"] > 1)
    general_lam = sel(lambda r: r["r"] not in (1, 2))
    return {
        "n_instances_total": len(records),
        "a_ge_3": {"n": len(a_ge3),
                   "n_verified": sum(1 for r in a_ge3 if r["verified"]),
                   "max_cycle_length": max((max(r["cycle_lengths"])
                                            for r in a_ge3), default=0),
                   "codes": sum(r["structural"]["n_codes_E2_N1"] for r in a_ge3)},
        "P_gt_1": {"n": len(p_gt1),
                   "n_verified": sum(1 for r in p_gt1 if r["verified"]),
                   "P_values": sorted({r["P"] for r in p_gt1}),
                   "codes": sum(r["structural"]["n_codes_E2_N1"] for r in p_gt1)},
        "a_ge_3_and_P_gt_1": {"n": len(both),
                              "n_verified": sum(1 for r in both if r["verified"]),
                              "codes": sum(r["structural"]["n_codes_E2_N1"]
                                           for r in both)},
        "general_lambda_r_not_1_or_2": {
            "n": len(general_lam),
            "n_verified": sum(1 for r in general_lam if r["verified"]),
            "r_values": sorted({r["r"] for r in general_lam}),
            "codes": sum(r["structural"]["n_codes_E2_N1"] for r in general_lam)},
    }


def cmd_novelty(args) -> int:
    """P3d -- the novelty sub-suite, reported from the full run's records."""
    print("=" * 78)
    print("P3d  NOVELTY SUB-SUITE -- the region inaccessible to PA-1")
    print("=" * 78)
    path = os.path.join(RESULTS, "phase3_summary.json")
    if not os.path.exists(path):
        print("P3c full run not found; run `full` first.")
        return 1
    with open(path) as fh:
        data = json.load(fh)
    records = data["records"]
    nov = novelty_stats(records)

    print(f"\ninstances in the full run : {nov['n_instances_total']}")
    print(f"{'regime':<34}{'n':>6}{'verified':>11}{'codes':>12}")
    for kname in ("a_ge_3", "P_gt_1", "a_ge_3_and_P_gt_1",
                  "general_lambda_r_not_1_or_2"):
        d = nov[kname]
        print(f"{kname:<34}{d['n']:>6}{d['n_verified']:>11}{d['codes']:>12}")
    print(f"\nmax cycle length seen with a>=3 : "
          f"{nov['a_ge_3']['max_cycle_length']}")
    print(f"P values with P>1              : {nov['P_gt_1']['P_values']}")
    print(f"r values with r not in {{1,2}}   : "
          f"{nov['general_lambda_r_not_1_or_2']['r_values']}")

    with open(os.path.join(RESULTS, "phase3_novelty.json"), "w") as fh:
        json.dump(nov, fh, indent=1)
    ok = (nov["a_ge_3"]["n_verified"] == nov["a_ge_3"]["n"] > 0
          and nov["a_ge_3_and_P_gt_1"]["n_verified"]
          == nov["a_ge_3_and_P_gt_1"]["n"] > 0)
    print(f"\nP3d NOVELTY SUB-SUITE VERDICT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


# --------------------------------------------------------------------------
def _print_row(rec: Dict[str, object]) -> None:
    agree = (rec["e2_distribution"] == rec["oracle_cycle_distribution"]
             == rec["oracle_flat_distribution"])
    print(f"{rec['label']:<7}{rec['q']:<5}{rec['n']:<4}{rec['k']:<3}"
          f"{rec['r']:<4}{rec['P']:<3}"
          f"{'|'.join(f'{a}x{d}' for a, d in rec['shapes']):<26}"
          f"{rec['structural']['n_codes_E2_N1']:<8}"
          f"{'OK' if agree else 'BAD':<6}"
          f"{'OK' if rec['verified'] else 'BAD'}")


def _print_aggregate(agg: Dict[str, object]) -> None:
    print(f"\n--- aggregate ---")
    print(f"instances run       : {agg['n_instances_run']}")
    print(f"verified            : {agg['n_verified']}")
    print(f"failed              : {agg['n_failed']}")
    print(f"coefficients cmp    : {agg['n_coefficients_compared']}")
    print(f"codes compared      : {agg['n_codes_compared']}")
    print(f"mismatched coeffs   : {agg['n_mismatched_coefficients']}")
    print(f"max |coeff diff|    : {agg['max_abs_coefficient_difference']}")
    print(f"structural failures : {agg['structural_failures']}")
    print(f"seconds             : {agg['seconds']}")
    print(f"peak RSS (MB)       : {agg['peak_rss_mb']}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("pilot", help="P3a pilot")
    sub.add_parser("estimate", help="P3b cost estimate")
    pf = sub.add_parser("full", help="P3c full validation")
    pf.add_argument("--resume", action="store_true",
                    help="resume from results/phase3_sweep_checkpoint.jsonl")
    sub.add_parser("novelty", help="P3d novelty sub-suite report")
    args = ap.parse_args(argv)
    return {"pilot": cmd_pilot, "estimate": cmd_estimate,
            "full": cmd_full, "novelty": cmd_novelty}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
