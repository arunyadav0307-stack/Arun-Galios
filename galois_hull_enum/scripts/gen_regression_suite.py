#!/usr/bin/env python3
"""Materialise the P1 regression instance list from regression_385.yaml.

Deterministic: no randomness, no dependence on the host. Expands the
`RANGE(a,b)` shorthand used in the YAML into explicit lambda exponent lists
and writes results/regression_manifest.csv.

It also VALIDATES every declared field of every (q, n) entry against values
computed from first principles (nu = v_p(n), n_prime = n/p^nu, P = p^nu),
so that a typo in the YAML cannot silently produce a wrong suite.

Usage:  python3 scripts/gen_regression_suite.py
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core.factor import v_p  # noqa: E402

YAML_PATH = ROOT / "instances" / "regression_385.yaml"
OUT_DIR = ROOT / "results"
OUT_CSV = OUT_DIR / "regression_manifest.csv"

_RANGE_RE = re.compile(r"^RANGE\(\s*(\d+)\s*,\s*(\d+)\s*\)$")


def expand_exponents(spec):
    """Accept an explicit list, or the RANGE(a,b) shorthand (inclusive)."""
    if isinstance(spec, list):
        return list(spec)
    m = _RANGE_RE.match(str(spec).strip())
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        return list(range(a, b + 1))
    raise ValueError(f"cannot parse lambda_exponents spec: {spec!r}")


def factor_q(q: int):
    """(p, e) with q = p^e, p prime."""
    for p in range(2, int(q ** 0.5) + 2):
        if q % p == 0:
            e = 0
            m = q
            while m % p == 0:
                m //= p
                e += 1
            if m == 1:
                return p, e
    return q, 1


def main() -> int:
    doc = yaml.safe_load(YAML_PATH.read_text(encoding="utf-8"))
    suite = doc["regression_suite_v1"]
    entries = suite["instances"]

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    problems = []
    for ent in entries:
        q = ent["q"]
        n = ent["n"]
        p, e = factor_q(q)
        if p ** e != q:
            problems.append(f"{ent['id']}: q={q} is not a prime power")
            continue
        nu = v_p(n, p)
        n_prime = n // (p ** nu)
        P = p ** nu

        # validate every declared derived field
        for key, computed in (("p", p), ("e", e), ("nu", nu),
                              ("n_prime", n_prime), ("P", P)):
            declared = ent.get(key)
            if declared != computed:
                problems.append(
                    f"{ent['id']}: declared {key}={declared} but computed {computed}"
                )
        declared_char = ent.get("char")
        computed_char = "even" if p == 2 else "odd"
        if declared_char != computed_char:
            problems.append(
                f"{ent['id']}: declared char={declared_char} but computed {computed_char}"
            )

        exps = expand_exponents(ent["lambda_exponents"])
        if len(exps) != ent["lambda_count"]:
            problems.append(
                f"{ent['id']}: lambda_count={ent['lambda_count']} "
                f"but {len(exps)} exponents given"
            )
        if exps != list(range(q - 1)):
            problems.append(
                f"{ent['id']}: exponents are not the complete 0..q-2 set "
                f"(lambda rule says ALL of F_q^*)"
            )

        has_nu_gt_e = nu > e
        declared_defects = ent.get("defects") or []
        if has_nu_gt_e and "defect_b_nu_gt_e" not in declared_defects:
            problems.append(
                f"{ent['id']}: nu={nu} > e={e} but defect_b_nu_gt_e not declared"
            )
        if not has_nu_gt_e and "defect_b_nu_gt_e" in declared_defects:
            problems.append(
                f"{ent['id']}: defect_b_nu_gt_e declared but nu={nu} <= e={e}"
            )
        if p == 2 and "defect_a_even_char" not in declared_defects:
            problems.append(f"{ent['id']}: p=2 but defect_a_even_char not declared")

        for i in exps:
            rows.append({
                "grid_id": ent["id"],
                "q": q, "p": p, "e": e, "n": n,
                "nu": nu, "n_prime": n_prime, "P": P,
                "char": computed_char,
                "lambda_exponent": i,
                "lambda": None,          # filled by run_phase1.py
                "structural_category": ent["category"],
                "defects": "|".join(declared_defects),
                "provenance": ent["provenance"],
            })

    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_grid = len(entries)
    n_even = sum(1 for r in rows if r["char"] == "even")
    n_odd = sum(1 for r in rows if r["char"] == "odd")
    n_nu = sum(1 for r in rows if r["nu"] > r["e"])
    n_rep = sum(1 for r in rows if r["P"] > 1)

    print(f"suite            : regression_suite_v1")
    print(f"grid entries     : {n_grid}   ((q,n) pairs)")
    print(f"TOTAL INSTANCES  : {len(rows)}")
    print(f"  even character.: {n_even}")
    print(f"  odd character. : {n_odd}")
    print(f"  nu > e         : {n_nu}   (defect (b))")
    print(f"  repeated roots : {n_rep}   (P > 1)")
    print(f"written          : {OUT_CSV.relative_to(ROOT)}")

    if problems:
        print(f"\nVALIDATION PROBLEMS ({len(problems)}):")
        for pr in problems:
            print("  -", pr)
        return 1
    print("\nYAML self-validation: OK (every declared field matches computation)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
