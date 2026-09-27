# P1 CLOSURE AUDIT — status correction (no new research)

**Phase:** P1 (of `EXECUTION_PLAN.md` @ `96c7885`)
**Frozen source:** `BLUEPRINT.md` @ `076f910` (FINAL) — **untouched**
**Supersedes:** the status wording in `docs/P1_REPORT.md` @ commit `1a70750`
**Executed:** 2026-09-27
**Branch:** `arena/01a0e2d5-arun-galios`

This document records a **targeted correction to the P1 status metadata only**.
No new computation, no new research, no re-run of the 5 670-instance suite, and
no change to the frozen research direction.

---

## 1. Why this correction was needed

The P1 report as first committed (`1a70750`) described OPEN-22 as
"**RESOLVED by redefinition**". That wording overstated the position: the
original inherited 385-instance regression suite was **not** recovered, and
cannot be, because its exact $\lambda$ and $k$ assignments and its
odd-characteristic $q$ and $n$ assignments are absent from the frozen
Blueprint. Replacing the suite is not the same as resolving the open item.

The correction below narrows the claim to what is actually true.

---

## 2. Corrected status (authoritative)

| Item | Status |
|---|---|
| **P1 Core Implementation** | **PASS** |
| **New 5 670-instance validation suite** | **PASS** |
| **Original 385-instance regression reproducibility** | **OPEN / NOT RECONSTRUCTIBLE** |
| **OPEN-28** | **retained** as the provenance label for newly chosen odd-$q$ inputs |

**OPEN-22 is NOT resolved.** It remains OPEN.

---

## 3. What is absent from the frozen Blueprint (unchanged finding)

| Missing field | Scope | Consequence |
|---|---|---|
| $\lambda$ | all 385 instances | no instance can be named |
| $k$ | all 385 instances | no instance can be named; and $k$ is not even an input to the two A.7.1 checks |
| $q$ | the 176 odd instances | the odd half cannot be enumerated at all |
| $n$ | the 176 odd instances | the odd half cannot be enumerated at all |
| per-instance expected results | all 385 instances | only the aggregate "385 instances, 0 failures" is recorded |

The only concrete odd instance named anywhere in the Blueprint is $q=3,n=9$
(the $\nu>e$ crash example).

Reconstructing $\lambda$ or $k$ from the totals 209/176 would be inference from
a reported total — forbidden. Inventing them is forbidden. Neither was done.
**No instance count, $\lambda$, $k$, $q$, $n$ or expected result was invented.**

---

## 4. Instance categories — kept strictly separate

Three categories exist and are **never mixed**. A result from one category is
never reported as a result from another.

| Category | Definition | Provenance | Status |
|---|---|---|---|
| **(a) Inherited / probe-derived** | computed during blueprinting by the throwaway prototype | `BLUEPRINT.md` Annex A — counts and outcomes only, no per-instance data | **NOT RECONSTRUCTIBLE / OPEN-22**. Never used as a P1 acceptance gate. |
| **(b) Newly constructed validation** | defined in P1 by an explicit closed-form rule over stated $q$/$n$ grids and all of $\mathbb F_q^*$ | `regression_385.yaml` → `regression_suite_v1`; odd-$q$ inputs labelled **OPEN-28** | **PASS 5 670 / 5 670** |
| **(c) Published-paper benchmark** | expected outcome printed in a published paper | P1 Ex. 4 / Ex. 5 / Ex. 6, P3 Ex. 3.6 (Blueprint §5.5, CP1) | **PASS 4 / 4** |

The 5 670 / 5 670 figure is described **only** as the newly defined validation
suite (category b). It is never described as recovery, reconstruction, or
reproduction of the original 385-instance suite (category a).

Cross-checks against category (a) — the 7 §5.5 cycle instances re-agreeing with
the A.3 [PROBE-DONE] reference — are reported as **corroboration of the new
implementation**, not as reproduction of the 385-instance suite. The A.3
reference values remain **[PROBE-DONE]** and are not paper results.

---

## 5. Retained OPEN items (explicitly preserved)

| ID | Item | Status |
|---|---|---|
| **OPEN-22** | Original 385-instance suite not reconstructible | **OPEN / NOT RECONSTRUCTIBLE** |
| **OPEN-28** | Provenance label for the newly chosen odd-$q$ inputs in `regression_suite_v1` | **retained** |
| **OPEN-29** | CP1-3 / CP1-4 are **shape**-comparison only: the source papers' primitive elements $\omega\in\mathbb F_{25},\mathbb F_{81}$ are not recoverable from the PDFs (extracted text mangles the field-definition displays), and $\lambda$ is only determined up to Frobenius conjugacy. Two of four CP1 checks are shape-level, not coefficient-level — which matches how `BLUEPRINT.md` A.1 records them ("identical shape"). | **OPEN** (owner P4) |
| **OPEN-30** | `galois` cannot factor $x^{21}+1$ over $\mathbb F_{16}$ / $\mathbb F_{256}$ (`RuntimeError` after 1000 tries), while our factoriser succeeds and is certified intrinsically. 2 oracle tests skipped. | **OPEN** — oracle limitation, not a defect of ours |

---

## 6. Files changed by this closure audit

| File | Change |
|---|---|
| `docs/P1_REPORT.md` | status header block; §B.4 rewritten as OPEN-22 STATUS; new §B.5 instance-provenance categories; §C.2 heading; §F results table split by category; §G restructured (G.1 closed / G.2 still-open / G.3 carried); §H.1 row 5–6; §H.2 verdict; closure-audit note in header |
| `galois_hull_enum/instances/regression_385.yaml` | new top-level `status` block recording the four status items and the three provenance categories |
| `galois_hull_enum/results/phase1_status.json` | new machine-readable status record |
| `docs/P1_CLOSURE_AUDIT.md` | this document |

**Unchanged:** all of `core/`, `tests/`, `scripts/`, and every results artefact
from the executed runs — `results/phase1_factorisation.csv` (5 670 rows),
`results/phase1_summary.json`, `results/regression_manifest.csv`,
`results/phase1_tests.log`, `results/phase1_regression.log`,
`results/phase1_pilot_estimate.json`, `results/independence_record.json`.
The 5 670-instance suite was **not re-run**; no number was altered.

---

## 7. Frozen-state verification

| File | sha256 | vs freeze commit |
|---|---|---|
| `BLUEPRINT.md` | `2e02572a0ef970a93dcf242c725f290938d9b180787beaccd8acedf6a609214f` | identical to `076f910` ✅ |
| `EXECUTION_PLAN.md` | `6026536e3df10ac80a35225ddc42fc3b5f52f845be92da98531aca4bc8f65b32` | identical to `96c7885` ✅ |
| `Galios Hull.zip` | `7d503417be98539096df819dfe1258786dc72e27d209dbfa4716309fefe98985` | unchanged ✅ |

---

## 8. Verdict

> ## **P1 Core Implementation: PASS**
> ## **New 5 670-instance validation suite: PASS**
> ## **Original 385-instance regression reproducibility: OPEN / NOT RECONSTRUCTIBLE**
> ## **OPEN-28: retained**
>
> The P1 implementation and the newly defined validation suite stand as PASS.
> The original 385-instance suite remains OPEN and is not claimed as
> reproduced, recovered, or reconstructed. Instance categories (a), (b) and (c)
> are kept strictly separate, and the 5 670 / 5 670 figure is attributed only
> to category (b). OPEN-29 and OPEN-30 are retained and documented.
>
> **STOPPING here. P2 has not been started.**
