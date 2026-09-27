# P1 REPORT — Independent Core Implementation

**Phase:** P1 (of `EXECUTION_PLAN.md` @ `96c7885`)
**Frozen source:** `BLUEPRINT.md` @ `076f910` (FINAL)
**Builds on:** P0 @ `076f977` / tag `p0`
**Executed:** 2026-09-27
**Branch:** `arena/01a0e2d5-arun-galios`
**Scope:** `galois_hull_enum/core/` — fields, polynomials, factorisation, the
frozen $f^{\#}$ map and Algorithm E1 cycle detection. **No P2+ work.**
**Closure audit applied:** metadata/status correction only (see
`docs/P1_CLOSURE_AUDIT.md`). No computation re-run, no research added.

> ## P1 STATUS (authoritative)
>
> | Item | Status |
> |---|---|
> | **P1 Core Implementation** | **PASS** |
> | **New 5 670-instance validation suite** | **PASS** |
> | **Original 385-instance regression reproducibility** | **OPEN / NOT RECONSTRUCTIBLE** |
> | **OPEN-28** | **retained** as the provenance label for newly chosen odd-$q$ inputs |
>
> **OPEN-22 is NOT resolved.** The original inherited 385-instance regression
> suite remains unreconstructible: its exact $\lambda$ and $k$ assignments, and
> its odd-characteristic $q$ and $n$ assignments, are absent from the frozen
> Blueprint. The 5 670 / 5 670 result below is the **newly defined** validation
> suite and is **not** a recovery of the original 385 suite. See §B.

---

## A. P1 IMPLEMENTATION REPORT

### A.1 What was built

| Module | Lines | Contents | Blueprint basis |
|---|---|---|---|
| `core/field.py` | 370 | $\mathbb F_{p^e}$: log/antilog + Zech-tables, digitwise add, Frobenius tables, order arithmetic | §3.1 (standing data) |
| `core/poly.py` | 270 | Dense polynomials over $\mathbb F_q$: add/mul/divmod/gcd/egcd/powmod/deriv/eval, Frobenius on coefficients and on $\mathbb F_q[x]/(f)$ | §3.1 ($f^\#$), §5.1 (E1) |
| `core/factor.py` | 395 | Musser squarefree decomposition, DDF, EDF, irreducibility, and $x^n-\lambda$ | §5.1 (E1 steps 1–3), A.7.1 (defects a,b) |
| `core/cycles.py` | 215 | The frozen $f^{\#}$ map, and E1 step 3–4 cycle detection with $(a,d)$ shapes | §3.1 ($f^\#$), §5.1 (E1 steps 3–4) |
| `tests/test_field.py` | 150 | field kernel | — |
| `tests/test_poly.py` | 175 | polynomial kernel | — |
| `tests/test_factor.py` | 285 | factorisation + **both defect regressions** | A.7.1 |
| `tests/test_cycles.py` | 285 | $f^\#$, E1 identities, A.3 cross-check | §5.1, A.3 |
| `tests/test_oracle.py` | 175 | independent `galois` oracle | P1 brief §4 |
| `tests/test_independence.py` | 145 | machine-checked independence proof | P1 brief §2/§3 |
| `scripts/gen_regression_suite.py` | 175 | materialises the instance manifest, self-validates the YAML | OPEN-22 |
| `scripts/run_phase1.py` | 260 | pilot → estimate → checkpoint → full run | plan §6.2 |

**Total: ~2 800 lines of implementation and tests.**

### A.2 Encoding decision, and why it matters for CP1

Elements are encoded as $a=\sum c_i p^i$, i.e. the base-$p$ digits of the
encoding **are** the coefficient vector in the basis $(1,\alpha,\dots,\alpha^{e-1})$.
Two deliberate consequences:

- $0$ is the zero element, $1$ the one element;
- **an integer $c\in[0,p)$ encodes the prime-subfield element $c$ directly.**

The second point is what makes **CP1-1 and CP1-2 exact rather than
shape-only**: the source papers' $\mathbb F_9$ and $\mathbb F_{27}$ polynomials
have all their coefficients in $\mathbb F_3$, so the comparison carries **no
primitive-element convention at all**. For CP1-3/CP1-4 the source's primitive
element ($\omega$ in $\mathbb F_{25}$, $\mathbb F_{81}$) is **not recoverable
from the PDFs** (the extracted text mangles the field-definition displays), so
those two are compared by **shape**, which is exactly what
`BLUEPRINT.md` A.1 records ("identical shape") and which is
Frobenius-invariant. **No λ was invented to close this gap.**

### A.3 A defect found in my own factoriser, and fixed

`test_factor_poly_is_correct` over $\mathbb F_2$ caught a genuine bug:
**Yun's squarefree-decomposition algorithm is invalid in small
characteristic.** On $f=x^4(x+1)^3$ it returned $[(x{+}1,1)]$ instead of
$[(x,4),(x{+}1,3)]$, because Yun's multiplicity counter is an ordinary integer
that gets reduced mod $p$. This is precisely the class of failure the
Blueprint's A.7.1 warns about.

**Fix:** replaced Yun with **Musser's algorithm**, which tracks multiplicities
by repeated gcd against the $\gcd(f,f')$ accumulator and recurses on $p$-th
roots — valid for small characteristic. Verified on the exact failing case and
on $(x+1)^6$ over $\mathbb F_3$ (now correctly reports multiplicity 6 = 2·3).

This was found by P1's own unit tests, **not** reported as a Blueprint defect.

---

## B. `regression_385.yaml` AUDIT (OPEN-22)

### B.1 The critical finding

**The inherited 385-instance suite is NOT reconstructible from the frozen
Blueprint.** I searched the entire frozen text. What it actually records:

| Half | Count | $q$ | $n$ | $\lambda$ | $k$ | expected results |
|---|---|---|---|---|---|---|
| even | 209 | ✅ $\{2,4,8,16,32,64,256\}$ | ✅ 10 values | ❌ **nowhere** | ❌ **nowhere** | ❌ aggregate only |
| odd | 176 | ❌ **nowhere** | ❌ **nowhere** | ❌ **nowhere** | ❌ **nowhere** | ❌ aggregate only |

The Blueprint's *only* concrete odd instance is $q=3,n=9$ (the $\nu>e$ crash
example). Reconstructing $\lambda$ and $k$ from the totals 209/176 would be
**inference from a reported total**, which is forbidden, and inventing them is
forbidden. So I did neither.

### B.2 What the file therefore contains

The file has **three** clearly separated sections:

1. **`inherited_probe_suite` — status `OPEN`.** A faithful audit of the
   385-instance suite, with each unrecoverable field explicitly `value: null`
   and a `reason:` explaining exactly why it cannot be recovered. Marked
   `epistemic_class: [PROBE-DONE]` and `reproducible_from_this_repo: false`.
   **Nothing inferred, nothing invented.**

2. **`regression_suite_v1` — status `DEFINED`.** The suite P1 actually runs,
   fully specified by a closed-form rule:
   - $q\in\{2,3,4,5,8,9,16,25,27,32,49,64,81,256\}$ (the 7 even values are
     **recorded** in A.7.1; the 7 odd values are **newly defined here and
     labelled as such**, chosen so that $\nu>e$, $\nu=e$ and $\nu<e$ all occur
     inside the recorded $n$-grid);
   - $n\in\{3,5,7,9,15,17,21,25,31,35\}$ — **recorded** in A.7.1;
   - $\lambda$ = **all** of $\mathbb F_q^*$, enumerated as $g^i$, $i=0..q-2$;
   - $k$ = `NOT_AN_INPUT`, because both A.7.1 checks (product identity,
     irreducibility) are $k$-independent. $k$ is exercised separately, with
     frozen values, in `cycle_suite`;
   - every derived field ($p,e,\nu,n',P$, character, defect coverage) is
     declared **and cross-checked against first principles** by the generator;
   - `expected_results_asserted: false` — expectations are **computed**, never
     asserted from a reported aggregate.

3. **`cycle_suite` — status `DEFINED`.** Seven instances with **frozen**
   $k,j,e-k$ and the standing hypotheses (H1)/(H2) checked, taken from
   `BLUEPRINT.md` §5.5 X1–X7.

### B.3 Resulting suite size

```
grid entries     : 140   ((q,n) pairs)
TOTAL INSTANCES  : 5670
  even character.: 3750
  odd characteristic. : 1920
  nu > e         : 6   (defect (b))
  repeated roots : 720   (P > 1)
```

**This is 5 670, not 385, and was not tuned towards 385.** It is a superset in
*coverage* (14 fields, all $\lambda$, all defect cases) chosen by an explicit
rule, because the inherited suite cannot be recovered.

### B.4 OPEN-22 STATUS — **OPEN / NOT RECONSTRUCTIBLE** (not resolved)

> **The original 385-instance regression suite is NOT recovered and NOT
> reconstructible.** It stays OPEN. What was absent, and remains absent:
>
> | Missing field | Scope | Why unrecoverable |
> |---|---|---|
> | $\lambda$ | all 385 instances | recorded nowhere in the frozen Blueprint |
> | $k$ | all 385 instances | recorded nowhere; and not even an input to the two A.7.1 checks |
> | $q$ | the 176 odd instances | no $q$ value recorded for any of them |
> | $n$ | the 176 odd instances | no $n$ value recorded for any of them |
> | expected results | all 385 instances | only the aggregate "0 failures" is recorded |
>
> Reconstructing $\lambda$ or $k$ from the totals 209/176 would be inference
> from a reported total, which is forbidden; inventing them is forbidden. So
> neither was done. **No instance count, $\lambda$, $k$, $q$, $n$ or expected
> result was invented, and none was inferred from the reported totals.**
>
> P1's validation gate is the **newly defined** `regression_suite_v1`, which is
> a *replacement* for the unreconstructible suite, **not** a recovery of it.
> The odd $q$-grid inside it is a *choice of test inputs* by me, retained under
> the provenance label **OPEN-28**, and is not a recovery of the inherited odd
> half.

### B.5 Instance provenance categories — kept strictly separate

Three categories of instance exist in this phase. **They are never mixed, and a
result from one category is never reported as a result from another.**

| Category | Definition | Provenance | Where recorded | Status |
|---|---|---|---|---|
| **(a) Inherited / probe-derived** | instances computed during blueprinting by the throwaway prototype | `BLUEPRINT.md` Annex A — counts and outcomes only, no per-instance data | `regression_385.yaml` → `inherited_probe_suite` (audit, fields `null`) | **NOT RECONSTRUCTIBLE / OPEN-22**; never used as a P1 acceptance gate |
| **(b) Newly constructed validation** | instances defined in P1 by an explicit closed-form rule over stated $q$/$n$ grids and all of $\mathbb F_q^*$ | `regression_385.yaml` → `regression_suite_v1`; odd-$q$ inputs labelled **OPEN-28** | `results/regression_manifest.csv`, `results/phase1_factorisation.csv` | **PASS 5 670 / 5 670** |
| **(c) Published-paper benchmark** | instances whose expected outcome is printed in a published paper | P1 Ex. 4 / Ex. 5 / Ex. 6, P3 Ex. 3.6 (Blueprint §5.5, CP1) | §C.3 of this report | **PASS 4 / 4** |

Cross-checks against category (a) — the 7 §5.5 cycle instances re-agreeing with
the A.3 [PROBE-DONE] reference (§C.6) — are reported as **corroboration of the
new implementation**, never as a reproduction of the 385-instance suite.

The A.3 reference values themselves remain **[PROBE-DONE]**: real sandbox
computations, but not paper results, not reproducible from this repository, and
not to be quoted in the manuscript until regenerated by committed code.

---

## C. TEST / VERIFICATION REPORT

### C.1 Full test suite — **375 passed, 5 skipped, 0 failed**

```
tests/test_field.py         78 passed
tests/test_poly.py          73 passed
tests/test_factor.py        78 passed   (incl. both defect regressions)
tests/test_cycles.py        50 passed, 3 skipped (k out of range for F_2/F_3)
tests/test_oracle.py        90 passed, 2 skipped (galois oracle limitation)
tests/test_independence.py   6 passed
                             -------------------------------
                             375 passed, 5 skipped, 0 failed   (53.6 s)
```

The 3 `k out of range` skips are parametrisation artefacts ($k>e-1$); the 2
oracle skips are documented below.

### C.2 Newly defined validation suite — **5 670 / 5 670 PASS**

> This is `regression_suite_v1` (category (b) of §B.5), defined in P1 by an
> explicit rule. It is **not** the original 385-instance suite (category (a)),
> which remains **OPEN / NOT RECONSTRUCTIBLE** (OPEN-22), and it is **not** a
> recovery of it.

```
PASS: 5670    FAIL: 0        (77.0 s, 73.7 inst/s)
  by q: 2(10) 3(20) 4(30) 5(40) 8(70) 9(80) 16(150) 25(240) 27(260)
        32(310) 49(480) 64(630) 81(800) 256(2550)   -- all 100%
  even : 3750/3750      odd  : 1920/1920
  defect (a) even char  : 3750/3750
  defect (b) nu > e     : 6/6
  repeated roots (P>1)  : 720/720
  64 distinct factor-degree profiles observed
```

Four independent checks run per instance:
`ok_structure` ($n=n'P$, $\gcd(n',p)=1$, $P=p^\nu$), `ok_mu` ($\mu^P=\lambda$,
the defect-(b) check), `ok_product` ($\prod f_i^{m_i}=x^n-\lambda$),
`ok_irreducible` (every factor irreducible, monic, multiplicity exactly $P$).

**Long-run discipline honoured:** pilot (200 inst, 0.4 s) → estimate committed
to `results/phase1_pilot_estimate.json` → checkpoints every 500 → full run.
Note the pilot's projected 11.5 s was **optimistically biased** (it sampled
only the smallest $q$); the real cost was 77.0 s. The estimate was committed
*before* the full run, as required.

**Reproducibility:** two independent full runs produce **byte-identical**
results (5 670 rows, ignoring per-instance timing). Verified by diff.

### C.3 CP1 gate — the four published factorisations — **4/4 PASS**

| ID | Source | Field | $n$ | $\lambda$ | Mode | Got | Result |
|---|---|---|---|---|---|---|---|
| CP1-1 | P1 Ex. 4 | $\mathbb F_9$ | 3 | 2 | **EXACT** | $[(x{+}1,3)]$ | ✅ |
| CP1-2 | P1 Ex. 6 | $\mathbb F_{27}$ | 15 | 2 | **EXACT** | $[(x{+}1,3),(x^4{+}2x^3{+}x^2{+}2x{+}1,3)]$ | ✅ |
| CP1-3 | P1 Ex. 5 | $\mathbb F_{25}$ | 7 | order 3 | SHAPE | $[1,3,3]$ | ✅ |
| CP1-4 | P3 Ex. 3.6 | $\mathbb F_{81}$ | 21 | order 10 | SHAPE | $[1,3,3]$ | ✅ |

CP1-1/CP1-2 are **string-exact** against the printed polynomials. CP1-3/CP1-4
are shape-exact, matching how A.1 records them.

### C.4 The two Blueprint defect cases — both exercised, both pass

**Defect (a) — even characteristic, absolute trace must include $i=0$.**
`test_defect_a_absolute_trace_includes_i_equals_zero` checks five things,
including: the fast transitive implementation equals the **literal**
$\sum_{i=0}^{ed-1}h^{2^i}$; the $i=1$-omitting version differs; $\operatorname{Tr}(h)$
satisfies $y^2=y$ mod $f$ (is genuinely $\mathbb F_2$-valued); and across 24
random $h$ the correct trace splits $f$ non-trivially while the buggy one is
provably **not** $\mathbb F_2$-valued. **3 750 instances** in the regression
sweep carry `defect_a_even_char`.

**Defect (b) — $\nu>e$: $\mu=\sigma^{e-(\nu\bmod e)}(\lambda)$.**
`test_defect_b_nu_gt_e` runs $(\mathbb F_3,e{=}1,n{=}9)$, $(\mathbb F_3,1,27)$,
$(\mathbb F_5,1,25)$, $(\mathbb F_3,2,27)$, $(\mathbb F_2,1,4)$, $(\mathbb F_2,1,8)$
— all $\nu>e$, all passing, with an internal assertion that
$\mu^{p^\nu}=\lambda$ so a wrong $\mu$ cannot slip through.
**6 instances** in the regression sweep carry `defect_b_nu_gt_e`, including the
Blueprint's named crash case $q=3,n=9$.

### C.5 Independent oracle checks — **90 passed, 2 skipped**

`tests/test_oracle.py` uses `galois` 0.4.11 (Hostetter, MIT) **only** as an
oracle. Since our field and galois's field use *different* modulus polynomials,
comparisons use isomorphism-invariant quantities:
- **O1** exact degree-profile agreement for $\lambda$ in the prime subfield;
- **O2** the profile **multiset over all of $\mathbb F_q^*$**, which is
  isomorphism-invariant even though individual $\lambda$ are not;
- **O3** irreducibility agreement.

Every oracle test **also** runs an *intrinsic* certificate — product identity
plus irreducibility of every factor, which is a complete certificate since
factorisation into irreducibles is unique. So an oracle limitation cannot mask
a real failure.

**The 2 skips are the oracle's fault, not ours:** `galois` raises
`RuntimeError: Failed to find a non-trivial factor after 1000 tries` on
$x^{21}+1$ over $\mathbb F_{16}$ and $\mathbb F_{256}$. Our factoriser succeeds
on both, and the results are certified intrinsically. Recorded rather than
hidden.

### C.6 Agreement with the [PROBE-DONE] A.3 reference

The 7 §5.5 instances were cross-checked against the cycle data recorded in
`BLUEPRINT.md` A.3 during blueprinting. **All 7 agree exactly** on $B$, on the
$(a,d)$ cycle shapes, and on the predicted code count:

| Instance | $B$ | Shapes | $(P{+}1)^{\sum a}$ |
|---|---|---|---|
| C1/X1 $\mathbb F_9,n{=}3,k{=}1$ | 1 | $[(1,1)]$ | 4 |
| C2/X2 $\mathbb F_{25},n{=}7,k{=}1$ | 3 | $[(1,1),(1,3),(1,3)]$ | 8 |
| C3/X3 $\mathbb F_{27},n{=}15,k{=}2$ | 2 | $[(1,1),(1,4)]$ | 16 |
| C4/X4 $\mathbb F_{81},n{=}21,k{=}2$ | 2 | $[(1,1),(2,3)]$ | 64 |
| C5/X5 $\mathbb F_{81},n{=}5,k{=}3$ | 2 | $[(1,1),(4,1)]$ | 32 |
| C6/X6 $\mathbb F_{81},n{=}17,k{=}3$ | 2 | $[(1,1),(4,4)]$ | 32 |
| C7/X7 $\mathbb F_{81},n{=}25,k{=}3$ | 3 | $[(1,1),(4,1),(4,5)]$ | 512 |

This is a genuine independent cross-check: a from-scratch re-implementation
landing on the same cycle data. It is **not** the P1 acceptance gate, and it is
**not** a paper result.

### C.7 A structural result verified by computation (Class C support)

`test_k_zero_gives_the_classical_reciprocal_hence_cycles_of_length_at_most_two`
verifies, computationally rather than by assertion, that for $k=0$ we have
$j=e$, $\sigma^e=\mathrm{id}$, hence $f^\#=f^*$ (the classical reciprocal) and
$(f^*)^*=f$ — **so for $k=0$ every cycle has length $\le 2$**, on every
instance tested. This is the structural reason PA-1 (whose maps $f\mapsto f^*$,
$f\mapsto f^{\dagger}$ are involutions) cannot reach cycle length $a\ge3$, and
it is what underpins Blueprint Class C1. Conversely,
`test_a_ge_three_cycles_are_reachable` confirms this implementation produces
cycles of length $\ge3$ (observed lengths 4 at $\mathbb F_{81}$, $k=3$).

---

## D. PROOF OF INDEPENDENCE FROM `/tmp/probe/`

Independence is **machine-checked**, not asserted. `tests/test_independence.py`
(6 tests, all passing):

- **I1** No source file in `core/`, `scripts/` or `tests/` (other than the
  proof file itself) mentions the prototype directory or imports any prototype
  module (`gf`, `ff`, `pa1_check`, `pa1_table2`, `scan_a3*`).
- **I2** Every `core.*` module resolves **inside the repository root**.
- **I3** In a **fresh subprocess**, importing all four core modules loads
  **zero** prototype modules; the prototype directory is not on `sys.path`.
- **I3b** The prototype's sha256s are recorded in
  `results/independence_record.json`, and **no core file is byte-identical to
  any prototype file**.
- **I4** `galois` is imported **nowhere** in `core/`.
- **I5** No core module was loaded from a stale bytecode cache.

Recorded prototype identity (for later diffing):

| File | sha256 | bytes |
|---|---|---|
| `ff.py` (abandoned) | `91604e34df6b745d52795db2ecc85ee1cd550e24298d015409fab984cb68fcfb` | 5 400 |
| **`gf.py`** | `674f7cbd550ec9618fc54f2e376b832b71fda53d250ad4e6f94363b447f13171` | 12 883 |
| `pa1_check.py` | `6d36cb7d8f901229e26a82c7af8ca8b9c48d6707f0c149442d87541c6f65fb05` | 1 433 |
| `pa1_table2.py` | `faaafd600a9ebc29434ba715dc082e251324d9b0c118dae1e9f1dece286075d2` | 2 838 |
| `scan_a3.py` (superseded) | `a75d6bd420b2c151f96a2ce00964490cb00ffbfdcc824b1e41e9ee45f06d5956` | 2 108 |
| `scan_a3b.py` | `7dda4c2594eaf90c61dcf85c8f7b0c3775559a67554c29076bddad2887283140` | 2 282 |
| `scan_a3P.py` | `a8e76e36f8b9a2658c665f41def1cc520ce085f7ca4c2576894056549fa86ce1` | 2 712 |

**`galois` 0.4.11 was used only as an oracle**, in `tests/test_oracle.py`
only, never inside `core/` (test I4).

---

## E. FILES CREATED / MODIFIED

**Created — implementation**

```
galois_hull_enum/__init__.py
galois_hull_enum/core/__init__.py
galois_hull_enum/core/field.py
galois_hull_enum/core/poly.py
galois_hull_enum/core/factor.py
galois_hull_enum/core/cycles.py
```

**Created — tests**

```
galois_hull_enum/tests/test_field.py
galois_hull_enum/tests/test_poly.py
galois_hull_enum/tests/test_factor.py
galois_hull_enum/tests/test_cycles.py
galois_hull_enum/tests/test_oracle.py
galois_hull_enum/tests/test_independence.py
```

**Created — scripts**

```
galois_hull_enum/scripts/gen_regression_suite.py
galois_hull_enum/scripts/run_phase1.py
```

**Created — instances and results**

```
galois_hull_enum/instances/regression_385.yaml     (OPEN-22 deliverable)
galois_hull_enum/results/regression_manifest.csv   (5 670 instances)
galois_hull_enum/results/phase1_factorisation.csv  (5 670 result rows)
galois_hull_enum/results/phase1_summary.json
galois_hull_enum/results/phase1_pilot_estimate.json
galois_hull_enum/results/independence_record.json
galois_hull_enum/results/phase1_tests.log         (verbose pytest log)
galois_hull_enum/results/phase1_regression.log    (full regression log)
```

**Modified**

```
requirements.txt      (+ pyyaml==6.0.3; note that galois is deliberately absent)
docs/P1_REPORT.md     (this file)
```

**Untouched:** `BLUEPRINT.md` (sha256 `2e02572a…`), `EXECUTION_PLAN.md`
(sha256 `6026536e…`), `Galios Hull.zip`, `sources/`, P0 report. Verified by
`git status` before commit.

---

## F. EXACT NUMERICAL RESULTS ACTUALLY OBTAINED

Every number below came from an executed run in this session. Nothing is
carried over from memory or from the probe.

**Category (b) — newly defined validation suite** (NOT the original 385 suite):

| Quantity | Value |
|---|---|
| Validation-suite instances run | **5 670** |
| Validation-suite instances passed | **5 670** |
| Validation-suite failures | **0** |
| Per-instance checks passed | 5 670 × 4 (structure, $\mu^P$, product, irreducibility) |
| Even-characteristic instances (defect a) | 3 750 / 3 750 |
| Odd-characteristic instances | 1 920 / 1 920 |
| $\nu>e$ instances (defect b) | 6 / 6 |
| Repeated-root instances ($P>1$) | 720 / 720 |
| Distinct factor-degree profiles seen | 64 |
| Validation-suite runtime | 77.0 s (73.7 inst/s) |
| Pilot | 200 inst / 0.4 s / 493.6 inst/s → projected 11.5 s (biased low; real 77.0 s) |
| Reproducibility | two runs byte-identical over 5 670 rows |
| **Original 385-instance suite reproduced** | **0 / 385 — OPEN / NOT RECONSTRUCTIBLE (OPEN-22)** |

**Category (c) — published-paper benchmark (CP1 quartet):**

| Quantity | Value |
|---|---|
| CP1 quartet reproduced | **4 / 4** (2 EXACT, 2 SHAPE — see OPEN-29) |

**Cycle detection on the 7 frozen §5.5 instances** (Blueprint-frozen $k$;
X1–X4 are the published-paper benchmarks, X5–X7 are Blueprint-frozen):

| Quantity | Value |
|---|---|
| Cycle shapes, 7 §5.5 instances | agree with A.3 on $B$, shapes, and $(P{+}1)^{\sum a}$ |
| Maximum cycle length observed | **4** ($\mathbb F_{81}$, $n=5,17,25$, $k=3$) |
| Fields constructed | $\mathbb F_{2,3,4,5,8,9,16,25,27,32,49,64,81,256}$ |
| Unit + oracle tests | **375 passed**, 5 skipped, 0 failed |

Cycle shapes for the 7 frozen instances are in §C.6. Full per-instance data for
the validation suite: `results/phase1_factorisation.csv`.

---

## G. OPEN ITEMS

### G.1 Closed in P1

| ID | Item | Disposition |
|---|---|---|
| — | Squarefree decomposition in small characteristic | ✅ Closed — new defect found by P1's own tests (Yun invalid for small $p$); fixed with Musser's algorithm. Not a Blueprint defect. |

### G.2 **STILL OPEN** — the original 385-instance suite

| ID | Item | Status | Impact | Owner |
|---|---|---|---|---|
| **OPEN-22** | The original 385-instance regression suite is **not reconstructible**: its exact $\lambda$ and $k$ assignments, and its odd-characteristic $q$ and $n$ assignments, are absent from the frozen Blueprint. Only the aggregate "385 instances, 0 failures" is recorded. | **OPEN / NOT RECONSTRUCTIBLE** | The original suite cannot be re-run or verified from committed inputs. Mitigated (not resolved) by `regression_suite_v1`, which is a *replacement*, not a recovery. | **P4/P9** — would need the original instance list; if it cannot be recovered, the paper must not cite the 385-instance figure as reproducible |

### G.3 Carried / newly raised

| ID | Item | Impact | Owner |
|---|---|---|---|
| **OPEN-28** | **Retained as the provenance label for newly chosen odd-$q$ inputs.** The 7 odd $q$-values in `regression_suite_v1` are a *choice of test inputs* by me, not a recovery of the inherited suite. They are labelled as such in the YAML and must stay labelled. | None on correctness; provenance only. | retained (documented) |
| **OPEN-29** | CP1-3 / CP1-4 use **shape** comparison only: the source papers' primitive elements $\omega\in\mathbb F_{25},\mathbb F_{81}$ are not recoverable from the PDFs (extracted text mangles the field-definition displays), and $\lambda$ is only determined up to Frobenius conjugacy. | Two of four CP1 checks are shape-level, not coefficient-level. Matches how A.1 records them ("identical shape"). | **P4** (re-read the source PDFs; XC-11 needs a second transcription of PA-1 Table 2 anyway) |
| **OPEN-30** | `galois` cannot factor $x^{21}+1$ over $\mathbb F_{16}$ / $\mathbb F_{256}$ (`RuntimeError` after 1000 tries). Our factoriser succeeds and is certified intrinsically. | 2 oracle tests skipped. | none — oracle limitation |
| **OPEN-23** | `instances/instances.yaml` (X1–X15 grid) not yet machine-readable | blocks P3/P4/P6/P8 | P3 |
| **OPEN-24** | `benchmarks/targets.yaml` (§7.1 targets) not yet built | blocks P4 | P4 |
| **OPEN-25** | `docs/SOURCE_INVENTORY.md` not yet built (needs XC-1) | P0 remnant | — |
| **OPEN-26** | `references.bib` not yet built (needs XC-2) | P9 | P9 |
| **OPEN-27** | `scripts/setup_env.sh` untested in a pristine sandbox | future phases | P2 |
| **OPEN-00/01/02/03/04/05/06/07/08/09/10/11/12/16/17/19/20** | Inherited, unchanged | per P0 report | per P0 report |

---

## H. P1 VERDICT

### H.1 Against the brief's acceptance criteria

| # | Criterion | Result |
|---|---|---|
| 1 | `regression_385.yaml` complete and provenance-backed, **or** every unresolved field explicitly OPEN with no invented values | ✅ **the second** — 5 fields OPEN with reasons, 0 invented values |
| 2 | Independent implementation exists under `core/` | ✅ 4 modules, ~1 250 lines |
| 3 | No dependency on `/tmp/probe/` | ✅ **machine-checked** (§D) |
| 4 | Core tests pass | ✅ 279 core tests (field+poly+factor+cycles) |
| 5 | Independent oracle checks pass | ✅ 90 passed, 2 skipped (OPEN-30, documented oracle limitation) |
| 6 | All required P1 identities verified by actual computation | ✅ 5 670 × 4 checks (validation suite) + CP1 4/4 + A.3 7/7 corroboration |
| 7 | Results reproducible from a clean invocation | ✅ two runs byte-identical |
| 8 | No claim marked PASS without executed evidence | ✅ every PASS traces to a logged run |

### H.2 Verdict

> ## **P1 Core Implementation: PASS**
> ## **New 5 670-instance validation suite: PASS**
> ## **Original 385-instance regression reproducibility: OPEN / NOT RECONSTRUCTIBLE**
>
> The core is implemented independently from scratch, is free of any dependency
> on the throwaway prototype (machine-checked), and passes **5 670 / 5 670**
> instances of the **newly defined** validation suite with **0 failures**,
> reproduces all **four** published CP1 factorisations, agrees with the A.3
> probe reference on all 7 cycle instances, and passes **375** unit and oracle
> tests.
>
> **OPEN-22 is NOT resolved and remains OPEN.** The original inherited
> 385-instance regression suite is not reconstructible from the frozen
> Blueprint — its exact $\lambda$ and $k$ assignments, and its
> odd-characteristic $q$ and $n$ assignments, are absent — so it is recorded
> as OPEN with per-field reasons. The 5 670 / 5 670 result is the **newly
> defined** validation suite and is **not** a recovery of the original 385
> suite. No instance count, $\lambda$, $k$, $q$, $n$ or expected result was
> invented, and none was inferred from the reported totals. The newly chosen
> odd-$q$ inputs carry the provenance label **OPEN-28**.
>
> **Three instance categories are kept strictly separate** (§B.5):
> (a) inherited/probe-derived — **NOT RECONSTRUCTIBLE**;
> (b) newly constructed validation — **PASS 5 670 / 5 670**;
> (c) published-paper benchmark — **PASS 4 / 4**. A result from one category
> is never reported as a result from another.
>
> **Retained OPEN items:** **OPEN-29** (CP1-3/CP1-4 are shape-only, because the
> sources' primitive elements in $\mathbb F_{25}$/$\mathbb F_{81}$ are not
> recoverable from the mangled PDFs) and **OPEN-30** (`galois` fails on
> $x^{21}+1$ over $\mathbb F_{16}$/$\mathbb F_{256}$; 2 oracle tests skipped).
>
> **New defect found and fixed in P1:** Yun's squarefree decomposition is
> invalid in small characteristic; replaced with Musser's algorithm.
>
> **Constraints honoured.** No P2+ computation was performed — no cycle
> polynomials $W_a^{(P)}$, no enumerator $N$, no hull dimensions, no moments.
> Both frozen documents are byte-identical to `076f910` / `96c7885`. No
> manuscript text was written. `galois` was used only as an oracle.

### H.3 Ready for P2

P2 (transfer matrix $T_P$, local weights $W_a^{(P)}$, Definitions 1–3, T3/T4)
can proceed. It needs `core/poly.py` and `core/field.py`, which are now
verified, plus the F2/F3 cross-checks of the source formulas.

**STOPPING after P1 as instructed.**
