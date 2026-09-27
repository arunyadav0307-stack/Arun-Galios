# P3 REPORT — GLOBAL ENUMERATOR VALIDATION (CP3)

**Phase:** P3 = **CP3**, the decisive gate of `BLUEPRINT.md` §6.2
**Frozen source:** `BLUEPRINT.md` @ `076f910` · `EXECUTION_PLAN.md` @ `96c7885`
**Builds on:** P1 @ `7d46765` · P2 @ `bde43fc` (tag `p2`)
**P3a/P3b checkpoint:** `0fa28f5` (estimate committed before the full run)
**P3c driver fixes:** `1d13d3e`
**Executed:** 2026-09-27/28
**Scope:** the global hull-dimension enumerator and its validation. **No P4+ work.**

> ## P3 VERDICT: **PASS**
>
> **OPEN-31 is CLOSED with actual evidence.** The global enumerator
> $N(y)=\prod_{\mathfrak c}W_{a(\mathfrak c)}^{(P)}(y^{d(\mathfrak c)})$ agrees,
> coefficient by coefficient, with two logically independent brute-force
> oracles on **4,377 / 4,377** instances, comparing **29,803 coefficients** and
> **14,446,825 codes**, with **zero mismatches**.
>
> **P1/P2 status distinctions preserved:** OPEN-22 original 385-instance suite =
> OPEN / NOT RECONSTRUCTIBLE; the 5,670-instance suite = validated new suite;
> OPEN-28 / OPEN-29 / OPEN-30 remain open. Frozen `BLUEPRINT.md` and
> `EXECUTION_PLAN.md` byte-identical throughout.

---

## A. WHAT WAS BUILT

| Module | Lines | Contents |
|---|---|---|
| `enumeration/enumerator.py` | 200 | **Algorithm E2** — $z\leftarrow y^{d}$ substitution, the per-cycle factor, the global product, the CP3 structural predicates |
| `oracle/brute.py` | 285 | **Algorithm E3**, twice, plus an independent $\#$ map and a factorisation certificate |
| `scripts/p3lib.py` | 210 | frozen §5.5 grid → instances; standing-hypothesis filters; X14 sweep construction |
| `scripts/run_phase3.py` | 630 | driver: `pilot` / `estimate` / `full` / `novelty`, checkpointed and resumable |
| `tests/test_phase3.py` | 260 | 39 tests for E2, both oracles, the $\#$ map, the certificate, and the stale-cache regression |

**Directory-name deviation (recorded, not silent).** BLUEPRINT §6.1 specifies
`enum/`. That name cannot be used in this repository: Python 3.11's `typing`
(imported by every core module) does `import enum` internally, and a top-level
package named `enum` on `sys.path` shadows the standard-library one, breaking
`import typing` with `AttributeError: module 'enum' has no attribute
'global_enum'`. The package is named `enumeration/`. Contents, names and
behaviour are exactly what the Blueprint specifies — a packaging constraint,
not a mathematical change.

---

## B. P3a — PILOT REPORT

### B.1 The three routes

| route | what it computes | shares with E2 |
|---|---|---|
| **E2** | $N(y)=\prod_{\mathfrak c}W_{a(\mathfrak c)}^{(P)}(y^{d(\mathfrak c)})$ | — |
| **oracle_cycle** | E3 *exactly as §5.3 writes it*: $u\in\prod_{\mathfrak c}[0,P]^{a(\mathfrak c)}$; build $g=\prod f_m^{u_m}$, $h=(x^n-\lambda)/g$, $h^{\#}$, $\dim=n-\deg\operatorname{lcm}(g,h^{\#})$; histogram | the cycle partition, used as an **indexing device only** |
| **oracle_flat** | the same mathematics indexed by a **flat** exponent vector over the distinct irreducible factors of $x^n-\lambda$ | **no cycle code at all** |

The frozen E3 oracle's dimension formula never mentions a cycle: it builds real
polynomials and takes a gcd. `oracle_flat` removes even the indexing dependency,
so if all three agree, the agreement cannot be an artefact of the cycle
partition — which is the direct answer to the brief's "hidden implementation
dependency" requirement.

### B.2 Oracle-independence apparatus

- **`hash_map_independent`** — a fresh implementation of
  $f^{\#}(x)=\sum_i f_0^{-p^{j}}f_i^{p^{j}}x^{m-i}$ written straight from §3.1,
  applying $\sigma^{j}$ one coefficient at a time. It does **not** import
  `core.cycles.hash_map`, so it is a separate code path; the test suite requires
  the two to agree on 10 $(p,e,n,\lambda,k)$ combinations, and requires
  $(f^{\#})^{\#}=f$ at $k=0$ (where $\sigma^{e}=\mathrm{id}$).
- **`certify_factorisation`** — re-derives from scratch, on every instance, that
  $\prod_i f_i^{P}=x^n-\lambda$, each $f_i$ is monic and irreducible with
  $f_i(0)\ne0$, the factors are pairwise distinct, and all multiplicities equal
  $P$. A test confirms it *rejects* a deliberately perturbed factorisation.

### B.3 Shared-dependency analysis

| shared component | E2 | oracles | verdict |
|---|---|---|---|
| `core/field.py` | no | yes | CP1-verified kernel; encodes no hull theory |
| `core/poly.py` | no | yes | same |
| `core/zz.py` | yes | no | same |
| `core/cycle_poly.py` | **yes** | no | P2-verified; the object under test |
| `factor_xn_minus_lambda` | yes (via E1) | yes | **the one genuine shared dependency** |
| `core.cycles.hash_map` | yes (via E1) | **no** | removed |
| cycle construction | yes (via E1) | `oracle_cycle` only | **removed from `oracle_flat`** |
| transfer matrix $T_P$ | **yes** | no | never appears in either oracle |

The single remaining shared dependency is `factor_xn_minus_lambda`. It is not
part of the enumerator — it is the CP1 field/factorisation kernel, verified in
P1 against the factorisations printed in the sources and **re-certified inside
the oracle on every instance**. No hull-dimension logic passes through it.

### B.4 Pilot (10 instances, all frozen §5.5 rows — nothing invented)

10/10 verified. All regimes the brief requires: $a=1$, $a=2$, genuine $a\ge3$
(4-cycles), $P=1$, $P>1$ ($P=3,4$), odd and even characteristic.
47 coefficients and 33,593 codes compared, 0 mismatches, peak RSS 14 MB.
Full detail in `docs/P3_PILOT_AND_ESTIMATE.md`.

---

## C. P3b — COST ESTIMATE (committed before the full run)

Measured from the pilot (pessimistic per-code rate = max over pilot instances):

| component | estimated (s) | **actual (s)** | ratio |
|---|---|---|---|
| `oracle_flat` | 3,767.5 | **5,343.6** | 1.42× |
| `oracle_cycle` | 2,409.5 | **5,670.7** | 2.35× |
| E2 | 1.6 | **2.2** | 1.34× |
| sweep construction | 645.3 | 571.7 + 645.3 | — |
| **total** | **6,823.8** | **11,016.5** | **1.61×** |

The estimate was **optimistic by 1.6×**, and honestly so: the pessimistic
per-code rate was taken from the pilot's largest instance ($\mathbb F_{16}$,
$n=15$), while the sweep reaches $n=60$ and $\lvert\mathscr C\rvert=10^{5}$, where
per-code cost is higher. The estimate was committed *before* the run and is not
retrospectively adjusted. Verdict was PROCEED, and the run completed in 3 h 04 m.

**Checkpoint boundaries (defined before the run):** CP3-1 … CP3-5, formed by
descending instance cost. In the event the run was checkpointed per instance to
`results/phase3_sweep_checkpoint.jsonl` (flush + `fsync`), and the sandbox
terminated the process at 2,643/4,377 — **all 2,643 completed instances were
recovered by `--resume` and the run finished with zero recomputation of verified
work.** That is the checkpoint mechanism doing its job.

**Complete-run definition (fixed before the run):** every instance from the
frozen X14 filter is run through all three routes; all three agree on every
coefficient; $\lvert\mathscr C\rvert=(P+1)^{\sum a}$ equals the flat oracle's
code count; $\#\{\dim=0\}=2^B$; all coefficients nonnegative integers.
Instances with $\lvert\mathscr C\rvert>10^{5}$ are **not** verified — recorded as
enumerator-only per the frozen validation-scope amendment.

---

## D. P3c — FULL VALIDATION REPORT

### D.1 Headline result

```
instances run       : 4,377
verified            : 4,377        (E2 == oracle_cycle == oracle_flat, coefficient by coefficient)
failed              : 0
coefficients cmp    : 29,803
codes compared      : 14,446,825
mismatched coeffs   : 0
max |coeff diff|    : 0
structural failures : []
seconds             : e2 2.2 | oracle_cycle 5,670.7 | oracle_flat 5,343.6 | total 11,016.5
peak RSS            : 65.7 MB
```

**This is not a spot check.** Every instance in the sweep was run in full; the
`n_verified` count equals the sweep size.

### D.2 The frozen §5.5 grid

| row | $q$ | $n$ | $k$ | $r$ | $P$ | shapes | $\lvert\mathscr C\rvert$ | result |
|---|---|---|---|---|---|---|---|---|
| X1 | 9 | 3 | 1 | 2 | 3 | $1{\times}1$ | 4 | verified |
| X2 | 25 | 7 | 1 | 3 | 1 | $1{\times}1,1{\times}3,1{\times}3$ | 8 | verified |
| X3 | 27 | 15 | 2 | 2 | 3 | $1{\times}1,1{\times}4$ | 16 | verified |
| X4 | 81 | 21 | 2 | 10 | 3 | $1{\times}1,\mathbf{2{\times}3}$ | 64 | verified |
| X5 | 81 | 5 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}1}$ | 32 | verified |
| X6 | 81 | 17 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}4}$ | 32 | verified |
| X7 | 81 | 25 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}1},\mathbf{4{\times}5}$ | 512 | verified |
| X8a | 16 | 5 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}1}$ | 32 | verified |
| X8b | 16 | 15 | 3 | 1 | 1 | $1{\times}1{\times}3,\mathbf{4{\times}1}{\times}3$ | 32,768 | verified |
| X8c | 16 | 17 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}2}{\times}2$ | 512 | verified |
| X9 | 64 | 20 | 3 | 9 | 4 | $1{\times}1,\mathbf{2{\times}2}$ | 125 | verified |
| X10a | 25 | 65 | 1 | 6 | 5 | $1{\times}1,2{\times}2{\times}3$ | 279,936 | **ENUMERATOR ONLY** |
| X10b | 25 | 65 | 1 | 3 | 5 | $1{\times}1,2{\times}2{\times}3$ | 279,936 | **ENUMERATOR ONLY** |

Bold marks the cycles that matter for the novelty sub-suite. **X10a/X10b exceed
the frozen E3 bound $\lvert\mathscr C\rvert\le10^{5}$ and are flagged, not
verified** — exactly as the Blueprint's validation-scope amendment requires.
They are *not* counted as passes anywhere in this report.

### D.3 The X14 sweep — coverage

| dimension | range |
|---|---|
| fields $q$ | all 14 frozen fields: 4, 8, 9, 16, 25, 27, 32, 49, 64, 81, 121, 125, 128, 256 |
| $n$ | 1 – 60 |
| $k$ | 0 – 7 |
| characteristic $p$ | 2, 3, 5, 7, 11 |
| $P=p^{\nu}$ | 1, 2, 3, 4, 5, 7, 8, 9, 11, 16, 25, 27, 32, 49 |
| cycle lengths $a$ | 1 – 7 |
| $r=\operatorname{ord}(\lambda)$ | 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 17 |
| $B$ | 1 – 16 |
| max hull dimension seen | 30 |
| verified per field | F_4:158 · F_8:174 · F_9:200 · F_16:351 · F_25:217 · F_27:255 · F_32:285 · F_49:214 · F_64:452 · F_81:440 · F_121:248 · F_125:249 · F_128:406 · F_256:728 |

**376 candidates rejected**, every one written to
`results/phase3_full_rejected.csv` with its reason (hypothesis failure, or
$\lvert\mathscr C\rvert>10^{5}$). Coverage is auditable, not implied.

**$\lambda$ sampling (recorded).** $\lambda$ is sampled **by multiplicative
order** — one representative per order — rather than exhaustively over
$\mathbb F_q^{*}$, which for $q=256$ would be $255\times60\times e$ candidates.
The cycle structure depends on $\lambda$ through $r=\operatorname{ord}(\lambda)$
and the concrete factorisation, and every candidate is factorised and certified
from scratch, so nothing about the factorisation is assumed. This is a
documented restriction of the sweep's breadth, not a hidden one.

### D.4 Structural checks — all 4,377 instances

| predicate | result |
|---|---|
| $\lvert\mathscr C\rvert=(P+1)^{\sum a}$ = the flat oracle's own code count | 4,377 / 4,377 |
| $\#\{\dim=0\}=2^{B}$ (T5 / LCD count) | 4,377 / 4,377 |
| all coefficients nonnegative integers | 4,377 / 4,377 |
| $\sum_{\mathfrak c}a(\mathfrak c)$ = number of distinct irreducible factors | 4,377 / 4,377 |

### D.5 Largest verified instances

| $q$ | $n$ | $k$ | $P$ | $B$ | shapes | $\lvert\mathscr C\rvert$ | coeffs | LCD | $2^{B}$ |
|---|---|---|---|---|---|---|---|---|---|
| 81 | 45 | 1 | 9 | 2 | $1{\times}1,\mathbf{4{\times}1}$ | 100,000 | 23 | 4 | 4 |
| 81 | 45 | 2 | 9 | 5 | $1{\times}1{\times}5$ | 100,000 | 21 | 32 | 32 |
| 81 | 45 | 3 | 9 | 2 | $1{\times}1,\mathbf{4{\times}1}$ | 100,000 | 23 | 4 | 4 |
| 81 | 45 | 0 | 9 | 3 | $1{\times}1,2{\times}1{\times}2$ | 100,000 | 23 | 8 | 8 |

These sit exactly on the frozen feasibility bound and were certified in full.

---

## E. P3d — NOVELTY SUB-SUITE REPORT

The region PA-1 cannot reach, selected because it is **the frozen Blueprint's
own validation grid**, not because it is easy.

| regime | instances | verified | codes | max cycle length | $P$ / $r$ values |
|---|---|---|---|---|---|
| genuine $a\ge3$ | **530** | **530** | 5,772,336 | **7** | — |
| repeated-root $P>1$ | **1,731** | **1,731** | 7,177,871 | — | $P\in\{2,3,4,5,7,8,9,11,16,25,27,32,49\}$ |
| **$a\ge3$ AND $P>1$** | **192** | **192** | 3,240,144 | 6 | $P=2$ |
| general $\lambda$, $r\notin\{1,2\}$ | **984** | **984** | 3,238,429 | — | $r\in\{3,4,5,6,8,9,10,12,17\}$ |

Every instance in every regime passes. The hardest cell — **192 instances with
simultaneously $a\ge3$ and $P>1$** (repeated roots *and* cyclically coupled
exponents, i.e. exactly where PA-1's product-of-independent-one-parameter-sums
form cannot express the count) — is fully verified, with 3,240,144 codes
compared and zero mismatches.

Representative $a\ge3 \wedge P>1$ instances:

| $q$ | $n$ | $k$ | $P$ | shapes | $\lvert\mathscr C\rvert$ |
|---|---|---|---|---|---|
| 8 | 14 | 1 | 2 | $1{\times}1,\mathbf{6{\times}1}$ | 2,187 |
| 8 | 18 | 1 | 2 | $1{\times}1,1{\times}2,\mathbf{3{\times}2}$ | 243 |
| 8 | 26 | 1 | 2 | $1{\times}1,\mathbf{3{\times}4}$ | 81 |
| 8 | 36 | 1 | 4 | $1{\times}1,1{\times}2,\mathbf{3{\times}2}$ | 3,125 |

**Verdict: P3d PASS.** The $a\ge3$ suite is traceable to the frozen grid (X5–X8
and the X14 sweep), not selected for convenience.

---

## F. COEFFICIENT-LEVEL COMPARISON STATISTICS

| statistic | value |
|---|---|
| instances with **complete** distribution agreement | 4,377 / 4,377 |
| coefficients compared (distinct dimensions, over all instances) | 29,803 |
| codes compared | 14,446,825 |
| mismatched coefficients | **0** |
| maximum absolute coefficient difference | **0** |
| coefficients per instance (min / median / max) | 1 / 4 / 31 |
| instances where E2 = `oracle_cycle` only | 0 |
| instances where E2 = `oracle_flat` only | 0 |
| structural failures | 0 |

Machine-readable: `results/phase3_full_coefficients.csv` (one row per
instance-dimension, with the three routes' coefficients side by side and an
explicit agreement flag).

### F.1 Independent re-verification

To rule out any artefact of the checkpointed aggregation, a **fresh process**
re-ran a stratified sample of 60 instances (12 each from: $a\ge3$; $a\ge3\wedge
P>1$; general $\lambda$; $\lvert\mathscr C\rvert\ge20{,}000$; uniformly random)
from scratch, recomputing E2 and **both** oracles and comparing against the
stored records. **60 / 60 agreed**, including the code counts.

---

## G. FAILURES AND THEIR RESOLUTIONS

Six defects were found by *running* the code. **None was a formula error; none
was repaired by assumption.** Each is recorded in full rather than quietly
fixed.

| # | where | symptom | diagnosis | resolution |
|---|---|---|---|---|
| 1 | `oracle_flat` | wrong histograms on exactly the $a\ge3$ instances ($\mathbb F_{81}$, $n=5$ and $n=25$, $k=3$) | the incremental `suffix[i]=∏_{l≥i}f_l^{e_l}` cache went **stale on an odometer carry**: a wrapped digit was reset to $0$ without refreshing its suffix entry. The two oracles' enumerated **code sets were verified identical** (32=32, 512=512, no code in one and absent from the other), so only the cache was wrong | **cache removed**, not patched; `g` is rebuilt from precomputed powers each step |
| 2 | `cmd_full` | would have pushed X10a/X10b (279,936 codes each) through the oracle | over-bound grid rows were not separated before the oracle loop, breaching the frozen $\lvert\mathscr C\rvert\le10^{5}$ feasibility bound | over-bound rows routed to the **flagged enumerator-only** list |
| 3 | `--resume` | `KeyError: 'key'` on the first checkpoint line | `run_instance` did not put `"key"` in the record it returns, though the resume set was built from `inst.key` | `"key"` added to the record; fold-in made robust with `.get()` |
| 4 | `aggregate()` | `KeyError: 'max_abs_coefficient_difference'` | the comparison dict uses `max_abs_coefficient_diff` | key name corrected |
| 5 | `enumeration/enumerator.py` | `AssertionError: trace != brute force` for $(P,a)=(3,1)$ | `cycle_poly_bruteforce` returns a **dict** `{weight: count}`; I passed it straight to `zz.ztrim`, which iterated its *keys* | routed through `dict_to_poly`; the cross-check then passed and is kept as free insurance |
| 6 | `tests/test_phase3.py` | 5 test failures | all five were **my test errors**: two out-of-range `k` values; one key belonging to the *driver's* `structural()` rather than the module's `structural_checks()`; one wrong $\lambda$ order (`from_int(2)` has order 4 in $\mathbb F_{25}$ and fails (H), where X10 needs order 6); one wrong claim that $g=1$ gives hull dimension $n$ — in fact $g=1$ is the whole space, whose $k$-Galois dual is $\{0\}$, so the hull is $\{0\}$ and the dimension is **0** | all five corrected with the arithmetic written out |

**No frozen formula was altered.** The Blueprint and execution plan are
byte-identical to `076f910` / `96c7885` before and after.

### G.1 Environment changes (recorded, not mathematical)

The sandbox lost two packages mid-session: `pytest` and `galois`. Both were
reinstalled (`pytest 9.1.1`, `galois 0.4.11`, plus `numpy 2.4.6`). The full
suite returns to **605 passed / 5 skipped**, matching P2's 566 + 39 new P3
tests, with the *same* 5 skips as P2 (3 "k out of range", 2 `galois`
$x^{21}+1$). No P3 number depends on `galois`; it is used only by P1's
independent-oracle tests.

---

## H. UPDATED OPEN LIST

### H.1 Closed by P3

| ID | item | disposition |
|---|---|---|
| **OPEN-31** | *The global enumerator $N$ (T2) has not been assembled or checked against a brute-force oracle over all codes.* | **CLOSED with evidence.** Assembled as Algorithm E2 and validated on 4,377 instances against two independent oracles: 29,803 coefficients, 14,446,825 codes, zero mismatches; 530 $a\ge3$ and 192 $a\ge3\wedge P>1$ instances included. |

### H.2 Still open (carried forward, unchanged)

| ID | item | status |
|---|---|---|
| **OPEN-22** | Original 385-instance suite not reconstructible (exact $\lambda$, $k$, odd $q$, $n$ absent from the frozen Blueprint) | **OPEN / NOT RECONSTRUCTIBLE** |
| **OPEN-28** | Provenance label for the newly chosen odd-$q$ inputs in `regression_suite_v1` | **retained** |
| **OPEN-29** | CP1-3/CP1-4 shape-only (sources' primitive elements in $\mathbb F_{25}$/$\mathbb F_{81}$ not recoverable from the mangled PDFs) | **OPEN** (P4) |
| **OPEN-30** | `galois` fails **flakily** on $x^{21}+1$ over $\mathbb F_{16}$/$\mathbb F_{256}$; 1–2 oracle tests skipped depending on the run (its EDF is randomised) | **OPEN** — oracle limitation, ours not implicated |
| **OPEN-23/24/25/26/27** | `instances.yaml`, `targets.yaml`, `SOURCE_INVENTORY.md`, `references.bib`, untested `setup_env.sh` | carried |
| **OPEN-32/33** | §3.5's optional structural reductions (Hadamard/rank-one form of $T_P$; LDL$^\mathsf T$ with $N[i,j]=(P+1-i-j)^+$) unverified | carried, optional |
| **OPEN-34** | $T_6$ (self-orthogonal / dual-containing counts) not implemented | **P7** |
| **OPEN-35** | Restricted/unrestricted mean difference not recomputed | **P6** |
| **OPEN-36** | $T_8$ (CLT) untouched | **P8** |

### H.3 New in P3

| ID | item | impact |
|---|---|---|
| **OPEN-37** | **X10a/X10b are not oracle-verified.** $\lvert\mathscr C\rvert=279{,}936>10^{5}$, so E3 cannot certify them. E2's distribution is computed and the structural checks pass, but there is **no independent confirmation**. | Flagged ENUMERATOR ONLY in every artefact; must not be cited as validated |
| **OPEN-38** | **Ring instances X11–X13 not run.** They need the idempotent decomposition of $A$ / $R_{m,q}$, which the frozen Blueprint assigns to CP6. | Deferred to CP6; recorded as not-attempted, not passed |
| **OPEN-39** | **$\lambda$ sampled by order, not exhaustively.** One representative per multiplicative order. | Breadth restriction of the X14 sweep; documented, coverage auditable via the rejected-CSV |

---

## I. ACCEPTANCE CRITERIA

| # | criterion | result |
|---|---|---|
| 1 | Complete-distribution agreement between independent methods | ✅ 4,377/4,377, 29,803 coefficients, 14,446,825 codes, 0 mismatches; E2 = `oracle_cycle` = `oracle_flat` |
| 2 | Global enumerator passes all required structural checks | ✅ $\lvert\mathscr C\rvert=(P+1)^{\sum a}$; $\#\{\dim=0\}=2^{B}$; nonnegative integer coefficients; $\sum a$ = factor count — all 4,377/4,377 |
| 3 | Genuine $a\ge3$ cases pass | ✅ 530 instances, max cycle length 7, 5,772,336 codes |
| 4 | Required $P>1$ cases pass | ✅ 1,731 instances, $P$ up to 49; 192 with $a\ge3$ **and** $P>1$ |
| 5 | Full run checkpointed and reproducible | ✅ per-instance checkpoint with `fsync`; `--resume` recovered all 2,643 instances after the sandbox killed the run; 60/60 stratified re-verification in a fresh process |
| 6 | Pilot → estimate → full-run evidence preserved | ✅ `docs/P3_PILOT_AND_ESTIMATE.md` + `results/phase3_pilot.json` + `results/phase3_estimate.json`, committed at `0fa28f5` **before** the full run |
| 7 | No unresolved mathematical mismatch hidden | ✅ six defects found by running the code, all diagnosed and documented in §G; none a formula error |
| 8 | OPEN-31 closed with actual evidence, or remains explicitly OPEN | ✅ **CLOSED** with the evidence above |

---

## J. REPRODUCIBILITY

```
python3 scripts/run_phase3.py pilot      # P3a
python3 scripts/run_phase3.py estimate   # P3b
python3 scripts/run_phase3.py full       # P3c  (~3 h; checkpointed)
python3 scripts/run_phase3.py novelty    # P3d
```

Exact integer arithmetic throughout; no floating point in the mathematical path
(`mean_exact` is emitted as a fraction string, never a float). The only
randomness anywhere is inside `galois`'s EDF, used solely by P1's independent
oracle tests.

**Artifacts:** `results/phase3_summary.json`, `phase3_full_instances.csv`,
`phase3_full_coefficients.csv`, `phase3_full_rejected.csv`,
`phase3_sweep_checkpoint.jsonl`, `phase3_pilot*.{json,csv}`,
`phase3_estimate.json`, `phase3_novelty.json`, `phase3_full_run.log`,
`phase3_tests.log`.

**Frozen state verified before commit:** `BLUEPRINT.md` sha256
`2e02572a0ef970a93dcf242c725f290938d9b180787beaccd8acedf6a609214f`,
`EXECUTION_PLAN.md` sha256
`6026536e3df10ac80a35225ddc42fc3b5f52f845be92da98531aca4bc8f65b32`.

---

## K. VERDICT

> ## **P3 = PASS**
>
> The global hull-dimension enumerator, assembled from P2's verified
> cycle-local polynomials, is validated against two logically independent
> brute-force oracles on 4,377 instances spanning all 14 frozen fields,
> $n\le60$, $k\le7$, $P\le49$, cycle lengths 1–7 and $r\le17$. **29,803
> coefficients and 14,446,825 codes compared; zero mismatches.**
>
> The $a\ge3$ region inaccessible to PA-1 is covered by 530 instances
> (192 of them with $P>1$ simultaneously), all passing.
>
> **OPEN-31 is closed with actual evidence.** OPEN-37 (X10a/X10b over the
> frozen E3 bound), OPEN-38 (ring instances) and OPEN-39 ($\lambda$ sampled by
> order) are recorded as explicit limitations rather than presented as passes.
>
> Six defects were found by executing the code and are documented in §G; none
> was a formula error, and no frozen formula was altered. P1 and P2 results are
> unchanged: 605 tests pass, 5 skipped.
>
> **Scope honoured:** no P4+ work, no moments, no limit law, no manuscript.
>
> **STOPPING after P3 as instructed.**
