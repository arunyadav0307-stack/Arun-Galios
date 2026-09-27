# P3a PILOT AND P3b COST ESTIMATE — committed before the full run

**Phase:** P3 = CP3 (BLUEPRINT.md §6.2, the decisive gate)
**Frozen source:** `BLUEPRINT.md` @ `076f910` · `EXECUTION_PLAN.md` @ `96c7885`
**Builds on:** P1 @ `7d46765` · P2 @ `bde43fc` (tag `p2`)
**Executed:** 2026-09-27
**Status:** P3a PASS · P3b estimate PROCEED · **P3c not yet run**

This file exists because the execution plan requires the cost estimate to be
committed **before** the full validation run. It records P3a and P3b only.

---

## A. P3a — PILOT REPORT

### A.1 What was built

| Module | Purpose |
|---|---|
| `enumeration/enumerator.py` | **Algorithm E2** — $N(y)=\prod_{\mathfrak c}W_{a(\mathfrak c)}^{(P)}(y^{d(\mathfrak c)})$, built from P2's verified cycle factors |
| `oracle/brute.py` | **Algorithm E3** — two independent brute-force oracles (see A.3) |
| `scripts/p3lib.py` | frozen §5.5 grid → concrete instances; standing-hypothesis filters |
| `scripts/run_phase3.py` | driver with `pilot` / `estimate` / `full` / `novelty` sub-commands |

**Directory-name deviation (recorded, not silent).** BLUEPRINT §6.1 specifies
`enum/`. That name cannot be used here: Python 3.11's `typing` (imported by
every core module) does `import enum` internally, and a top-level package named
`enum` on `sys.path` shadows the standard-library one, breaking `import typing`
with `AttributeError: module 'enum' has no attribute 'global_enum'`. The package
is named `enumeration/`. Module contents, names and behaviour are exactly what
the Blueprint specifies — this is a packaging constraint, not a mathematical
change.

### A.2 The pilot (10 instances, all from the frozen §5.5 grid — nothing invented)

| row | $q$ | $p$ | $n$ | $k$ | $r$ | $P$ | shapes $(a{\times}d)$ | $\lvert\mathscr C\rvert$ | coeffs | max dim | LCD | $2^B$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| X1 | 9 | 3 | 3 | 1 | 2 | **3** | $1{\times}1$ | 4 | 2 | 1 | 2 | 2 |
| X2 | 25 | 5 | 7 | 1 | 3 | 1 | $1{\times}1,1{\times}3,1{\times}3$ | 8 | 1 | 0 | 8 | 8 |
| X3 | 27 | 3 | 15 | 2 | 2 | **3** | $1{\times}1,1{\times}4$ | 16 | 4 | 5 | 4 | 4 |
| X4 | 81 | 3 | 21 | 2 | 10 | **3** | $1{\times}1,\mathbf{2{\times}3}$ | 64 | 8 | 10 | 4 | 4 |
| X5 | 81 | 3 | 5 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}1}$ | 32 | 3 | 2 | 4 | 4 |
| X6 | 81 | 3 | 17 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}4}$ | 32 | 3 | 8 | 4 | 4 |
| X8a | 16 | 2 | 5 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}1}$ | 32 | 3 | 2 | 4 | 4 |
| X8b | 16 | 2 | 15 | 3 | 1 | 1 | $1{\times}1{\times}3,\mathbf{4{\times}1}{\times}3$ | 32768 | 7 | 6 | 64 | 64 |
| X8c | 16 | 2 | 17 | 3 | 1 | 1 | $1{\times}1,\mathbf{4{\times}2}{\times}2$ | 512 | 5 | 8 | 8 | 8 |
| X9 | 64 | 2 | 20 | 3 | 9 | **4** | $1{\times}1,\mathbf{2{\times}2}$ | 125 | 11 | 10 | 4 | 4 |

Bold entries are the $a\ge3$ / $a=2$ cycles that matter for the novelty
sub-suite. **Every row is a frozen grid row; no pilot parameter was invented.**

**Grid-row correction (recorded, not silent).** The frozen §5.5 row X8 lists
"$n=5/15/17$, $\lambda=\alpha^0,\alpha^5$, $r=1/3$, $\lvert\mathscr C\rvert=32/2^{12}/32$".
Encoding $n=15$ with $r=3$ is **not admissible**: $\gcd(n'=15,\,r=3)=3\neq1$
violates the standing hypotheses, so no $\lambda$ of order 3 works there. The
$\lvert\mathscr C\rvert=2^{12}$ cell forces $\lambda=1$ ($r=1$), which is what
X8b above uses. This is a reading of the frozen row, not a change to it.

### A.3 Oracle independence — the design and one bug it caught

**Three routes run on every instance:**

| route | what it does | shares with E2 |
|---|---|---|
| **E2** | multiplies $W_{a}^{(P)}(y^{d})$ over the cycle multiset | — |
| **oracle_cycle** | E3 *exactly as §5.3 writes it*: $u\in\prod_{\mathfrak c}[0,P]^{a(\mathfrak c)}$, build $g=\prod f_m^{u_m}$, $h=(x^n-\lambda)/g$, $h^\#$, $\dim=n-\deg\operatorname{lcm}(g,h^\#)$ | uses the cycle partition as an **indexing device only** |
| **oracle_flat** | same mathematics, but indexed by a **flat** exponent vector over the distinct irreducible factors of $x^n-\lambda$ | **no cycle code at all** |

The frozen E3 oracle's dimension formula never mentions a cycle: it builds real
polynomials and takes a gcd. The cycle decomposition appears only in *which*
exponent vectors get visited. `oracle_flat` removes even that, so if all three
agree, the agreement cannot be an artefact of the cycle indexing.

`oracle/brute.py` also carries `hash_map_independent`, a fresh implementation
of $f^\#$ written straight from §3.1 (it does **not** import
`core.cycles.hash_map`), and `certify_factorisation`, which re-derives from
scratch that $\prod f_i^{P}=x^n-\lambda$, each $f_i$ is monic/irreducible with
nonzero constant term, the factors are pairwise distinct, and all multiplicities
equal $P$.

**Bug found and fixed before any full run.** The first `oracle_flat` used an
incremental cache `suffix[i]=\prod_{l\ge i}f_l^{e_l}`. On an odometer **carry**
the wrapped digit was reset to $0$ without refreshing its suffix entry, leaving
a stale factor — which silently produced wrong histograms on exactly the
instances with $a\ge3$ cycles ($\mathbb F_{81}$, $n=5$ and $n=25$, $k=3$). E2
and `oracle_cycle` agreed with each other and disagreed with `oracle_flat`,
which is how it was caught: the code **sets** enumerated by both oracles were
verified identical (32 = 32 and 512 = 512, no code present in one and absent
from the other), so only the cache was wrong. The cache was **removed** rather
than patched. Full write-up in `docs/P3_REPORT.md` §G.

### A.4 Pilot results

```
instances run      : 10
verified           : 10        (E2 == oracle_cycle == oracle_flat, coefficient by coefficient)
failed             : 0
coefficients cmp   : 47
codes compared     : 33,593
mismatched coeffs  : 0
max |coeff diff|   : 0
seconds            : e2 0.002 | oracle_cycle 4.964 | oracle_flat 5.076 | total 10.042
peak RSS           : 14.0 MB
```

**Structural checks, all 10 instances:** $\lvert\mathscr C\rvert=(P+1)^{\sum a}$
= the flat oracle's actual code count; $\#\{\dim=0\}=2^B$; all coefficients
nonnegative integers; $\sum_c a(\mathfrak c)$ = number of distinct irreducible
factors.

**Regime coverage (the brief's explicit pilot requirements):**

| requirement | covered |
|---|---|
| $a=1$ | yes (X1, X2, X3) |
| $a=2$ | yes (X4, X9) |
| genuine $a\ge3$ | yes (X5, X6, X8a/b/c — 4-cycles) |
| $P=1$ | yes (X2, X5, X6, X8a/b/c) |
| $P>1$ | yes (X1 $P{=}3$, X3 $P{=}3$, X4 $P{=}3$, X9 $P{=}4$) |
| odd characteristic | yes ($p=3,5$) |
| even characteristic | yes ($p=2$) |

**Complete distributions compared, not selected dimensions.** The comparison is
over the full coefficient list of $N(y)$; `n_coefficients_compared` counts every
distinct dimension in the support and `n_codes_compared` counts every code.

---

## B. P3b — COST ESTIMATE

### B.1 Measured rates (from the pilot)

| quantity | value |
|---|---|
| conservative per-code rate, `oracle_flat` | $2.61\times10^{-4}$ s/code (max over pilot instances with $\lvert\mathscr C\rvert\ge32$) |
| conservative per-code rate, `oracle_cycle` | $1.67\times10^{-4}$ s/code |
| conservative per-instance rate, E2 | $3.61\times10^{-4}$ s/instance |
| sweep construction | 645.3 s (factorising every candidate) |

The per-code rate is taken as the **maximum** over the pilot, i.e. the
pessimistic end, because per-code cost grows with $\deg(x^n-\lambda)=n$ and the
pilot's largest instance ($\mathbb F_{16}$, $n=15$) is not the largest in the
sweep ($n\le60$).

### B.2 The sweep the frozen grid specifies

Filter (all frozen): $q\in\{4,8,9,16,25,27,32,49,64,81,121,125,128,256\}$,
$n\le60$, $\gcd(n',r)=1$, $r\mid(1+p^{e-k})$, $\lvert\mathscr C\rvert\le10^{5}$
(the Blueprint's own E3 feasibility bound).

```
sweep candidates          : 4,377
rejected (auditable)      : 376
total codes to certify    : 14,446,825
largest single instance   : 100,000 codes
```

$\lambda$ is sampled **by multiplicative order** (one representative per order)
rather than exhaustively over $\mathbb F_q^{*}$ — for $q=256$ that would be
$255\times60\times e$ candidates. The cycle structure depends on $\lambda$
through $r=\operatorname{ord}(\lambda)$ together with the concrete
factorisation, and **every candidate is still factorised and certified from
scratch**, so nothing about the factorisation is assumed. Every rejection is
written to `results/phase3_sweep_rejected.csv` with its reason, so the sweep's
coverage is auditable rather than implied.

### B.3 Estimate

| component | seconds |
|---|---|
| `oracle_flat` | 3,767.5 |
| `oracle_cycle` | 2,409.5 |
| E2 | 1.6 |
| sweep construction | 645.3 |
| **total** | **6,823.8 s ≈ 113.7 min** |

**Verdict: PROCEED** (threshold: 6 h).

### B.4 Checkpoint boundaries

| checkpoint | instances | codes | est. seconds |
|---|---|---|---|
| CP3-1 | 34 | 2,862,411 | 746.5 |
| CP3-2 | 45 | 2,871,276 | 748.8 |
| CP3-3 | 80 | 2,874,847 | 749.7 |
| CP3-4 | 104 | 2,874,546 | 749.6 |
| CP3-5 | 4,114 | 2,963,745 | 772.9 |

Chunks are formed by descending instance cost so that no single checkpoint
dominates. The run is checkpointed to
`results/phase3_sweep_checkpoint.jsonl` (one JSON record per instance,
flushed and `fsync`ed), and `run_phase3.py full --resume` continues from it.

### B.5 What constitutes a complete run

> A **COMPLETE RUN** is: every instance produced by the frozen X14 filter is run
> through E2, `oracle_cycle` and `oracle_flat`; all three distributions agree on
> **every** coefficient; $\lvert\mathscr C\rvert=(P+1)^{\sum a}$ equals the flat
> oracle's actual code count; $\#\{\dim=0\}=2^B$; all coefficients are
> nonnegative integers. Instances with $\lvert\mathscr C\rvert>10^5$ are **not**
> verified — they are recorded as enumerator-only, exactly as the frozen
> Blueprint's validation-scope amendment requires. The run is complete only when
> the recorded `n_verified` equals the number of sweep instances and
> `n_failed = 0`.

### B.6 Scope exclusions recorded now

- **X11–X13 are ring instances** (the affine algebra $A$ and $R_{m,q}$). The
  frozen Blueprint assigns them to CP6. They are **not** run in P3 and are
  recorded as deferred, not as passed.
- **X10a/X10b have $\lvert\mathscr C\rvert=279{,}936>10^5$**, so E3 cannot
  certify them. They are computed by E2 and **flagged ENUMERATOR ONLY**.
- No minimum distance (CP8), no moments (CP5), no asymptotics (CP7) — all
  outside P3.

---

## C. Shared-dependency analysis (E2 vs the oracles)

| shared component | used by E2 | used by oracles | verdict |
|---|---|---|---|
| `core/field.py` ($\mathbb F_q$ arithmetic) | no | yes | CP1-verified kernel; encodes no hull theory |
| `core/poly.py` (polynomial arithmetic over $\mathbb F_q$) | no | yes | same |
| `core/zz.py` (exact $\mathbb Z[z]$) | yes | no | same |
| `core/cycle_poly.py` ($W_a^{(P)}$) | **yes** | no | P2-verified; the object under test |
| `core/factor.py` `factor_xn_minus_lambda` | yes (via E1) | yes | **genuine shared dependency — identified** |
| `core/cycles.py` `hash_map` | yes (via E1) | **no** (`hash_map_independent`) | removed |
| `core/cycles.py` cycle construction | yes (via E1) | `oracle_cycle` only | **removed from `oracle_flat`** |
| transfer matrix $T_P$ | **yes** | no | never appears in either oracle |

The one genuinely shared dependency left is `factor_xn_minus_lambda`. It is not
part of the enumerator — it is the CP1 field/factorisation kernel, independently
verified in P1 against the factorisations printed in the sources and re-certified
inside `oracle/brute.py` on every instance via `certify_factorisation`. No hull-
dimension logic passes through it.

---

## D. Verdict

> **P3a = PASS** (10/10 instances, all regimes covered, 33,593 codes compared,
> zero mismatches).
>
> **P3b = PROCEED.** Estimated 113.7 min for 4,377 instances and 14,446,825
> codes, checkpointed in five chunks and resumable.
>
> **The estimate above is committed before the full run**, as the execution plan
> requires. P3c has **not** been executed at the time of this commit.
