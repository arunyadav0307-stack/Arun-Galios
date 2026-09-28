# AUDIT REPORT — `locked_blueprint.docx` vs frozen `BLUEPRINT.md` / `EXECUTION_PLAN.md`

**Audit type:** design-fidelity audit (read-only). No file was modified. The frozen
Blueprint and Execution Plan were not touched.

**Auditor note up front:** the document named in the request could not be found.
Sections A–F below therefore report only what is *verifiable from the frozen
references*, and the audit verdict is **BLOCKED — NOT ASSESSABLE**, not
"faithful" and not "not faithful". Section G gives the exact inputs required to
unblock. Section H.1 is a complete expected-value register so that the audit can
be re-run mechanically the moment the document is supplied.

---

## 0. Blocker: the document under audit is not present in this workspace

| Required input | Status | Evidence |
|---|---|---|
| `locked_blueprint.docx` | **ABSENT** | `find / -iname "*.docx"` → no results; `find / -iname "*locked_blueprint*"` → no results |
| reference PDF ("mechanics reference PDF") | **ABSENT** | `find / -iname "*.pdf"` → no results anywhere on the filesystem |
| `BLUEPRINT.md` @ `076f910` | **PRESENT, verified** | sha256 `2e02572a0ef970a93dcf242c725f290938d9b180787beaccd8acedf6a609214f` |
| `EXECUTION_PLAN.md` @ `96c7885` | **PRESENT, verified** | sha256 `6026536e3df10ac80a35225ddc42fc3b5f52f845be92da98531aca4bc8f65b32` |

Additional environment fact recorded for provenance: the sandbox had been
re-cloned and local git history reset to the single grafted commit `00fefba`.
The history was re-fetched from `origin` (authentication works) and the two cited
commits were confirmed to exist:

- `076f910` — "Blueprint v4/FINAL: close OPEN-18 on full-text evidence; freeze contributions"
- `96c7885` — "Add phase-by-phase execution roadmap (EXECUTION_PLAN.md)"

Both frozen documents in the working tree are **byte-identical** to the blobs at
those commits, so the reference side of this audit is sound. Only the audited
side is missing.

**No claim about the content of `locked_blueprint.docx` is made anywhere in this
report.** Nothing has been assumed, inferred or invented about it.

---

## A. PASS items

These are the checks that could actually be executed.

| # | Check | Result |
|---|---|---|
| A-1 | Frozen reference commits `076f910` and `96c7885` exist and are reachable | **PASS** |
| A-2 | Working-tree `BLUEPRINT.md` is byte-identical to `076f910:BLUEPRINT.md` | **PASS** (sha256 `2e02572a…`) |
| A-3 | Working-tree `EXECUTION_PLAN.md` is byte-identical to `96c7885:EXECUTION_PLAN.md` | **PASS** (sha256 `6026536e…`) |
| A-4 | Frozen Blueprint is internally consistent on the novelty boundary (§2.0.5/§2.0.6/§8.1 Classes A–E are mutually consistent) | **PASS** |
| A-5 | Frozen Blueprint's OPEN register and the Execution Plan's OPEN register are reconcilable (EP §F.3 ⊆ BP §11, plus execution-phase items OPEN-22…36 raised later) | **PASS** |
| A-6 | No scientific content was imported from any external reference during this audit (the reference PDF is absent; nothing was fetched) | **PASS** (vacuous) |

---

## B. MISMATCHES

**None assessable.** Every one of the eight requested checks compares the
document against the frozen references, and the document is absent. Recording
"no mismatches" here would be an unsupported claim, so it is not recorded.

One discrepancy *within the request itself* is worth flagging now, because it
will affect the audit when the document arrives:

- **The request asks for "hypotheses H/H1/H2".** The frozen Blueprint defines
  exactly **one** standing hypothesis, **(H)** `λ^(1+p^(e-k)) = 1` (§2.0.4),
  plus two standing *restrictions* in §3.1 — `r ∣ (1+p^(e-k))` and
  `gcd(n',r)=1` — which are recorded as Limitation **L2**, not as hypotheses
  with H-labels. There is no "H1" or "H2" anywhere in the frozen Blueprint.
  If `locked_blueprint.docx` introduces H1/H2 labels, those labels are **not in
  the frozen source** and must either be removed or explicitly justified as the
  document's own naming of the §3.1 restrictions. This must be checked, not
  assumed.

---

## C. MISSING CONTENT

| # | Item | Note |
|---|---|---|
| C-1 | `locked_blueprint.docx` | The audited document. Absent. |
| C-2 | Reference PDF | Absent, so check 6 (design-only contribution) cannot be executed at all. |

---

## D. OVERCLAIMS

**Not assessable** (document absent). Section H.1.7 below lists the exact
permitted and forbidden wording that must be matched against the document's
claim statements C1–C5 and its Class A–E table, so this check can be completed
without re-reading the Blueprint.

---

## E. STATUS ERRORS

**Not assessable** (document absent). Section H.1.3 gives the authoritative
status of P0–P3 that the document must reproduce, and H.1.4 gives the OPEN
register.

---

## F. REFERENCE-PDF CONTAMINATION

**Not assessable.** The reference PDF is not present in this workspace, so it is
impossible to determine whether any scientific content was imported from it.
This is recorded as **OPEN**, not as PASS. No contamination can be ruled in or
out. When the PDF is supplied, the check must be run against the frozen source
inventory (`sources/` contains only P1/P2/P3 extracted text) and the
"no invented data" rule.

---

## G. REQUIRED CORRECTIONS

No corrections are proposed to `locked_blueprint.docx`, because the document was
not available to audit. The corrections required are to the **audit inputs**:

1. **Supply `locked_blueprint.docx`** into the workspace (e.g. the repository
   root or `docs/`). Without it, checks 1–8 cannot run.
2. **Supply the reference PDF** if contamination check 6 is to be executed;
   otherwise mark that check permanently OPEN.
3. **Confirm the H1/H2 question** (see B): either the document uses only the
   frozen label **(H)** and the two §3.1 restrictions (L2), or it introduces
   H1/H2 of its own — in which case those labels need an explicit provenance
   note, since they are absent from the frozen Blueprint.
4. Re-run this audit using the register in H.1 as the checklist.

**No change is made automatically**, per instruction.

---

## H. FINAL VERDICT

### **AUDIT BLOCKED — NOT ASSESSABLE**

The frozen reference side is verified sound (A-1…A-6). The audited document is
absent, so no statement of design fidelity — positive or negative — can be
supported by evidence. The verdict is deliberately **not** "DESIGN-FAITHFUL WITH
CORRECTIONS": that would assert a comparison that was never performed.

---

# H.1 EXPECTED-VALUE REGISTER (the audit checklist)

Everything below is quoted or derived **verbatim** from the frozen documents.
Each row is a predicate the document must satisfy.

### H.1.1 Mathematical fidelity (check 1)

**Definitions (§3.2) — exactly five, no more, no fewer:**

| Def | Content that must appear |
|---|---|
| **Def 1 (cycle word)** | exponent word `u^(c)=(u_0,…,u_{a-1}) ∈ [0,P]^a`, `u_m` = exponent of `f_m` in `g(x)` |
| **Def 2 (local weight)** | `w_a^(P)(u) := Σ_{m=0}^{a-1} min{u_{m-1}, P-u_m}`, indices mod `a`, range `[0, ⌊aP/2⌋]` |
| **Def 3 (cycle polynomial)** | `W_a^(P)(z) := Σ_{u∈[0,P]^a} z^{w_a^(P)(u)} ∈ Z[z]`, `c_{a,i}^(P) = [z^i]W_a^(P)` |
| **Def 4 (hull enumerator)** | `N_{n,q,λ,k}(y) := Σ_C y^{dim hull_k(C)}` |
| **Def 5 (restricted cycle polynomial)** | for `a≥2`, `W̃_a^(P) := tr(T̃_P^a)` with `T̃_P := T_P∘(J-I)` (diagonal zeroed); `W̃_1^(P) := W_1^(P)`; encodes P3's sample space `C(n,q,λ)` |

**Theorems (§3.3) — exactly T1–T9, with these status tags frozen:**

| Thm | Statement | Frozen status tag |
|---|---|---|
| T1 | Additive decomposition `dim hull_k(C) = Σ_c ord_{j(c)}(q)·w_{a(c)}^(P)(u^(c))` | `[EXPECTED; verified in probe on 7 instances]` |
| T2 | Enumerator factorisation `N = Π_c W_{a(c)}^(P)(y^{ord_{j(c)}(q)})` | `[EXPECTED; verified in probe on 7 instances]` |
| T3 | Transfer matrix `W_a^(P) = tr(T_P(z)^a)`, `T_P[z][x,y]=z^{min{x,P-y}}`; rational GF denominator `det(I-tT_P(z))` of t-degree `P+1`; order-`(P+1)` recurrence; `deg W_a^(P)=⌊aP/2⌋` | `[VERIFIED-IN-PROBE for P≤4, a≤5]` |
| T4 | Semisimple closed form, `P=1`: `W_a^(1)=(1+√z)^a+(1-√z)^a = 2Σ_{i=0}^{⌊a/2⌋} C(a,2i) z^i` | `[VERIFIED-IN-PROBE, a≤6]` |
| T5 | LCD count `c_{a,0}^(P)=2`, extremal words `0^a, P^a`; hence `#{dim=0}=2^B` | `[VERIFIED-IN-PROBE; proof §4.2]` |
| T6 | Self-orthogonal / dual-containing counts via `S_≥`, `S_≤` | `[EXPECTED]` |
| T7 | Moments: `μ_a^(P)=(W')/(P+1)^a`, `v_a^(P)=(W''+W')/(P+1)^a−μ²` | `[EXPECTED]` |
| T8 | CLT, conditional | `[EXPECTED — CONDITIONAL, highest risk]` |
| T9 | Algorithm, polynomial in `log|C|`, `n`, `log q` | `[EXPECTED]` |

**Standing data (§3.1) that must be reproduced unchanged:**
`q=p^e`; `k∈{0,…,e-1}`; `n=n'p^ν`, `gcd(n',p)=1`; `λ∈F_q*`; `r=ord(λ)`;
`r ∣ (1+p^{e-k})`; `gcd(n',r)=1`; `P:=p^ν`; `σ=(·)^p`;
`f^#(x)=Σ f_0^{-p^{e-k}} f_i^{p^{e-k}} x^{m-i}`; the `l'`/`l`/`a_t`/`A`/`B_i`/`D_t`/
`β_t(j)`/`Δ_t`/`Σ_t a_t Δ_t = n'` indexing convention; `|C|=B=Σ_{t,j}β_t(j)`.

**Hypothesis:** exactly one, **(H)** `λ^{1+p^{e-k}} = 1` (§2.0.4), and **all of
T1–T6 must be stated under (H)**. Evidence quoted: 45 mismatches when (H) was
ignored, **all** at (H)-violating instances; 0 mismatches in 63 further instances
after imposing (H). Note that (H) is automatic for `k=0`, `λ=±1` (PA-1's scope).

**Scope items that must be present, not dropped:**
- `a ≥ 3` (cycle lengths beyond PA-1's `a∈{1,2}`)
- `P > 1` (repeated roots), including combined with `a≥3`
- general `λ` under (H)
- **even characteristic in scope** — the equal-degree-splitting defect is
  *repaired*; 209 even-characteristic factorisation instances across
  `q∈{2,4,8,16,32,64,256}` pass with 0 failures; grid rows X8 and X14-even are
  therefore valid
- the `ν>e` correction `σ^ν = σ^{ν mod e}` in the `p^ν`-th-root step
  `μ^{p^ν}=λ` (176 odd-characteristic instances pass, 0 failures)
- exclusions: `|C| > 10^5` cannot be E3-certified (enumerator only, flagged);
  ring instances X11–X13 not yet run

### H.1.2 Novelty boundary (check 2)

**PA-1 boundary must be preserved:** PA-1 (Sangwisut–Jitman–Ling–Udomkavanich,
*Finite Fields Appl.* 33 (2015) 232–257) reaches only `λ=±1`, Euclidean/Hermitian,
`a∈{1,2}`. Within that scope **our T1/T2/T5 are reformulations of PA-1 Thm 5 /
Thm 11 / Cor 12 / Rem 13(2)**, and T4 for `a≤2, P=1` reduces to PA-1 Eq. (26).
This must be stated, not hidden.

**Mean must NOT be a novelty claim.** §2.0.6 item 3 and Class E1: the first
moment is fully occupied — PA-10 Thm 10 (Euclidean cyclic), PA-12 Thm 3
(Euclidean negacyclic), PA-11 (Hermitian constacyclic, `ord(λ)∣(q+1)`), P3
(`k`-Galois). Only **second and higher moments** and the **distribution** may be
claimed. E1 is "Removed" from headline contributions.

**Novelty targets that must be preserved:** `a≥3` (B2), repeated-root `P>1`
combined with `a≥3` (B3), higher moments/variance (C3), distributional limit law
(C4), ring enumerators over `A` and `R_{m,q}` (C5), algorithmic whole-distribution
single-pass enumeration (D1).

**Permitted wording (verbatim, §8.1):**
- *"generalises Sangwisut–Jitman–Ling–Udomkavanich (2015) to $k$-Galois duality"*
- *"to the best of our knowledge, the first distributional limit law for a hull dimension"*
- *"the first enumeration over the affine algebra ring $A$"* (subject to OPEN-16)

**Forbidden wording (verbatim, §8.1):**
- *"we introduce the first enumerator for hulls"*
- *"a novel LCD count $2^B$"*
- *"we determine the average hull dimension"*
- *"a new generating function"*
- *"first asymptotic study of hull dimensions"*

C4 additionally **must** be worded as the first *distributional* limit law, never
as the first asymptotic result — PA-10 Thms 26–27 already give
`limsup/liminf` of the normalised **mean**.

**Class E items must remain removed/held:** E1 mean (removed), E2
Hermitian-case variance/limit law (removed pending OPEN-20), E3 any claim resting
on absence of prior art (held pending OPEN-16 cited-by crawl).

### H.1.3 Validation fidelity (check 3)

| Phase | Required status in the document |
|---|---|
| **P0** | **PASS conditional** |
| **P1** core implementation | **PASS** |
| **P1** new 5,670-instance validation suite | **PASS** (described only as the newly defined validation suite, **not** as recovery of the original 385 suite) |
| **P1** original 385-instance suite | **OPEN / NOT RECONSTRUCTIBLE** |
| **P2** | **PASS** (32 checks, 0 failures) |
| **P3** | **OPEN unless actual P3 evidence exists** — actual P3 evidence **does** exist: `results/phase3_summary.json`, 4,377/4,377 verified, 0 failed, 29,803 coefficients and 14,446,825 codes compared, 0 mismatches, max \|coeff diff\| = 0, peak RSS 65.7 MB, total 11,016.5 s. So P3 may be reported as **PASS** with that evidence cited. It must not be upgraded further. |

**No status may be upgraded without evidence.** Instance categories must stay
distinct: (a) inherited/probe-derived — NOT RECONSTRUCTIBLE / OPEN-22, never used
as an acceptance gate; (b) newly constructed validation — PASS 5,670/5,670;
(c) published-paper benchmark instances. These three must never be mixed.

### H.1.4 OPEN register (check 4)

Every item below must be preserved. Statuses as of the latest evidence:

| ID | Item | Status |
|---|---|---|
| **OPEN-22** | Original 385-instance suite not reconstructible | **OPEN / NOT RECONSTRUCTIBLE** — must **not** be described as resolved |
| **OPEN-28** | Provenance label for newly chosen odd-`q` inputs in `regression_suite_v1` | **retained** |
| **OPEN-29** | CP1-3/CP1-4 shape-comparison only (primitive elements in F_25/F_81 not recoverable; λ determined up to Frobenius conjugacy) | **OPEN** (owner P4) |
| **OPEN-30** | `galois` fails **flakily** on `x^21+1` over F_16/F_256; 1–2 oracle tests skipped depending on run | **OPEN** — oracle limitation, ours not implicated |
| **OPEN-31** | Global enumerator `N` (T2) not assembled/checked against brute force | **CLOSED with evidence** (P3: 4,377 instances, two independent oracles, 0 mismatches) |
| OPEN-23/24/25/26/27 | `instances.yaml`, `targets.yaml`, `SOURCE_INVENTORY.md`, `references.bib`, untested `setup_env.sh` | carried |
| OPEN-32/33 | §3.5 structural reductions (Hadamard/rank-one form; `M=LDL^T`, `T_P=MR`, `tr(T_P^a)=tr((DN)^a)`; Hankel form `N[i,j]=(P+1-i-j)^+`) — marked `[EXPECTED — needs verification]`, **not verified in P2** | **OPEN** (P2 extension / P7) |
| OPEN-34 | `T_6` self-orthogonal/dual-containing counts | **OPEN** (P7) |
| OPEN-35 | Restricted/unrestricted mean difference (§3.6, `P=3`: `34/16` vs `24/12`) not recomputed | **OPEN** (P6) |
| OPEN-36 | `T_8` CLT untouched | **OPEN** (P8) |
| OPEN-37 | X10a/X10b `|C|=279,936 > 10^5` — over the frozen E3 bound, **enumerator only**, never citable as validated | **OPEN** |
| OPEN-38 | Ring instances X11–X13 not run | **OPEN** (CP6) |
| OPEN-39 | `λ` sampled **by multiplicative order**, not exhaustively | **OPEN** |
| OPEN-16/17/19/20 | Cited-by crawl of PA-1 not run; PA-2 §3–§5 paywalled; second prior-art search not run; PA-11 §3–§6 paywalled | carried from the frozen Blueprint §11 / EP §F.3 |
| OPEN-00/01/02/04 | Workflow v6 doc missing; no Magma; no SageMath; no quartile metrics | carried |

### H.1.5 No invention (check 5)

Every number in the document must trace to one of:
- the frozen Blueprint (sha256 `2e02572a…`) — e.g. 45 pre-(H) mismatches, 63
  post-(H) instances, 40/40 PA-1 Table 2 rows / 529 cells, 40 `a≥3` instances,
  23 `a≥3 ∧ P>1` instances, 209 even-char instances, 176 odd-char instances,
  385 internal factorisation instances, 32,768-code instance, `a=1..12` closed
  form, symbolic agreement `P=1..7`;
- the frozen Execution Plan (sha256 `6026536e…`);
- a committed artefact under `results/` with a sha256 — e.g. P3's
  `phase3_summary.json` (4,377 / 29,803 / 14,446,825 / 0).

Any number with no such trace is **invented** and must be struck.

### H.1.6 Registers (check 8)

The frozen documents contain no artefact literally titled "equation register",
"table register", "figure register", "validation register", "reproducibility
register" or "production-plan register". The corresponding authoritative content
lives in:

| Register requested | Authoritative frozen location |
|---|---|
| Equation register | BP §3.1 standing data, §3.2 Defs 1–5, §3.3 T1–T9, §4.1–4.7 proof sketches, Annex A |
| Table register | BP §0.1, §2.0.3, §2.0.8 prior-art matrices; §5.5 experiment grid X1–X15; §8.1 Classes A–E; §9 venue table |
| Figure register | BP §10 paper structure (no figures are pre-specified; Q–Q plots are P8, pre-registered per EP XC-20) |
| Validation register | BP §6.2 CP0–CP8 gates; EP §C PHASE ACCEPTANCE MATRIX; EP §D CROSS-CHECK MATRIX (XC-1…XC-20) |
| Reproducibility register | BP §6.0 environment, §6.1 repository layout; EP §E FINAL DELIVERABLES TREE; EP §G.3 |
| Production plan | EP §P0–§P11 phase definitions, §B dependency graph and schedule, §G DEFINITION OF "PAPER COMPLETE" |

**CP3 gate wording that must be preserved exactly:** *"For every instance in
§5.5 with `|C| ≤ 10^5`: E2 output ≡ E3 histogram, and `|C|=(P+1)^{Σ a(c)}`, and
`#{dim=0}=2^B`"*; plus the gate *"≥300 verified instances across the X14 sweep
with zero mismatches"*.

**X14 sweep definition that must be preserved:** all `(q,n,λ,k)` with
`q∈{4,8,9,16,25,27,32,49,64,81,121,125,128,256}`, `n≤60`, `gcd(n',r)=1`,
`r ∣ (1+p^{e-k})`, as many as E3 can certify.

### H.1.7 Claims C1–C5 (check 7)

The frozen Blueprint's §8.1 Class C table is the authority. C1–C5 as the document
states them must match:

| # | Frozen claim | Overclaim trap to check for |
|---|---|---|
| **C1** | Transfer-matrix/trace structure `W_a^(P)=tr(T_P^a)`, `T_P[x,y]=z^{min{x,P-y}}` | Must state that **for `a≤2` it agrees with PA-1 and the paper must say so**; load-bearing only for `a≥3` |
| **C2** | Closed form `W_a^(1)=(1+√z)^a+(1-√z)^a` **for `a≥3`** | Must state it **reduces to PA-1 Eq. (26) for `a≤2`**; verified `a=1..12` |
| **C3** | Second and higher moments / variance | Must **not** extend to the mean; must note PA-10/PA-12 read in full with no second moment |
| **C4** | Distributional limit law (conditional CLT) at fixed `n`; `Θ(n)` vs `O(1)` contrast with Sendrier (PA-6) | Must be worded as first **distributional** limit law; must distinguish from PA-10 Thms 26–27 |
| **C5** | Ring-level enumerators over `A` (P2's open problem) and `R_{m,q}` | Must note PA-4/5/7/8 cover chain rings, Z_4, a non-chain ring (dimensions only), double-cyclic — **none covers `A`** |

Also verify Class B1–B5 are framed as **generalisation, never as first**, and
Class A1–A10 appear as **known/reformulated, explicitly not claimed as novelty**.
Class D1 must keep the single-pass vs PA-1 Money-Changing `O((s+t)|h(ℓ)|)` per
dimension contrast.

---

## Summary

The frozen reference side of this audit is **verified sound and byte-identical**.
The audited document `locked_blueprint.docx` **is not present in this
workspace**, and neither is the reference PDF. Therefore:

- checks 1–5, 7, 8 could not be executed against the document;
- check 6 could not be executed at all;
- **no mismatch, overclaim, status error or contamination is asserted**, because
  none was observed — nothing was available to observe.

**Verdict: AUDIT BLOCKED — NOT ASSESSABLE.** Re-run once the document (and,
ideally, the reference PDF) is supplied; section H.1 is the complete checklist.
