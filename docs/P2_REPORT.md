# P2 REPORT — Cycle-Polynomial / Transfer-Matrix Foundation

**Phase:** P2 (of `EXECUTION_PLAN.md` @ `96c7885`)
**Frozen source:** `BLUEPRINT.md` @ `076f910` (FINAL)
**Builds on:** P0 `076f977` · P1 `7d46765` (closure-audited)

> **Checkpoint bookkeeping (factual, no history rewritten).** The annotated,
> pushed tag `p1` points at commit `1a70750`, the P1 *implementation*
> checkpoint. The P1 **closure audit** commit `7d46765` -- the final, accepted P1
> state -- sits directly on top of that tag and is what this report builds on.
> The published tag was **not** moved; the history is left as-is and this note
> records the relationship so a reader checking out tag `p1` knows it is one
> commit short of the audited state.

**Executed:** 2026-09-27
**Branch:** `arena/01a0e2d5-arun-galios`
**Scope:** the cycle-local enumeration framework — local weights, cycle
polynomials, the transfer matrix and its trace representation, the $P=1$
closed form, and the PA-1 boundary. **No P3+ work.**

> **Gate result: P2 = PASS.** 32/32 verification checks pass; 191 new unit
> tests pass, taking the whole suite to 566 passed / 5 skipped /
> 0 failed (375 of those passing tests come from P1). Every core identity is
> derived and then
> confirmed both by an independent brute-force enumeration of local states and
> by sympy over $\mathbb Q[z]$.
>
> **Status distinctions carried forward from P1, unchanged:**
> original 385-instance suite = **OPEN / NOT RECONSTRUCTIBLE** (OPEN-22);
> the 5 670-instance suite = **validated new suite** (category (b));
> **OPEN-28 / OPEN-29 / OPEN-30 remain open.**

---

## A. P2 MATHEMATICAL RESULTS

### A.1 What was built

| Module | Lines | Contents |
|---|---|---|
| `core/zz.py` | 266 | Exact $\mathbb Z[z]$ polynomials and matrices (no floats): add/mul/pow, determinant by Laplace, characteristic polynomial via Newton's identities on $\operatorname{tr}(A^i)$ |
| `core/cycle_poly.py` | 304 | Definition 2 (local weight), Definition 3 (cycle polynomial by word enumeration), Definition 5 (restricted), T3 (transfer matrix, bivariate GF, $\det(I-tT)$, recurrence), T4 (closed form), degree/LCD statements |
| `core/hull_dim.py` | 151 | T1 verification instrument: $h=(x^n-\lambda)/g$, $h^\#$, $\operatorname{lcm}$, $\dim=n-\deg\operatorname{lcm}$ |
| `tests/test_zz.py` | 97 | $\mathbb Z[z]$ kernel (6 tests) |
| `tests/test_cycle_poly.py` | 385 | all P2 identities, edge cases, PA-1 boundary, T1 (185 tests) |
| `scripts/run_phase2.py` | 626 | full verification driver with machine-readable outputs |

**Scope discipline.** `core/hull_dim.py` computes a *single* code's hull
dimension from polynomials. It is **not** Algorithm E2 (the enumerator $N$ that
multiplies the $W$ factors and histograms over all codes) — that is P3's work
and is deliberately absent. No moments (P6), no limit law (P8).

### A.2 The frozen hull-dimension identity (task 1)

**T1** (Blueprint §3.3, proof §4.1):
$$\dim\operatorname{hull}_k(C)=\sum_{\mathfrak c\in\mathfrak C}\operatorname{ord}_{j(\mathfrak c)}(q)\cdot w_{a(\mathfrak c)}^{(P)}\big(u^{(\mathfrak c)}\big).$$

The frozen proof's engine is the pointwise identity
$$P-\max\{u_m,P-u_{m-1}\}=\min\{P-u_m,u_{m-1}\},$$
which I verified exhaustively for $P=1..6$ and all $u_m,u_{m-1}\in[0,P]$ -- that is
$\sum_{P=1}^{6}(P+1)^2=154$ ordered pairs, all holding. Since $\min$ is
commutative the RHS **is** $\phi(u_{m-1},u_m)=\min\{u_{m-1},P-u_m\}$, the frozen
local weight of Definition 2. Then
$$\sum_m\big(P-\max\{u_m,P-u_{m-1}\}\big)=aP-\sum_m\max\{u_m,P-u_{m-1}\}=aP-\deg\operatorname{lcm}\big|_{\mathfrak c},$$
$\deg f_m=d$ is constant along a cycle, $n=\sum_{\mathfrak c}d(\mathfrak c)a(\mathfrak c)P$, and subtracting gives T1.

**Verified computationally, not just argued:** on 12 instances × 3 choices of
$\lambda$ × 7 word assignments = **126 checks**, the polynomial route
($n-\deg\operatorname{lcm}(g,h^\#)$) and the local-weight route agree in
**126/126**. 14 further (instance, $\lambda$) pairs were **skipped** because
hypothesis **(H)** $\lambda^{1+p^{e-k}}=1$ fails there, so $\#$ does not permute
the irreducible factors — exactly as the Blueprint's §2.0.4 predicts.

### A.3 Local weight for a single Galois-reciprocal cycle (task 2)

**Definition 2.** For $a\ge1$, $P\ge1$, $u\in[0,P]^a$ (indices mod $a$):
$$w_a^{(P)}(u):=\sum_{m=0}^{a-1}\min\{u_{m-1},P-u_m\},\qquad w_a^{(P)}(u)\in\big[0,\lfloor aP/2\rfloor\big].$$

The decisive structural observation (Blueprint §2.3) is that this is a
**nearest-neighbour cyclic interaction**:
$$w_a^{(P)}(u)=\sum_{m=0}^{a-1}\phi(u_{m-1},u_m),\qquad \phi(x,y)=\min\{x,P-y\}.$$
This is verified directly in the tests (`test_local_weight_is_a_cyclic_nearest_neighbour_sum`).
The range bound $\lfloor aP/2\rfloor$ is verified for $P=1,2,3,4,5$ and
$a=1..5$, and the degree claim $\deg W_a^{(P)}=\lfloor aP/2\rfloor$ is verified
on 21 $(P,a)$ pairs with an explicit maximiser exhibited.

For $a=2$ the weight collapses to a closed form, derived and verified:
$$w_2^{(P)}(v,w)=\min\{v+w,\,2P-v-w\},$$
by the case split $v+w\le P$ (then $\min\{v,P-w\}=v$, $\min\{w,P-v\}=w$) versus
$v+w>P$ (then both mins become $P-w$ and $P-v$).

### A.4 Derivation of $W_a^{(P)}(z)=\operatorname{tr}(T_P(z)^a)$ (task 3)

**Definition 3.** $W_a^{(P)}(z):=\sum_{u\in[0,P]^a}z^{w_a^{(P)}(u)}\in\mathbb Z[z]$.

**Derivation (T3).** Expand the trace as a sum over closed walks:
$$\operatorname{tr}(T_P^a)=\sum_{u_0,\dots,u_{a-1}\in[0,P]}\prod_{m=0}^{a-1}T_P[u_{m-1},u_m]$$
with cyclic indexing $u_{-1}\equiv u_{a-1}$. Since
$T_P[x,y]=z^{\min\{x,P-y\}}=z^{\phi(x,y)}$, each summand is
$$z^{\sum_m\phi(u_{m-1},u_m)}=z^{w_a^{(P)}(u)},$$
and the sum over all $(u_0,\dots,u_{a-1})\in[0,P]^a$ is exactly $W_a^{(P)}(z)$ by
Definition 3. $\square$

The $(P+1)^a$ summands of the trace expansion are in bijection with the
$(P+1)^a$ words of Definition 3 — this is the whole content of the theorem, and
it is why no approximation or limit is involved.

**Computationally verified** on 24 $(P,a)$ pairs: brute-force enumeration of all
local states (which never mentions a matrix) equals the trace (which never
enumerates words). And independently by sympy over $\mathbb Q[z]$ on 16 pairs.

### A.5 $T_P$ defined precisely (task 4)

$$T_P(z)=\big(z^{\min\{x,P-y\}}\big)_{0\le x,y\le P}\in\mathbb Z[z]^{(P+1)\times(P+1)}.$$

- **Index set:** $\{0,1,\dots,P\}$ — the possible exponents $u_m$ of a single
  irreducible factor, i.e. the state space of Definition 1.
- **Row index** $x$ = the *previous* state $u_{m-1}$; **column index** $y$ = the
  *current* state $u_m$.
- **Exponent entry** $\min\{x,P-y\}=\phi(x,y)$, the frozen local weight term.
- Shape $(P+1)\times(P+1)$; verified entry-by-entry for $P=0..5$.
- At $z=1$, $T_P(1)=J_{P+1}$ (all-ones), hence $\det(I-tT_P(1))=1-(P+1)t$ —
  verified for $P=0..5$.

### A.6 Product decomposition over independent cycles (task 5)

**T2** (Blueprint §3.3, proof §4.3):
$$N_{n,q,\lambda,k}(y)=\prod_{\mathfrak c\in\mathfrak C}W_{a(\mathfrak c)}^{(P)}\big(y^{\operatorname{ord}_{j(\mathfrak c)}(q)}\big).$$

**Derivation.** (i) P1's factorisation of $x^n-\lambda$ is into *pairwise
coprime* monic irreducibles, so $C\mapsto\big(u^{(\mathfrak c)}\big)_{\mathfrak c\in\mathfrak C}$
is a **bijection** $\mathscr C_{n,q,\lambda}\to\prod_{\mathfrak c}[0,P]^{a(\mathfrak c)}$.
(ii) By T1 the hull dimension is a **sum over cycles** of a function of that
cycle's own word alone. (iii) Summing a product of independent factors gives the
product of the sums, each sum being $W_{a(\mathfrak c)}^{(P)}$ with $z$ replaced
by $y^{\operatorname{ord}_{j(\mathfrak c)}(q)}$. $\square$

**Corollaries verified:** $N(1)=|\mathscr C_{n,q,\lambda}|=(P+1)^{\sum_{\mathfrak c}a(\mathfrak c)}$
— this is the code count P1's Algorithm E1 predicts, so the two phases are
consistent; and $\#\{C:\dim\operatorname{hull}_k(C)=0\}=2^{B}$ (T5).

**Scope note.** The *global* product $N$ is assembled and histogrammed in P3
(Algorithm E2). P2 establishes and verifies the *per-cycle* factors $W_a^{(P)}$
and the decomposition argument; the end-to-end enumerator against a brute-force
oracle is P3's gate.

### A.7 Arbitrary validated cycle length $a$, especially $a\ge3$ (task 6)

$W_a^{(P)}$ is computed for arbitrary $a$ with no assumption on $a$: the trace
formula is uniform in $a$, and brute force was run up to $a=12$ at $P=1$,
$a=5$ at $P=2$, $a=4$ at $P=3,4$, $a=3$ at $P=5$.

The $a\ge3$ regime is where the $a$ exponents become **cyclically coupled**:
$w_a$ is a cyclic chain of $\phi$-interactions, and no variable can be
eliminated the way PA-1 eliminates one variable at $a=2$. This coupling is
precisely what the transfer matrix handles.

### A.8 Distinction from PA-1 (task 7) — no unsupported novelty claim

**What PA-1 already has (Class A, explicitly not claimed):**
per the frozen §2.0.1, for $\lambda=\pm1$ with the Euclidean or Hermitian inner
product — PA-1's entire scope — our T1, T2 and T5 are **equivalent
reformulations** of PA-1's Theorem 5 / Theorem 11 / Corollary 12 / Remark 13(2),
and our T4 collapses to PA-1's Eq. (26). **The transfer-matrix/trace idea is
therefore NOT claimed as new in general.** The frozen Class C1 wording is
respected: *"For $a\le2$ it agrees with PA-1 and the paper must say so; the
structure is load-bearing precisely for $a\ge3$."*

**Verified boundary, computationally:**
- $W_1^{(P)}=\sum_{b}\|\{b,P-b\}\|z^b$ = PA-1's per-cycle object, **$P=1..7$, exact**.
- $W_2^{(P)}=\sum_{a=0}^{P}2^{1-\lfloor a/P\rfloor}(a+1)z^a$ = PA-1's per-cycle
  object, **$P=1..7$, exact**.
- For $k=0$ every cycle has length $\le2$ on **8/8** instances tested — because
  for $k=0$, $j=e$ and $\sigma^e=\mathrm{id}$, so $f^\#=f^*$ and $(f^*)^*=f$.
  This is the *structural reason* PA-1's involutions force $a\le2$, and it is
  verified by computation, not asserted.
- For $k=3$ ($k\notin\{0,e/2\}$) cycle lengths $\ge3$ do occur (lengths seen:
  $\{1,4\}$).

**The genuine boundary (frozen §2.0.1, restated):** PA-1's Eq. (26) is a
product of **independent one-parameter** sums; for $a\ge3$ the $a$ exponents are
cyclically coupled and that product form cannot express the count. The
transfer matrix is the object that can.

### A.9 The $P=1$ specialisation (task 8)

**T4** (Blueprint §3.3, derivation §4.5): for $P=1$,
$$W_a^{(1)}(z)=(1+\sqrt z)^a+(1-\sqrt z)^a=2\sum_{i=0}^{\lfloor a/2\rfloor}\binom{a}{2i}z^i.$$

**Derivation.** $T_1(z)=\begin{pmatrix}1&1\\z&1\end{pmatrix}$; characteristic
polynomial $\lambda^2-2\lambda+(1-z)$ (verified in the tests); eigenvalues
$1\pm\sqrt z$; $\operatorname{tr}(T_1^a)=(1+\sqrt z)^a+(1-\sqrt z)^a$; the
binomial expansion keeps only even powers of $\sqrt z$, giving the even-part sum.

**Verified** against the trace for $a=1..14$, and the coefficient lists
$[2],[2,2],[2,6],[2,12,2],[2,20,10],[2,30,30,2]$ for $a=1..6$ match the frozen
Blueprint §2.4 **verbatim**. Also confirmed symbolically by sympy for $a=1..10$.

**Exactly which previously known results it reduces to.** At $P=1$:
- $W_1^{(1)}=2$ and $W_2^{(1)}=2(1+z)$, hence
  $$N(y)=2^{B}\prod_{\mathfrak c:\,a(\mathfrak c)=2}\big(1+y^{d(\mathfrak c)}\big).$$
- PA-1's **Corollary 12** states $\#\{\dim=\ell\}=2^{s+t}|h(\ell)|$ where $|h(\ell)|$
  is the coefficient of $X^\ell$ in their **Eq. (26)**. At $P=1$ that equation's
  $\chi=0$ factors collapse to the single term $1$ (because $\lfloor p^\nu/2\rfloor=0$)
  and its $\chi=1$ factors are exactly $(1+X^{\operatorname{ord}_j(q)})$.
- Matching term by term: $s$ = number of length-1 cycles, $t$ = number of
  length-2 cycles, $B=s+t$, and the $\chi=1$ factors correspond one-to-one with
  the length-2 cycles. **The two expressions are identical.**

So at $P=1$ with $a\le2$ our framework **reproduces** PA-1 Cor 12 exactly — it
does not extend it there. That is stated plainly, per the frozen Class A rule.

### A.10 $P>1$ handled only within the frozen hypotheses (task 9)

$P=p^\nu>1$ arises only when $p\mid n$, and then every irreducible factor of
$x^n-\lambda$ occurs with multiplicity exactly $P$. All $P>1$ computations were
run under the frozen standing hypotheses: $n=n'p^\nu$, $\gcd(n',p)=1$,
$r=\operatorname{ord}(\lambda)$, $r\mid(1+p^{e-k})$, $\gcd(n',r)=1$, and
**(H)** $\lambda^{1+p^{e-k}}=1$. Where (H) fails, the computation is **skipped
and recorded as skipped**, never forced.

Verified: $P=2$ with $a=3$ ($W=z^3+18z^2+6z+2$), $P=3,4,5$ with $a\le3,4,3$;
$\deg W=\lfloor aP/2\rfloor$ and $c_0=2$ hold for all of them; and
$\mu=\sigma^{e-(\nu\bmod e)}(\lambda)$ with $\mu^{P}=\lambda$ is asserted
internally by P1's factoriser (defect (b) repair).

### A.11 Symbolic verification (task 10)

Every algebraic identity was checked symbolically with sympy, **independently of
the $\mathbb Z[z]$ kernel**:

| Identity | Range | Result |
|---|---|---|
| $\operatorname{tr}(T_P(z)^a)$ over $\mathbb Q[z]$ | 16 $(P,a)$ pairs | 16/16 |
| $W_a^{(1)}=(1+\sqrt z)^a+(1-\sqrt z)^a$ | $a=1..10$ | 10/10 |
| $\mathcal G_P(z,t)=\operatorname{tr}(tT_P(I-tT_P)^{-1})=\sum_{a=1}^{7}W_at^a$ | $P=1,2,3$ | 3/3 |
| $\deg_t\det(I-tT_P(z))=P+1$ | $P=1..5$ | 5/5 |
| $\operatorname{tr}(\tilde T_P^a)$ (Definition 5) | 5 pairs | 5/5 |
| Cayley–Hamilton recurrence | $P=1,2,3$ | 3/3 |

The recurrence was derived **from** the characteristic polynomial via
Cayley–Hamilton and Newton's identities on $\operatorname{tr}(T_P^i)$ — so the
order-$(P+1)$ statement is proved, not merely observed.

### A.12 The $a\ge3$ example, independently constructed (task 12)

Three independent constructions of the same object:

**(i) Purely combinatorial.** $P=1$, $a=4$: enumerate all $2^4=16$ words
$u\in\{0,1\}^4$, apply Definition 2 term by term, histogram. The 16 local states
and their weights:

| weight $w$ | 0 | 1 | 2 |
|---|---|---|---|
| # words | 2 | 12 | 2 |

The two weight-0 words are $0000$ and $1111$ (T5's extremal words); the two
weight-2 words are $0101$ and $1010$. This gives
$W_4^{(1)}(z)=2+12z+2z^2$.

**(ii) Transfer matrix.** $\operatorname{tr}(T_1^4)=2+12z+2z^2$. **Agrees.**

**(iii) Closed form.** $(1+\sqrt z)^4+(1-\sqrt z)^4=2+12z+2z^2$. **Agrees.**

Also verified at $P>1$: $P=2$, $a=3$ gives $W=z^3+18z^2+6z+2$ by both routes.

**(iv) A real $a=4$ cycle from the P1 core.** $\mathbb F_{81}$, $n=5$, $k=3$,
$\lambda=1$ produces cycle shapes $[(1,1),(4,1)]$ — a genuine 4-cycle. On that
instance T1 was verified on 7 word assignments (7/7 agree), i.e. the hull
dimension computed from actual polynomials $g,h,h^\#$ equals
$d\cdot w_4^{(1)}(u)$ for the 4-cycle.

### A.13 Edge cases (task 13)

| Case | Result |
|---|---|
| $a=1$ (self-reciprocal) | $W_1^{(P)}=\sum_x z^{\min\{x,P-x\}}$; verified $P=1..5$; matches PA-1 |
| $a=2$ (reciprocal pair) | $w_2=\min\{v+w,2P-v-w\}$; $W_2^{(P)}$ matches PA-1, $P=1..7$ |
| $P=1$ (semisimple) | closed form T4 verified $a=1..14$; reduces to PA-1 Cor 12 |
| $a\ge3$ | verified for $P=1$ ($a\le12$), $P=2$ ($a\le5$), $P=3,4$ ($a\le4$), $P=5$ ($a\le3$) |
| $P=0$ (degenerate, single state) | $W_a^{(0)}=1$; handled without error, recorded not asserted |
| restricted, odd $a$, $P=1$ | $\tilde W_a^{(1)}=0$ (a 2-symbol cycle with no adjacent equals is impossible for odd $a$) — verified on both routes |
| $\operatorname{tr}(\tilde T_P)$ | $=0$, hence Definition 5's explicit $\tilde W_1:=W_1$ — honoured and tested |

All 12 combined edge cases satisfy: trace == brute force, $W(1)=(P+1)^a$,
$c_0=2$, $\deg=\lfloor aP/2\rfloor$.

### A.14 Mismatches encountered, and how they were handled (task 14)

Four checks **failed on the first run**. Each was diagnosed before any change
was made. **All four were defects in my verification code, not in the frozen
formulas.** Recorded here in full rather than quietly fixed:

| # | Failed check | Diagnosis | Resolution |
|---|---|---|---|
| 1 | $\det(I-tT_1(1))=1-2t$, got $[1,-2,0]$ | The top $t$-coefficient of $\det(I-tT_P(z))$ is a power of $(1-z)$ and **vanishes at $z=1$**; my comparison list did not trim trailing zeros. Verified $\det(I-tT_1(z))=1-2t+t^2(1-z)$ exactly. | Representation artefact; trim after evaluating. Formula correct. |
| 2 | $W_2^{(P)}$ vs PA-1's $\sum_a2^{1-\lfloor a/P\rfloor}(a+1)z^a$ | I summed $a=0..2P$, producing fractional coefficients ($1.5$ at $a=2P$). The correct range is $a=0..P$, since $\deg W_2^{(P)}=\lfloor 2P/2\rfloor=P$. With the correct range it matches **$P=1..7$ exactly**. | My summation range was wrong; the frozen formula is confirmed. |
| 3 | pointwise identity $P-\max\{u,P-v\}=\min\{u,P-v\}$ | I paired $\max\{u_m,P-u_{m-1}\}$ with $\min\{u_m,P-u_{m-1}\}$ — same order, which is **false**. The frozen §4.1 pairs it with $\min\{P-u_m,u_{m-1}\}$, which holds. | My variable order was wrong; the frozen identity holds. |
| 4 | sympy closed form | My reduction of $(1+\sqrt z)^a+(1-\sqrt z)^a$ was half-finished (dead code, wrong substitution). | Rewrote as an exact reduction in $\mathbb Q(w)[z]$ with $w^2=z$; matches $a=1..10$. |

Two further failures appeared in the unit tests and were likewise my assertion
errors, both diagnosed by hand before fixing:
- `local_weight((0,1,0),1)`: I asserted 2; the term-by-term evaluation gives
  $0+0+1=1$. Test corrected to 1 with the arithmetic written out.
- $\operatorname{lcm}(x+1,h)$ over $\mathbb F_3$, $n=3$, $\lambda=2$: I asserted
  $g,h$ are coprime. They are not — $x^3+1=(x+1)^3$ by Freshman's dream, so
  $h=(x+1)^2$, $\gcd=x+1$, $\operatorname{lcm}=(x+1)^2$. The **computed lcm was
  correct**; my comment "coprime here" was wrong. Test now asserts the correct
  non-coprime behaviour and adds a genuinely coprime contrast case.

**No frozen formula was altered, and no mismatch was repaired by assumption.**

---

## B. THEOREMS / LEMMAS ESTABLISHED

**Lemma 2.1 (nearest-neighbour form).** $w_a^{(P)}(u)=\sum_m\phi(u_{m-1},u_m)$
with $\phi(x,y)=\min\{x,P-y\}$. *Verified by direct computation for
$P=1..4$, $a=1..5$, all words.*

**Lemma 2.2 (range).** $0\le w_a^{(P)}(u)\le\lfloor aP/2\rfloor$, and the upper
bound is attained. *Verified for $P=1..5$, $a=1..5$; a maximiser is exhibited by
brute force.*

**Lemma 2.3 (extremal words).** $w_a^{(P)}(u)=0$ iff $u=0^a$ or $u=P^a$.
*Verified exhaustively for $P=1,2,3$, $a=1..5$.* Proof as frozen in §4.2: if
every term vanishes and some $u_i\ne P$, then $m=i$ forces $u_{i-1}=0$; then
$u_{i-1}\ne P$ (as $P\ge1$) so $m=i-1$ forces $u_{i-2}=0$; induct around the
cycle to get $u\equiv0$. Otherwise $u\equiv P$.

**Lemma 2.4 (two-state closed form).** For $a=2$,
$w_2^{(P)}(v,w)=\min\{v+w,2P-v-w\}$. *Proved by the case split on $v+w\le P$;
verified for $P=1,2,3$, all $v,w$.*

**Theorem 2.5 (T3, trace representation).**
$W_a^{(P)}(z)=\operatorname{tr}(T_P(z)^a)$ with
$T_P(z)=(z^{\min\{x,P-y\}})_{0\le x,y\le P}$. *Proof in §A.4. Verified on 24
$(P,a)$ pairs by brute force and 16 pairs by sympy.*

**Corollary 2.6 (bivariate GF).**
$\mathcal G_P(z,t)=\sum_{a\ge1}W_a^{(P)}t^a=\operatorname{tr}(tT_P(I-tT_P)^{-1})$
is rational with $\deg_t\det(I-tT_P(z))=P+1$. *Verified: series agreement for
$P=1,2,3$ to order $t^7$; determinant degree for $P=0..5$ (ours) and $P=1..5$
(sympy).*

**Corollary 2.7 (recurrence).** $W_a^{(P)}$ satisfies an order-$(P+1)$ linear
recurrence in $a$ over $\mathbb Z[z]$, obtained from the characteristic
polynomial of $T_P$ by Cayley–Hamilton. *Proved (not observed) and verified for
$P=0..5$.*

**Corollary 2.8 (degree).** $\deg W_a^{(P)}=\lfloor aP/2\rfloor$. *Verified on
21 $(P,a)$ pairs.*

**Corollary 2.9 (LCD count, T5).** $c_{a,0}^{(P)}=2$, hence
$\#\{C:\dim\operatorname{hull}_k(C)=0\}=2^B$. *Verified on 13 $(P,a)$ pairs.*

**Theorem 2.10 (T4, semisimple closed form).** For $P=1$,
$W_a^{(1)}(z)=(1+\sqrt z)^a+(1-\sqrt z)^a=2\sum_{i=0}^{\lfloor a/2\rfloor}\binom{a}{2i}z^i$.
*Proof in §A.9; verified $a=1..14$ against the trace, $a=1..10$ by sympy, and
$a=1..6$ against the frozen Blueprint's verbatim coefficient lists.*

**Definition 5 implementation (restricted).** $\tilde W_a^{(P)}=\operatorname{tr}(\tilde T_P^a)$
for $a\ge2$ with $\tilde T_P=T_P\circ(J-I)$, and $\tilde W_1^{(P)}:=W_1^{(P)}$.
*Both branches verified against brute force over words with no cyclically
adjacent equal exponents; the $a=1$ special case is honoured because
$\operatorname{tr}(\tilde T_P)=0$.*

**Theorem 2.11 (T1, hull-dimension identity).** *Statement and proof in §A.2;
verified on 126 (instance, word) checks against the polynomial route.*

**Proposition 2.12 (PA-1 boundary, structural).** For $k=0$, $j=e$, so
$\sigma^e=\mathrm{id}$ and $f^\#=f^*$ with $(f^*)^*=f$; hence every cycle has
length $\le2$. *Verified on 8/8 instances; this is the structural reason PA-1's
involutions confine it to $a\le2$.*

---

## C. SYMBOLIC AND BRUTE-FORCE VERIFICATION

### C.1 Two independent routes per identity

| Route | Never uses | Purpose |
|---|---|---|
| **Brute force over local states** | the transfer matrix; enumerates all $(P+1)^a$ words and applies Definition 2 directly | structurally independent oracle for T3 |
| **Transfer matrix / trace** | word enumeration; exact matrix powering over $\mathbb Z[z]$ | the theorem's object |
| **sympy over $\mathbb Q[z]$, $\mathbb Q[z,t]$** | our $\mathbb Z[z]$ kernel entirely | independent symbolic confirmation |
| **Polynomial route for T1** | the cycle decomposition; builds $g,h,h^\#$ and $\operatorname{lcm}$ | independent oracle for T1 |

### C.2 Exact run figures

```
checks                    : 32 passed, 0 failed        (4.1 s)
T1                        : 126/126 agree, 14 skipped (hypothesis (H) fails)
trace vs brute force      : 24/24 (P,a) pairs
sympy symbolic            : 39/39 across 6 identity families
closed form T4            : a = 1..14 vs trace; a = 1..10 vs sympy; a = 1..6 vs frozen lists
PA-1 per-cycle objects    : W_1, W_2 for P = 1..7, exact
edge cases                : 12/12 combined
unit tests                : 185 new in test_cycle_poly.py + 6 in test_zz.py;
                            566 passed / 5 skipped / 0 failed for the whole
                            suite (375 passed coming from P1)
machine-readable outputs  : results/phase2_summary.json, results/phase2_cycle_polys.csv (51 rows)
```

### C.3 Reproducibility -- measured, not assumed

Two distinct determinism properties, each **tested** rather than asserted:

**(1) `scripts/run_phase2.py` is deterministic.** The script was run twice and
the outputs compared byte for byte: `phase2_cycle_polys.csv` is **bit-identical**,
and `phase2_summary.json` differs **only** in the `elapsed_seconds` field (4.47 s
vs 4.10 s). The mathematical path contains no RNG, so all 32 checks and every
polynomial are reproducible exactly.

**(2) The P1 test suite contains one flaky test, inherited from P1 and entirely
unrelated to P2.** Repeated runs of `tests/test_oracle.py` gave
`90 passed, 2 skipped` in 5 of 6 runs and `91 passed, 1 skipped` once. The cause
is that the `galois` library's equal-degree factorisation draws **random**
polynomials, so whether $x^{21}+1$ over $\mathbb F_{16}$ splits inside its
1000-try budget varies run to run. This **refines OPEN-30**: it is not a hard
oracle failure but a *flaky* one. Our intrinsic certificate passes in every
case, so the oracle is decorative at those two nodes and no P2 number depends on
it. P1's status is unchanged.

---

## D. PA-1 SPECIALISATION CHECK

| Object | Ours | PA-1 | Verified |
|---|---|---|---|
| length-1 per-cycle GF | $W_1^{(P)}=\sum_x z^{\min\{x,P-x\}}$ | $\sum_b\|\{b,P-b\}\|z^b$ | ✅ exact, $P=1..7$ |
| length-2 per-cycle GF | $W_2^{(P)}=\sum_{v,w}z^{\min\{v,P-w\}+\min\{w,P-v\}}$ | $\sum_{a=0}^{P}2^{1-\lfloor a/P\rfloor}(a+1)z^a$ | ✅ exact, $P=1..7$ |
| $P=1$ total count | $N=2^B\prod_{a(c)=2}(1+y^{d(c)})$ | Cor 12: $2^{s+t}\|h(\ell)\|$ from Eq. (26) | ✅ identical term-by-term |
| LCD count | $2^B$ | Rem 13(2): $2^{s+t}$ | ✅ equal since $B=s+t$ when $a\le2$ |
| scope | $k$-Galois, general $\lambda$, any $a$, any $P$ | $\lambda=\pm1$, Euclidean/Hermitian, $a\le2$ | ✅ boundary established |

**Stated plainly, per the frozen Class A rule:** within PA-1's scope our T1, T2,
T5 are equivalent reformulations of PA-1's Theorem 5 / Theorem 11 /
Corollary 12 / Remark 13(2), and our T4 collapses to PA-1's Eq. (26). **They are
not new there.** The transfer-matrix structure is claimed only as *load-bearing
for $a\ge3$*, which is exactly the frozen Class C1 wording.

---

## E. FILES CREATED / MODIFIED

**Created — implementation**

```
galois_hull_enum/core/zz.py            exact Z[z] polynomials and matrices
galois_hull_enum/core/cycle_poly.py    local weight, cycle polynomial, transfer matrix
galois_hull_enum/core/hull_dim.py      T1 verification instrument
```

**Created — tests**

```
galois_hull_enum/tests/test_zz.py
galois_hull_enum/tests/test_cycle_poly.py
```

**Created — scripts and results**

```
galois_hull_enum/scripts/run_phase2.py
galois_hull_enum/results/phase2_summary.json
galois_hull_enum/results/phase2_cycle_polys.csv
galois_hull_enum/results/phase2_tests.log
galois_hull_enum/results/phase2_run.log
```

**Modified**

```
core/zz.py            (+ mat_neg_like, used by the tests)
docs/P2_REPORT.md     (this file)
```

**Untouched:** `BLUEPRINT.md` (sha256 `2e02572a…`), `EXECUTION_PLAN.md`
(sha256 `6026536e…`), `Galios Hull.zip`, `sources/`, P0 and P1 reports, all P1
core modules and all P1 results. Verified by `git status` before commit.

---

## F. OPEN ISSUES

### F.1 Carried from P1, unchanged

| ID | Item | Status |
|---|---|---|
| **OPEN-22** | Original 385-instance suite not reconstructible (exact $\lambda$, $k$, odd $q$, $n$ absent from the frozen Blueprint) | **OPEN / NOT RECONSTRUCTIBLE** |
| **OPEN-28** | Provenance label for the newly chosen odd-$q$ inputs in `regression_suite_v1` | **retained** |
| **OPEN-29** | CP1-3/CP1-4 shape-only (sources' primitive elements in $\mathbb F_{25}$/$\mathbb F_{81}$ not recoverable from the mangled PDFs) | **OPEN** (P4) |
| **OPEN-30** | `galois` fails **flakily** on $x^{21}+1$ over $\mathbb F_{16}$/$\mathbb F_{256}$; 1-2 oracle tests skipped depending on the run (its EDF is randomised -- see section C.3(2)) | **OPEN** -- oracle limitation, ours not implicated |
| **OPEN-23/24/25/26/27** | instances.yaml, targets.yaml, SOURCE_INVENTORY.md, references.bib, untested setup_env.sh | carried |

### F.2 New in P2

| ID | Item | Impact | Owner |
|---|---|---|---|
| **OPEN-31** | The **global** enumerator $N$ (T2) has *not* been assembled or checked against a brute-force oracle over all codes. P2 verifies the per-cycle factors $W_a^{(P)}$ and the decomposition argument only. | Blocks nothing in P2; this is precisely P3's gate (CP3). | **P3** |
| **OPEN-32** | The structural reductions of Blueprint §3.5 — the Hadamard/rank-one form $z^{\min\{x,P-y\}}=\prod_{i=1}^P(1+(z-1)a_i(x)b_i(y))$ and the LDL$^\mathsf T$ form $M=LDL^\mathsf T$, $T_P=MR$, $\operatorname{tr}(T_P^a)=\operatorname{tr}((DN)^a)$ — are marked **[EXPECTED — needs verification]** in the frozen Blueprint and were **not verified in P2**. | Optional depth only; the frozen Blueprint states the paper does not depend on them. | **P2 extension / P7** |
| **OPEN-33** | The Hankel form $N[i,j]=(P+1-i-j)^+$ in §3.5(2) is likewise unverified. | as OPEN-32 | as OPEN-32 |
| **OPEN-34** | $T_6$ (self-orthogonal / dual-containing counts via $S_\ge,S_\le$) was **not** implemented or verified in P2. | P7 scope | **P7** |
| **OPEN-35** | The restricted/unrestricted mean difference quoted in Blueprint §3.6 ($P=3$: $34/16$ vs $24/12$) was **not** recomputed in P2; it belongs to T7/moments. | P6 scope | **P6** |
| **OPEN-36** | T8 (CLT) untouched — highest-risk item, deliberately not approached in P2. | P8 scope, 3-tier mitigation in Blueprint §4.6 | **P8** |

---

## G. P2 VERDICT

### G.1 Against the brief's acceptance criteria

| # | Criterion | Result |
|---|---|---|
| 1 | Every core identity has an explicit derivation | ✅ T1, T2, T3, T4, T5, Definition 5 — all derived in §A |
| 2 | Global product formula follows from the cycle decomposition | ✅ T2 derived from the bijection + T1 additivity (§A.6) |
| 3 | Trace representation is computationally verified | ✅ 24 brute-force pairs + 16 sympy pairs |
| 4 | $a\ge3$ is explicitly tested | ✅ 3 independent constructions + a real 4-cycle from P1 |
| 5 | $P=1$ reduction agrees with the frozen PA-1 specialisation | ✅ $W_1,W_2$ exact $P=1..7$; $P=1$ count identical to PA-1 Cor 12 |
| 6 | No unsupported novelty claim is introduced | ✅ PA-1's scope credited explicitly; only the $a\ge3$ load-bearing role claimed |
| 7 | All numerical results produced by actual executed tests | ✅ 32 checks + 191 tests, all logged |

### G.2 Verdict

> ## **P2 = PASS**
>
> The cycle-local enumeration framework is established. $W_a^{(P)}=\operatorname{tr}(T_P^a)$
> is derived from first principles and verified by two independent routes plus
> sympy; the $P=1$ closed form is derived and verified; the product
> decomposition follows from the cycle bijection and T1's additivity; T1 itself
> is verified on 126 checks against the polynomial route; the $a\ge3$ regime is
> constructed and verified three independent ways including on a real 4-cycle
> from the P1 core; and the PA-1 boundary is established computationally, with
> PA-1's scope credited rather than overlapped.
>
> **No frozen formula was altered.** Four first-run failures and two unit-test
> failures were diagnosed as defects in my verification code and are documented
> in §A.14 rather than silently repaired.
>
> **Scope honoured.** No enumerator over all codes (P3), no moments (P6), no
> limit law (P8), no manuscript text. The frozen Blueprint and execution plan
> are byte-identical to `076f910` / `96c7885`.
>
> **Ready for P3**, which must assemble $N$ from the now-verified $W$ factors
> and validate it against Algorithm E3's brute force (OPEN-31).

**STOPPING after P2 as instructed.**
