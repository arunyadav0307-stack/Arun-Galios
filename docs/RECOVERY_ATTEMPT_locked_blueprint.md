# RECOVERY ATTEMPT — `locked_blueprint.docx` and the design reference PDF

**Outcome: NOT RECOVERABLE. Neither file exists in this workspace, in git
history, on the remote, or anywhere on the filesystem.**

**Nothing was reconstructed, recreated, substituted or fabricated. Neither file
was modified (they do not exist to modify). No audit was performed.**

---

## 1. `locked_blueprint.docx` — ABSENT

Exhaustive search, all negative:

| # | Method | Result |
|---|---|---|
| 1 | `find / -iname "*.docx"` | **0 results** (whole filesystem, excluding `/proc`, `/sys`) |
| 2 | `find / -iname "*locked*"` | only `docs/AUDIT_locked_blueprint.md` (my own audit report, not the audited document) |
| 3 | `git rev-list --objects --all` filtered to `.docx/.pdf/.odt/.doc` | **0 results** — never committed to any reachable object |
| 4 | `git log --all --diff-filter=A --name-only` | **0** such paths ever added, on any branch or tag |
| 5 | `git ls-remote origin` | remote has only `main`, `arena/01a0e2d5-arun-galios`, tags `p0`–`p3`. No document blob. |
| 6 | `git stash list` | empty |
| 7 | `git fsck --dangling --lost-found` | **0** dangling/unreachable blobs |
| 8 | `git reflog --all` | only clone/checkout entries; no lost commits |

**Conclusion:** no `locked_blueprint.docx` was ever created in this repository,
committed, pushed, or left on disk.

## 2. Design reference PDF — ABSENT

| # | Method | Result |
|---|---|---|
| 1 | `find / -iname "*.pdf"` | **0 results** on the filesystem outside `/proc`, `/sys` |
| 2 | git history / remote / dangling objects | **0** PDFs ever committed or reachable |
| 3 | `Galios Hull.zip` contents | 3 PDFs — see §3 — all **scientific source papers**, not design references |

**Conclusion:** no design/mechanics reference PDF exists. The only PDFs in the
workspace are the three research source papers inside `Galios Hull.zip`, which
are *scientific* content already extracted to `sources/*.txt`. They are not a
design reference and must not be treated as one.

## 3. What PDFs *do* exist (identified, for the record)

Inside `Galios Hull.zip` (the only primary research/work source):

| Path in ZIP | Bytes | SHA-256 |
|---|---|---|
| `Galios Hull/Galios hulls of constacyclic codes over finite fields.pdf` | 2,597,643 | `aef8c9c8c8aedca3c17cbab8b466bd48dae5f7b6c0ebf0e420e7389684fdb765` |
| `Galios Hull/Average_dimensions_of_Galois_hulls_of_constacyclic.pdf` | 664,819 | `57561ee72adb6c7b4c5157698ea0b1174ece32b3eb7bc74b59d4f49fd7b6beed` |
| `Galios Hull/Galios Hull over Affin Algebra.pdf` | 291,546 | `0db90bfbeda69ad1af231b81d92f24be3589ed15ac2e0ade6605ba7ad0ddc039` |

Extracted text counterparts (already in the repo):

| File | SHA-256 |
|---|---|
| `sources/P1_galois_hulls_constacyclic.txt` | `53f2afd151b4f0418383e5a7b2c17ae3a3ac8b9d500deebf82293062361a034a` |
| `sources/P2_affine_algebra.txt` | `f3d6534f4ed3a8e035ac79ba0692e993c0395a716374a133bd5e7bb0475afa7e` |
| `sources/P3_average_dimensions.txt` | `f231f9da69ed9457ebbe171aa0f1d30702973507fd7145c94fccbc9d3b76079f` |

## 4. Provenance of the premise

The session record for this workspace covers P0, P1, P2, P3 and the blocked
audit. **There is no design-conversion step in that record, and no `.docx` was
ever produced.** Every artefact written in this session was a `.md`, `.py`,
`.json`, `.csv`, `.log` or `.yaml` file. The sandbox has been re-cloned at least
twice during this work (local history reset to `00fefba` both times), so a file
created *outside* the persisted workspace root would not have survived — but
searches 1–8 above show no trace of such a file either, on disk or in git.

It is therefore not possible to confirm:

1. the file path or SHA-256 of `locked_blueprint.docx` — **no such file exists**;
2. that it is the document generated from a previous design-conversion step —
   **no such step exists in the record**;
3. the identity or SHA-256 of the design reference PDF — **no such PDF exists**.

---

## Standing correction applied

The frozen `BLUEPRINT.md` (sha256 `2e02572a0ef970a93dcf242c725f290938d9b180787beaccd8acedf6a609214f`,
commit `076f910`) defines exactly **one** standing hypothesis:

> **(H)**  `λ^(1+p^(e-k)) = 1`  (§2.0.4)

The restrictions `r ∣ (1+p^{e-k})` and `gcd(n',r)=1` are **limitations /
restrictions**, recorded as **Limitation L2** in §8.2 and as standing data in
§3.1. They are **not** hypotheses and are **not** labelled H1/H2 anywhere in the
frozen source. No H1/H2 is treated as frozen in any report produced here.

## Next step

Both files must be supplied by the user (or their creation step re-run in a
session whose artefacts persist). Until then the design-fidelity audit remains
**BLOCKED — NOT ASSESSABLE**, and the checklist in
`docs/AUDIT_locked_blueprint.md` §H.1 stands ready to be applied.
