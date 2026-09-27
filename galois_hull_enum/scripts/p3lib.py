r"""Shared helpers for the P3 (CP3) validation driver.

Everything here is plumbing: turning the frozen Blueprint's section 5.5
experiment grid into concrete (q, n, lambda, k) tuples, applying the standing
hypotheses as filters, and running one instance through all three routes.
"""

from __future__ import annotations

import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import poly as P
from core.cycles import CycleData, e1_cycle_data
from core.factor import factor_xn_minus_lambda, v_p
from core.field import FiniteField

# The frozen section 5.5 grid.  lambda is given as an ORDER (r) wherever the
# Blueprint writes it as a power of a generator, because the element itself
# depends on the irreducible polynomial chosen for F_q -- which the frozen
# Blueprint does not fix.  Using r is exact and loses nothing: the standing
# hypotheses and the whole cycle structure depend on lambda only through
# r = ord(lambda) together with the actual factorisation, and the factorisation
# is recomputed from scratch for the chosen representative.
GRID: List[Dict[str, object]] = [
    {"id": "X1",  "q": 9,  "k": 1, "n": 3,  "r": 2,   "src": "P1 Ex. 4"},
    {"id": "X2",  "q": 25, "k": 1, "n": 7,  "r": 3,   "src": "P1 Ex. 5"},
    {"id": "X3",  "q": 27, "k": 2, "n": 15, "r": 2,   "src": "P1 Ex. 6"},
    {"id": "X4",  "q": 81, "k": 2, "n": 21, "r": 10,  "src": "P3 Ex. 3.6 / Table 1"},
    {"id": "X5",  "q": 81, "k": 3, "n": 5,  "r": 1,   "src": "new (l=4 regime)"},
    {"id": "X6",  "q": 81, "k": 3, "n": 17, "r": 1,   "src": "new (l=4 regime)"},
    {"id": "X7",  "q": 81, "k": 3, "n": 25, "r": 1,   "src": "new (l=4 regime)"},
    {"id": "X8a", "q": 16, "k": 3, "n": 5,  "r": 1,   "src": "new (even q, l=4)"},
    {"id": "X8b", "q": 16, "k": 3, "n": 15, "r": 3,   "src": "new (even q, l=4)"},
    {"id": "X8c", "q": 16, "k": 3, "n": 17, "r": 1,   "src": "new (even q, l=4)"},
    {"id": "X9",  "q": 64, "k": 3, "n": 20, "r": 9,   "src": "P3 Ex. 3.1"},
    {"id": "X10a", "q": 25, "k": 1, "n": 65, "r": 6,  "src": "P3 Ex. 4.13 / Table 2"},
    {"id": "X10b", "q": 25, "k": 1, "n": 65, "r": 3,  "src": "P3 Ex. 4.13 / Table 2"},
]

# X11-X13 are RING instances (the affine algebra A and R_{m,q}); the frozen
# Blueprint assigns them to CP6, which is outside this phase.
RING_ROWS = ["X11", "X12", "X13"]

# the X14 sweep fields (frozen, section 5.5)
X14_FIELDS = [4, 8, 9, 16, 25, 27, 32, 49, 64, 81, 121, 125, 128, 256]

MAX_CODES = 10 ** 5          # frozen E3 feasibility bound (pure Python)


# --------------------------------------------------------------------------
def element_of_order(F: FiniteField, r: int) -> Optional[int]:
    """A deterministic element of multiplicative order exactly r in F_q^*.

    Returns None when r does not divide q-1 (no such element exists).
    """
    if r == 1:
        return F.one
    if (F.q - 1) % r:
        return None
    # walk the nonzero elements in the kernel's deterministic order
    for a in F.nonzero_elements():
        if F.multiplicative_order(a) == r:
            return a
    return None


def hypotheses_ok(F: FiniteField, n: int, lam: int, k: int) -> Tuple[bool, str]:
    """The standing hypotheses of BLUEPRINT section 2.0.4 / 3.1.

    Checks r | (1 + p^{e-k}) and gcd(n', r) = 1, plus the extra identity
    mu^{p^nu} = lambda that E1 itself enforces.
    """
    r = F.multiplicative_order(lam)
    if (1 + F.p ** (F.e - k)) % r:
        return False, f"r={r} does not divide 1+p^(e-k)={1 + F.p ** (F.e - k)}"
    nu = v_p(n, F.p)
    n_prime = n // (F.p ** nu)
    if math.gcd(n_prime, r) != 1:
        return False, f"gcd(n'={n_prime}, r={r}) != 1"
    return True, ""


# --------------------------------------------------------------------------
class Instance:
    """One (q, n, lambda, k) together with everything computed from it."""

    def __init__(self, F: FiniteField, n: int, lam: int, k: int,
                 label: str = "", source: str = ""):
        self.F = F
        self.n = n
        self.lam = lam
        self.k = k
        self.label = label
        self.source = source
        self.cd: Optional[CycleData] = None
        self.factors = None
        self.skip_reason: Optional[str] = None

    # ---------------- construction / filtering ----------------
    def prepare(self) -> bool:
        """Run E1 and the standing hypotheses. False => not a valid instance."""
        ok, why = hypotheses_ok(self.F, self.n, self.lam, self.k)
        if not ok:
            self.skip_reason = "hypotheses: " + why
            return False
        try:
            self.cd = e1_cycle_data(self.F, self.n, self.lam, self.k)
        except Exception as exc:                      # noqa: BLE001
            self.skip_reason = f"E1 failed: {exc}"
            return False
        if not self.cd.permutation_ok:
            self.skip_reason = "hypothesis (H) fails; # does not permute the factors"
            return False
        self.factors = factor_xn_minus_lambda(self.F, self.n, self.lam)[4]
        return True

    # ---------------- descriptive ----------------
    @property
    def key(self) -> str:
        return f"F_{self.F.q}|n={self.n}|k={self.k}|lam={self.lam}"

    @property
    def tag(self) -> str:
        return (f"{self.label or ''}F_{self.F.q} n={self.n} k={self.k} "
                f"lam={self.lam}").strip()

    def summary_dict(self) -> Dict[str, object]:
        cd = self.cd
        return {
            "label": self.label,
            "source": self.source,
            "q": self.F.q,
            "p": self.F.p,
            "e": self.F.e,
            "n": self.n,
            "k": self.k,
            "lam": self.lam,
            "r": self.F.multiplicative_order(self.lam),
            "P": cd.P,
            "nu": cd.nu,
            "n_prime": cd.n_prime,
            "j": cd.j,
            "B": cd.B,
            "sum_a": cd.sum_a,
            "n_factors": len(self.factors),
            "shapes": [list(s) for s in cd.cycle_shapes],
            "cycle_lengths": cd.cycle_lengths,
            "n_codes_predicted": cd.predicted_code_count,
            "permutation_ok": cd.permutation_ok,
        }

    # ---------------- the three routes ----------------
    def run_all(self, crosscheck: bool = True):
        """Return (e2_dist, oracle_cycle_dist, oracle_flat_dist, timings)."""
        from enumeration.enumerator import enumerator_distribution
        from oracle.brute import (oracle_distribution_cycle,
                                  oracle_distribution_flat)

        t0 = time.time()
        e2 = enumerator_distribution(self.cd.cycle_shapes, self.cd.P,
                                     crosscheck=crosscheck)
        t1 = time.time()
        oc = oracle_distribution_cycle(self.F, self.cd, self.n, self.lam, self.k)
        t2 = time.time()
        of, st = oracle_distribution_flat(self.F, self.n, self.lam, self.k,
                                          self.factors)
        t3 = time.time()
        return e2, oc, of, {"e2": t1 - t0, "oracle_cycle": t2 - t1,
                            "oracle_flat": t3 - t2, "oracle_flat_stats": st}


def grid_instances(rows: Sequence[Dict[str, object]]) -> List[Instance]:
    """Materialise the frozen grid rows as Instances."""
    out: List[Instance] = []
    for row in rows:
        q = int(row["q"])                                  # type: ignore[arg-type]
        p, e = _pe(q)
        F = FiniteField(p, e)
        lam = element_of_order(F, int(row["r"]))            # type: ignore[arg-type]
        if lam is None:
            continue
        out.append(Instance(F, int(row["n"]), lam, int(row["k"]),      # type: ignore[arg-type]
                            label=str(row["id"]) + " ", source=str(row["src"])))
    return out


def _pe(q: int) -> Tuple[int, int]:
    for p in (2, 3, 5, 7, 11, 13):
        e = 0
        m = q
        while m % p == 0:
            m //= p
            e += 1
        if m == 1 and e:
            return p, e
    raise ValueError(f"q={q} is not a prime power")


def sweep_instances(fields: Sequence[int] = X14_FIELDS, n_max: int = 60,
                    max_codes: int = MAX_CODES, lam_per_order: int = 1,
                    k_filter: Optional[Sequence[int]] = None
                    ) -> Tuple[List[Instance], List[Dict[str, object]]]:
    """The X14 sweep: every (q, n, lambda, k) the frozen grid specifies,
    restricted to instances the E3 oracle can actually certify.

    lambda is sampled by multiplicative order (`lam_per_order` representatives
    per order) rather than exhaustively over F_q^*, which for q = 256 would be
    255 elements x 60 values of n x e values of k. The cycle structure depends
    on lambda through r = ord(lambda) and the concrete factorisation, and every
    candidate is still factorised and certified from scratch, so nothing is
    assumed.

    Returns (instances, rejected) where `rejected` records why each candidate
    was dropped -- so the sweep's coverage is auditable rather than implied.
    """
    instances: List[Instance] = []
    rejected: List[Dict[str, object]] = []
    for q in fields:
        p, e = _pe(q)
        F = FiniteField(p, e)
        orders = sorted({F.multiplicative_order(a) for a in F.nonzero_elements()})
        for n in range(1, n_max + 1):
            for k in (k_filter if k_filter is not None else range(e)):
                for r in orders:
                    if (1 + p ** (e - k)) % r:
                        continue
                    nu = v_p(n, p)
                    n_prime = n // (p ** nu)
                    if math.gcd(n_prime, r) != 1:
                        continue
                    reps = [a for a in F.nonzero_elements()
                            if F.multiplicative_order(a) == r][:lam_per_order]
                    for lam in reps:
                        inst = Instance(F, n, lam, k, source="X14 sweep")
                        if not inst.prepare():
                            rejected.append({"q": q, "n": n, "k": k, "lam": lam,
                                             "r": r,
                                             "reason": inst.skip_reason})
                            continue
                        if inst.cd.predicted_code_count > max_codes:
                            rejected.append({
                                "q": q, "n": n, "k": k, "lam": lam, "r": r,
                                "reason": (f"|C|={inst.cd.predicted_code_count} "
                                           f"> {max_codes}; E3 infeasible "
                                           f"(enumerator only, flagged)")})
                            continue
                        instances.append(inst)
    return instances, rejected


def peak_rss_mb() -> float:
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
