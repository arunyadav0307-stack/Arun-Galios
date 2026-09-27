r"""The #-map and Algorithm E1 cycle detection.

INDEPENDENT IMPLEMENTATION. Written from scratch for P1; nothing is imported,
copied or adapted from the throwaway blueprinting prototype, which lives
outside this repository.

The frozen construction (BLUEPRINT.md section 3.1, verbatim):

    For f(x) = sum_{i=0}^{m} f_i x^i with f_0 != 0,

        f^{#}(x) = sum_{i=0}^{m} f_0^{-p^{e-k}} f_i^{p^{e-k}} x^{m-i}

    equivalently   f^{#} = f_0^{-sigma^{e-k}} x^{m} f^{sigma^{e-k}}(1/x),

where sigma = (.)^p is the Frobenius and j := e - k.

Two immediate consequences, both used as unit tests below:

  * deg f^{#} = deg f, and the coefficient of x^{m} in f^{#} is
    f_0^{-p^j} f_0^{p^j} = 1, so f^{#} is ALREADY MONIC (the "monic
    normalisation" in E1 step 3 is therefore a no-op; we still apply it for
    safety).
  * f^{#} is not an involution in general: the cycle length is the order of
    the map (alpha -> alpha^{-p^{e-k}}) and for k not in {0, e/2} it can
    exceed 2. That is the whole point of contribution class C1/C2.

Algorithm E1 (BLUEPRINT.md section 5.1), steps 1-4:

  1. n' = n / p^nu, P = p^nu, mu with mu^{P} = lambda.
  2. Factor x^{n'} - mu (squarefree) into distinct monic irreducibles.
  3. For each monic irreducible f, iterate f -> monic(f^{#}) to close the
     cycle; record (a, d) with d = deg f = ord_j(q).
  4. Output the multiset of cycles; B = |C|; sum_c a(c) = number of distinct
     monic irreducible factors of x^n - lambda.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from .poly import deg, format_poly, frobenius_coeffs, pmonic, trim
from .factor import factor_xn_minus_lambda

__all__ = ["hash_map", "hash_map_monic", "cycle_data", "e1_cycle_data", "CycleData"]


def hash_map(F, f: Sequence[int], j: int) -> List[int]:
    r"""The frozen f^{#} map, with j = e - k.

    f^{#}(x) = sum_{i=0}^{m} f_0^{-p^j} f_i^{p^j} x^{m-i},
    where m = deg f and f_0 != 0.
    """
    f = trim(f)
    m = deg(f)
    if m < 0:
        raise ValueError("f^{#} is not defined for the zero polynomial")
    f0 = f[0]
    if f0 == 0:
        raise ValueError(r"f^{#} requires f(0) != 0 (irreducible factors of "
                         r"x^n - lambda always satisfy this)")

    # c = f_0^{-p^j} = sigma^j(f_0^{-1}) = (sigma^j(f_0))^{-1}
    c = F.inv(F.frobenius(f0, j))

    # apply sigma^j to every coefficient, then reverse the order
    fj = frobenius_coeffs(F, f, j)          # [f_i^{p^j}]_{i=0..m}
    out = [F.mul(c, coeff) for coeff in reversed(fj)]
    return trim(out)


def hash_map_monic(F, f: Sequence[int], j: int) -> List[int]:
    """monic(f^{#}) -- the exact iterate used by E1 step 3.

    f^{#} is already monic whenever f(0) != 0, so this is a no-op in exact
    arithmetic; it is applied anyway to mirror E1 literally and to guard
    against any representation slip.
    """
    return pmonic(F, hash_map(F, f, j))


class CycleData:
    """Container for the output of Algorithm E1."""

    def __init__(self, q, p, e, k, j, n, n_prime, nu, P, lam, mu,
                 factors, cycles, permutation_ok):
        self.q = q
        self.p = p
        self.e = e
        self.k = k
        self.j = j
        self.n = n
        self.n_prime = n_prime
        self.nu = nu
        self.P = P
        self.lam = lam
        self.mu = mu
        self.factors = factors                 # list of (poly, multiplicity)
        self.cycles = cycles                   # list of dicts, see cycle_data
        self.permutation_ok = permutation_ok   # does # permute the factors?

    # ---- quantities named in BLUEPRINT.md section 3.1 / 5.1 ----
    @property
    def distinct_irreducible_factors(self) -> int:
        return len(self.factors)

    @property
    def B(self) -> int:
        """|C| -- the number of cycles."""
        return len(self.cycles)

    @property
    def sum_a(self) -> int:
        """sum over cycles of a(c); must equal the number of distinct
        monic irreducible factors of x^n - lambda."""
        return sum(c["a"] for c in self.cycles)

    @property
    def predicted_code_count(self) -> int:
        r"""|C_{n,q,lambda}| = (P+1)^{sum_c a(c)} -- E1's stated output.

        Reported but NOT validated here: validating it against an enumeration
        is P3's job (CP3).
        """
        return (self.P + 1) ** self.sum_a

    @property
    def cycle_lengths(self) -> List[int]:
        return sorted(c["a"] for c in self.cycles)

    @property
    def cycle_shapes(self) -> List[Tuple[int, int]]:
        """[(a, d), ...] sorted -- the (length, degree) pairs of E1 step 3."""
        return sorted((c["a"], c["d"]) for c in self.cycles)

    def summary(self) -> str:
        return (f"F_{self.q} n={self.n} k={self.k} j={self.j} "
                f"n'={self.n_prime} nu={self.nu} P={self.P} | "
                f"factors={self.distinct_irreducible_factors} B={self.B} "
                f"sum_a={self.sum_a} |C|~={self.predicted_code_count} "
                f"shapes={self.cycle_shapes}")


def cycle_data(F, factors: Sequence[Sequence[int]], j: int):
    """Given the distinct monic irreducible factors of x^{n'} - mu and
    j = e - k, return (cycles, permutation_ok).

    permutation_ok is False if some f^{#} is not again a factor -- which is
    what happens when the standing hypothesis (H) lambda^{1+p^{e-k}} = 1 is
    violated. (H) is recorded in BLUEPRINT.md section 2.0.4 and is what makes
    hull_k(C) lambda-constacyclic.
    """
    flist = [pmonic(F, f) for f in factors]
    index: Dict[Tuple[int, ...], int] = {}
    for i, f in enumerate(flist):
        index[tuple(f)] = i

    # image of each factor under #
    perm: List[int] = []
    permutation_ok = True
    for f in flist:
        g = hash_map_monic(F, f, j)
        key = tuple(g)
        if key not in index:
            permutation_ok = False
            perm.append(-1)
        else:
            perm.append(index[key])

    cycles: List[dict] = []
    if permutation_ok:
        seen = [False] * len(flist)
        for start in range(len(flist)):
            if seen[start]:
                continue
            orbit = []
            cur = start
            while not seen[cur]:
                seen[cur] = True
                orbit.append(cur)
                cur = perm[cur]
            a = len(orbit)
            d = deg(flist[orbit[0]])
            # every element of a cycle has the same degree (# preserves degree)
            if any(deg(flist[i]) != d for i in orbit):
                raise AssertionError("a cycle contains factors of differing degree")
            cycles.append({
                "a": a,
                "d": d,
                "orbit": orbit,
                "representative": flist[orbit[0]],
                "poly": format_poly(F, flist[orbit[0]]),
            })
        cycles.sort(key=lambda c: (c["a"], c["d"]))
    return cycles, permutation_ok


def e1_cycle_data(F, n: int, lam: int, k: int) -> CycleData:
    """Full Algorithm E1: from (n, q, lambda, k) to the cycle multiset."""
    if not (0 <= k <= F.e - 1):
        raise ValueError(f"k must lie in [0, {F.e - 1}]")
    j = F.e - k
    n_prime, nu, P, mu, factors = factor_xn_minus_lambda(F, n, lam)
    distinct = [f for (f, _m) in factors]
    cycles, perm_ok = cycle_data(F, distinct, j)
    return CycleData(
        q=F.q, p=F.p, e=F.e, k=k, j=j, n=n, n_prime=n_prime, nu=nu, P=P,
        lam=lam, mu=mu, factors=factors, cycles=cycles, permutation_ok=perm_ok,
    )
