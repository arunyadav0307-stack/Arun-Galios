r"""The frozen hull-dimension identity T1, and the polynomial route to it.

This module exists ONLY to verify T1 (BLUEPRINT.md section 3.3, proof in
section 4.1):

    dim hull_k(C) = sum_{mathfrak c} ord_{j(mathfrak c)}(q) * w_{a(mathfrak c)}^{(P)}(u^{(mathfrak c)})

by computing the LEFT-HAND SIDE from polynomials -- g, h, h^{#} and
deg lcm(g, h^{#}) -- exactly as the Blueprint's validation principle
(section 4.7) requires. That route never uses the cycle decomposition, so
agreement with the RHS is a genuine check.

SCOPE NOTE. This is the T1 verification instrument. It is NOT Algorithm E2
(the enumerator N, which multiplies the W factors and histograms over ALL
codes) -- that is P3's work and is deliberately not implemented here.

INDEPENDENT IMPLEMENTATION. Nothing is imported from or adapted from the
throwaway blueprinting prototype, which lives outside this repository.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

from . import poly as P
from .cycles import hash_map
from .factor import factor_xn_minus_lambda
from .cycle_poly import local_weight

__all__ = [
    "parity_check", "hash_map_poly", "polylcm", "hull_dimension_poly",
    "build_generator_from_cycle_words", "cycle_words_from_instance",
    "verify_T1_on_instance",
]


def parity_check(F, n: int, lam: int, g: Sequence[int]) -> List[int]:
    """h(x) = (x^n - lambda) / g(x) -- the parity-check polynomial of g."""
    xn_minus_lam = [F.neg(lam)] + [0] * (n - 1) + [1]
    q, r = P.pdivmod(F, xn_minus_lam, g)
    if r:
        raise ValueError("g does not divide x^n - lambda; g is not a valid "
                         "generator polynomial")
    return q


def hash_map_poly(F, h: Sequence[int], j: int) -> List[int]:
    """h^{#} = h_0^{-p^j} x^{deg h} h^{p^j}(1/x)  -- the frozen map, section 3.1."""
    return hash_map(F, h, j)


def polylcm(F, a: Sequence[int], b: Sequence[int]) -> List[int]:
    """lcm(a, b) = a b / gcd(a, b), monic."""
    g = P.pgcd(F, a, b)
    if not g:
        return []
    return P.pmonic(F, P.pmul(F, a, P.pdivmod(F, b, g)[0]))


def hull_dimension_poly(F, n: int, lam: int, k: int, g: Sequence[int]) -> int:
    """dim hull_k(C) = n - deg lcm(g, h^{#})  -- the polynomial route.

    P1's Remark 2 / BLUEPRINT section 2.0.4: this identity is what makes
    hull_k(C) lambda-constacyclic, and it holds under hypothesis (H)
    lambda^{1+p^{e-k}} = 1.
    """
    j = F.e - k
    h = parity_check(F, n, lam, g)
    h_sharp = hash_map_poly(F, h, j)
    lcm = polylcm(F, g, h_sharp)
    return n - P.deg(lcm)


def build_generator_from_cycle_words(F, cycle_data_obj, words) -> List[int]:
    """Build g = prod over cycles of prod_m f_m^{u_m} from explicit words."""
    g = [1]
    for cyc, u in zip(cycle_data_obj.cycles, words):
        for idx, exp in zip(cyc["orbit"], u):
            if exp:
                f = cycle_data_obj.factors[idx][0] if isinstance(
                    cycle_data_obj.factors[idx], tuple) else cycle_data_obj.factors[idx]
                for _ in range(exp):
                    g = P.pmul(F, g, f)
    return P.pmonic(F, g)


def cycle_words_from_instance(F, n: int, lam: int, k: int):
    """Return (cycle_data, list_of_words) where the words are all-zero -- the
    canonical starting point, plus the cycle_data itself."""
    from .cycles import e1_cycle_data
    cd = e1_cycle_data(F, n, lam, k)
    words = [tuple([0] * c["a"]) for c in cd.cycles]
    return cd, words


def verify_T1_on_instance(F, n: int, lam: int, k: int, word_sets=None,
                          verbose: bool = False):
    """Verify T1 on one instance for a list of candidate word assignments.

    For each candidate, compute
        LHS = n - deg lcm(g, h^{#})      (polynomial route)
        RHS = sum_c d_c w_{a_c}(u_c)     (local-weight route, Definition 2)
    and require equality. Returns a list of result dicts.

    `word_sets` defaults to a deterministic spread of words including the
    all-zero word, the all-P word, and a few mixed words.
    """
    from .cycles import e1_cycle_data
    cd = e1_cycle_data(F, n, lam, k)
    if not cd.permutation_ok:
        return [{"instance": f"F_{F.q} n={n} k={k}",
                 "skipped": "hypothesis (H) fails; # does not permute the factors"}]

    if word_sets is None:
        word_sets = []
        # all-zero
        word_sets.append([tuple([0] * c["a"]) for c in cd.cycles])
        # all-P
        word_sets.append([tuple([cd.P] * c["a"]) for c in cd.cycles])
        # alternating 0/P where possible
        alt = []
        for c in cd.cycles:
            alt.append(tuple([cd.P if (i % 2) else 0 for i in range(c["a"])]))
        word_sets.append(alt)
        # a deterministic pseudo-spread
        import random
        rng = random.Random(n * 1000 + k)
        for _ in range(4):
            ws = []
            for c in cd.cycles:
                ws.append(tuple(rng.randrange(cd.P + 1) for _ in range(c["a"])))
            word_sets.append(ws)

    results = []
    for words in word_sets:
        g = build_generator_from_cycle_words(F, cd, words)
        lhs = hull_dimension_poly(F, n, lam, k, g)
        rhs_pairs = [(u, c["d"]) for u, c in zip(words, cd.cycles)]
        rhs = sum(d * local_weight(u, cd.P) for u, d in rhs_pairs)
        ok = (lhs == rhs)
        results.append({
            "instance": f"F_{F.q} n={n} k={k} lambda={lam}",
            "words": [list(w) for w in words],
            "lhs_polynomial_route": lhs,
            "rhs_local_weight_route": rhs,
            "agree": bool(ok),
        })
        if verbose:
            print(f"    words={[list(w) for w in words]}  lhs={lhs} rhs={rhs} "
                  f"{'OK' if ok else 'MISMATCH'}")
    return results
