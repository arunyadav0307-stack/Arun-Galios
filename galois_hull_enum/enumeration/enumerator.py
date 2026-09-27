r"""Algorithm E2 -- the global hull-dimension enumerator.

BLUEPRINT.md section 5.2, verbatim:

    dist <- {0 : 1};
    for c in C: dist <- dist * (W_{a(c)}^{(P)} with z <- y^{d(c)}),

where W_a^{(P)} is obtained either by (i) brute force over [0,P]^a when
(P+1)^a <= 10^6, or (ii) tr(T_P(z)^a) by binary powering in Z[z].
Polynomial multiplication in Z[y].

DIRECTORY-NAME DEVIATION (recorded, not silent).
BLUEPRINT section 6.1 lays this out as `enum/`. That name is NOT usable in
this repository: Python 3.11's `typing` module (imported by every core module)
does `import enum` internally, and a top-level package named `enum` on
sys.path shadows the standard-library one, which breaks `import typing` with
`AttributeError: module 'enum' has no attribute 'global_enum'`. The package is
therefore named `enumeration/`; the module contents, names and behaviour are
exactly what the Blueprint specifies. This is a packaging constraint, not a
change to the mathematics.

WHAT THIS MODULE IS
    The product of the per-cycle factors W_{a(c)}^{(P)} over the cycle
    multiset of Algorithm E1. Every factor comes from core/cycle_poly.py,
    which P2 derived and verified (Definition 3 = Definition 2 brute force =
    tr(T_P^a) = T4 closed form at P=1).

WHAT THIS MODULE IS NOT
    It does not compute any hull dimension. It never builds a polynomial g,
    never calls the #-map, and never touches a finite field. The only inputs
    it needs from the outside world are the cycle SHAPES (a, d) -- integers --
    plus P. That is deliberate: it is what makes the comparison against
    oracle/brute.py meaningful.

INDEPENDENCE FROM THE ORACLE
    E2 imports core/cycle_poly (verified in P2) and nothing else from the
    problem domain. oracle/brute.py imports core/poly, core/factor and an
    INDEPENDENT #-map implementation. The only code the two share is the
    exact-integer polynomial arithmetic primitives (core/zz.py, core/poly.py)
    and the finite-field kernel (core/field.py), neither of which encodes any
    part of the hull-dimension theory. See docs/P3_REPORT.md section B.4 for
    the full shared-dependency analysis.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from core import zz
from core.cycle_poly import (
    cycle_poly_bruteforce, cycle_poly_trace, degree_bound, dict_to_poly,
)

__all__ = [
    "substitute_z_power_y", "zy_mul", "zy_from_Zz",
    "cycle_factor", "cycle_factor_bruteforce",
    "enumerator_poly", "enumerator_distribution",
    "structural_checks", "E2Result", "run_e2",
]

# Z[y] is a tuple of ints, low degree first -- the same convention as Z[z].


# --------------------------------------------------------------------------
# substitution z <- y^d   (BLUEPRINT section 5.2)
# --------------------------------------------------------------------------
def substitute_z_power_y(W: Sequence[int], d: int) -> Tuple[int, ...]:
    """The Z[z] polynomial W evaluated at z = y^d, as a Z[y] polynomial.

    If W = sum_i c_i z^i then W(y^d) = sum_i c_i y^{i d}.
    """
    if d < 0:
        raise ValueError("d must be nonnegative")
    if not W:
        return ()
    n = zz.zdeg(W)
    out = [0] * (n * d + 1)
    for i, c in enumerate(W):
        if c:
            out[i * d] += c
    return zz.ztrim(out)


def zy_mul(a: Sequence[int], b: Sequence[int]) -> Tuple[int, ...]:
    """Multiplication in Z[y]."""
    return zz.zmul(a, b)


def zy_from_Zz(W: Sequence[int]) -> Tuple[int, ...]:
    """View a Z[z] polynomial as a Z[y] polynomial (identity on coefficients)."""
    return zz.ztrim(W)


# --------------------------------------------------------------------------
# the per-cycle factor
# --------------------------------------------------------------------------
def cycle_factor(a: int, P: int, d: int, crosscheck: bool = True) -> Tuple[int, ...]:
    """W_a^{(P)}(y^d) -- the contribution of one cycle of length a, degree d.

    Uses the TRACE route (P2's T3), which is polynomial time in a. When
    (P+1)^a <= 10^6 the brute-force route is also evaluated and the two are
    required to agree -- this is the Blueprint's own E2 clause (i)/(ii), and
    it is free insurance that the factor entering the global product is the
    object P2 verified.
    """
    if a < 1:
        raise ValueError("a >= 1 required")
    if P < 1:
        raise ValueError("P >= 1 required")
    if d < 1:
        raise ValueError("d >= 1 required (d = ord_j(q) >= 1)")

    W = cycle_poly_trace(P, a)                      # tr(T_P^a), Z[z]
    if crosscheck and (P + 1) ** a <= 10 ** 6:
        # NOTE: cycle_poly_bruteforce returns a DICT {weight: count}, not a
        # coefficient tuple; dict_to_poly converts it to Z[z].
        Wb = dict_to_poly(cycle_poly_bruteforce(P, a))
        if Wb != W:
            raise AssertionError(
                f"E2 internal inconsistency: trace != brute force for "
                f"(P,a)=({P},{a}): {zz.zformat(W)} vs {zz.zformat(Wb)}")
    return substitute_z_power_y(W, d)


def cycle_factor_bruteforce(a: int, P: int, d: int) -> Tuple[int, ...]:
    """The same factor, but by brute-force enumeration only (no matrix)."""
    return substitute_z_power_y(dict_to_poly(cycle_poly_bruteforce(P, a)), d)


# --------------------------------------------------------------------------
# the global enumerator
# --------------------------------------------------------------------------
def enumerator_poly(shapes: Sequence[Tuple[int, int]], P: int,
                    crosscheck: bool = True) -> Tuple[int, ...]:
    """N(y) = prod_c W_{a(c)}^{(P)}(y^{d(c)}) as a Z[y] coefficient tuple.

    `shapes` is the multiset of (a, d) cycle shapes from Algorithm E1.
    """
    dist: Tuple[int, ...] = (1,)
    for (a, d) in shapes:
        dist = zy_mul(dist, cycle_factor(a, P, d, crosscheck=crosscheck))
    return zz.ztrim(dist)


def enumerator_distribution(shapes: Sequence[Tuple[int, int]], P: int,
                            crosscheck: bool = True) -> Dict[int, int]:
    """N(y) as a histogram {hull dimension : number of codes}."""
    poly = enumerator_poly(shapes, P, crosscheck=crosscheck)
    if not poly:
        return {}
    return {i: c for i, c in enumerate(poly) if c}


# --------------------------------------------------------------------------
# structural checks that must hold for ANY instance (BLUEPRINT 6.2 CP3)
# --------------------------------------------------------------------------
def structural_checks(shapes: Sequence[Tuple[int, int]], P: int,
                      predicted_code_count: int, B: int) -> Dict[str, object]:
    """The three CP3 structural predicates, evaluated on the E2 output.

    1. |C| = (P+1)^{sum a}   -- total number of codes, N(1).
    2. #{dim = 0} = 2^B      -- the LCD count (T5).
    3. all coefficients nonnegative integers.
    """
    poly = enumerator_poly(shapes, P)
    total = zz.zval(poly, 1)
    lcd = poly[0] if poly else 0

    sum_a = sum(a for (a, _d) in shapes)
    return {
        "sum_a": sum_a,
        "B": B,
        "P": P,
        "n_codes_E2_N1": total,
        "n_codes_predicted": predicted_code_count,
        "n_codes_match": total == predicted_code_count,
        "lcd_count_E2": lcd,
        "lcd_count_2^B": 2 ** B,
        "lcd_match": lcd == 2 ** B,
        "all_coeffs_nonneg_int": all(isinstance(c, int) and c >= 0 for c in poly),
        "max_dim": zz.zdeg(poly),
        "degree_bound_ok": zz.zdeg(poly) <= sum(d * degree_bound(a, P)
                                                for (a, d) in shapes),
    }


class E2Result:
    """Container for one E2 run."""

    def __init__(self, shapes, P, poly, dist, checks):
        self.shapes = list(shapes)
        self.P = P
        self.poly = poly
        self.dist = dist
        self.checks = checks

    @property
    def n_codes(self) -> int:
        return zz.zval(self.poly, 1)

    def summary(self) -> str:
        flags = [v for k, v in self.checks.items()
                 if k.endswith("_match") or k.endswith("_ok")
                 or k.startswith("all_coeffs")]
        return (f"P={self.P} shapes={self.shapes} "
                f"|C|={self.n_codes} maxdim={zz.zdeg(self.poly)} "
                f"lcd={self.checks['lcd_count_E2']} "
                f"checks={'OK' if all(flags) else 'FAIL'}")


def run_e2(shapes, P, predicted_code_count=None, B=None) -> E2Result:
    """Convenience wrapper: build the enumerator and its structural checks."""
    poly = enumerator_poly(shapes, P)
    dist = {i: c for i, c in enumerate(poly) if c}
    checks = structural_checks(
        shapes, P,
        predicted_code_count if predicted_code_count is not None
        else (P + 1) ** sum(a for a, _ in shapes),
        B if B is not None else len(shapes))
    return E2Result(shapes, P, poly, dist, checks)
