"""Dense polynomial arithmetic over F_{p^e}.

INDEPENDENT IMPLEMENTATION. Written from scratch for P1; nothing is imported,
copied or adapted from the throwaway blueprinting prototype, which lives
outside this repository.

Representation
--------------
A polynomial is a Python list of field-element encodings, index = degree, in
LOW-to-HIGH order, with no trailing zero coefficients. The zero polynomial is
the empty list []. This normalisation is maintained by every function that
returns a polynomial, so equality is plain list equality.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

__all__ = [
    "deg", "trim", "is_zero", "const", "make_x",
    "padd", "psub", "pneg", "pmul", "pscal", "pmonic", "plead",
    "pdivmod", "pmod", "pgcd", "pegcd", "ppowmod", "pderiv", "peval",
    "frobenius_coeffs", "frobenius_basis", "apply_frobenius",
    "poly_from_int_coeffs", "format_poly",
]


# --------------------------------------------------------------------------
# basics
# --------------------------------------------------------------------------
def deg(a: Sequence[int]) -> int:
    """Degree of a; the zero polynomial has degree -1."""
    return len(a) - 1 if a else -1


def trim(a: Sequence[int]) -> List[int]:
    out = list(a)
    while out and out[-1] == 0:
        out.pop()
    return out


def is_zero(a: Sequence[int]) -> bool:
    return not a


def const(F, c: int) -> List[int]:
    """The constant polynomial c."""
    return [] if c == 0 else [c]


def make_x(F) -> List[int]:
    """The polynomial x."""
    return [0, 1]


def plead(F, a: Sequence[int]) -> int:
    """Leading coefficient (0 for the zero polynomial)."""
    return a[-1] if a else 0


# --------------------------------------------------------------------------
# ring operations
# --------------------------------------------------------------------------
def padd(F, a: Sequence[int], b: Sequence[int]) -> List[int]:
    add = F.add
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        v = 0
        if i < len(a):
            v = a[i]
        if i < len(b):
            v = add(v, b[i])
        out[i] = v
    return trim(out)


def pneg(F, a: Sequence[int]) -> List[int]:
    neg = F.neg
    return trim([neg(c) for c in a])


def psub(F, a: Sequence[int], b: Sequence[int]) -> List[int]:
    return padd(F, a, pneg(F, b))


def pscal(F, c: int, a: Sequence[int]) -> List[int]:
    if c == 0 or not a:
        return []
    mul = F.mul
    return trim([mul(c, x) for x in a])


def pmul(F, a: Sequence[int], b: Sequence[int]) -> List[int]:
    if not a or not b:
        return []
    add, mul = F.add, F.mul
    out = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    out[i + j] = add(out[i + j], mul(ca, cb))
    return trim(out)


def pmonic(F, a: Sequence[int]) -> List[int]:
    """Scale a so that its leading coefficient is 1. Returns [] for zero."""
    if not a:
        return []
    lc = a[-1]
    if lc == 1:
        return trim(list(a))
    inv = F.inv(lc)
    return trim([F.mul(inv, c) for c in a])


# --------------------------------------------------------------------------
# division
# --------------------------------------------------------------------------
def pdivmod(F, a: Sequence[int], b: Sequence[int]) -> Tuple[List[int], List[int]]:
    """Euclidean division a = q*b + r with deg r < deg b."""
    if not b:
        raise ZeroDivisionError("division by the zero polynomial")
    a = trim(list(a))
    b = trim(list(b))
    if not a:
        return [], []
    db = len(b) - 1
    if len(a) - 1 < db:
        return [], a
    sub, mul, inv = F.sub, F.mul, F.inv
    inv_lead = inv(b[-1])
    r = list(a)
    q = [0] * (len(a) - db)
    for i in range(len(a) - db - 1, -1, -1):
        coeff = r[i + db]
        if coeff:
            coeff = mul(coeff, inv_lead)
            q[i] = coeff
            for j in range(db + 1):
                r[i + j] = sub(r[i + j], mul(coeff, b[j]))
    return trim(q), trim(r)


def pmod(F, a: Sequence[int], b: Sequence[int]) -> List[int]:
    return pdivmod(F, a, b)[1]


def pgcd(F, a: Sequence[int], b: Sequence[int]) -> List[int]:
    """Monic gcd."""
    a = trim(list(a))
    b = trim(list(b))
    while b:
        a, b = b, pmod(F, a, b)
    return pmonic(F, a) if a else []


def pegcd(F, a: Sequence[int], b: Sequence[int]):
    """Extended gcd: (g, s, t) with g = s*a + t*b, g monic."""
    r0, r1 = trim(list(a)), trim(list(b))
    s0, s1 = [1], []
    t0, t1 = [], [1]
    while r1:
        q, r = pdivmod(F, r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, psub(F, s0, pmul(F, q, s1))
        t0, t1 = t1, psub(F, t0, pmul(F, q, t1))
    if not r0:
        return [], [], []
    lc = r0[-1]
    inv = F.inv(lc)
    g = [F.mul(inv, c) for c in r0]
    s = [F.mul(inv, c) for c in s0]
    t = [F.mul(inv, c) for c in t0]
    return trim(g), trim(s), trim(t)


def ppowmod(F, base: Sequence[int], exp: int, mod: Sequence[int]) -> List[int]:
    """base^exp mod mod, by square-and-multiply. exp >= 0."""
    if exp < 0:
        raise ValueError("ppowmod needs a non-negative exponent")
    if not mod:
        raise ZeroDivisionError("ppowmod modulo the zero polynomial")
    result = [1]
    b = pmod(F, base, mod)
    while exp > 0:
        if exp & 1:
            result = pmod(F, pmul(F, result, b), mod)
        exp >>= 1
        if exp:
            b = pmod(F, pmul(F, b, b), mod)
    return result


# --------------------------------------------------------------------------
# calculus / evaluation
# --------------------------------------------------------------------------
def pderiv(F, a: Sequence[int]) -> List[int]:
    """Formal derivative."""
    if len(a) <= 1:
        return []
    mul, add = F.mul, F.add
    out = [0] * (len(a) - 1)
    for i in range(1, len(a)):
        # coefficient i * a[i], where i is repeated addition in the field
        acc = 0
        for _ in range(i):
            acc = add(acc, a[i])
        out[i - 1] = acc
    return trim(out)


def peval(F, a: Sequence[int], x: int) -> int:
    """Horner evaluation at the field element x."""
    if not a:
        return 0
    add, mul = F.add, F.mul
    acc = a[-1]
    for c in reversed(a[:-1]):
        acc = add(mul(acc, x), c)
    return acc


# --------------------------------------------------------------------------
# Frobenius on coefficients and on the quotient ring
# --------------------------------------------------------------------------
def frobenius_coeffs(F, a: Sequence[int], times: int = 1) -> List[int]:
    """Apply sigma^times to every coefficient (sigma = Frobenius x -> x^p)."""
    fb = F.frobenius
    return trim([fb(c, times) for c in a])


def frobenius_basis(F, f: Sequence[int]) -> List[List[int]]:
    """[(x^i)^q mod f for i = 0 .. deg(f)-1].

    Because (x^i)^q = (x^q)^i, this is computed by iterating one
    multiplication by x^q mod f. It lets the q-power map on F_q[x]/(f) --
    which is F_q-LINEAR, since coefficients lie in F_q -- be applied to any
    residue class in O(deg^2) field operations instead of O(log q) full
    modular squarings.
    """
    n = len(f) - 1
    if n <= 0:
        return [[]]
    xq = ppowmod(F, [0, 1], F.q, f)
    basis: List[List[int]] = [[1]]
    cur = [1]
    for _ in range(1, n):
        cur = pmod(F, pmul(F, cur, xq), f)
        basis.append(list(cur))
    return basis


def apply_frobenius(F, h: Sequence[int], basis: Sequence[Sequence[int]], f: Sequence[int]) -> List[int]:
    """h^q mod f, using the precomputed Frobenius basis.

    h^q = sum_i h_i * (x^i)^q  mod f, because h_i^q = h_i for h_i in F_q.
    """
    if not h:
        return []
    add, mul, zero = F.add, F.mul, F.zero
    n = len(basis)
    out: List[int] = []
    for i, c in enumerate(h):
        if c and i < n:
            term = [mul(c, x) for x in basis[i]]
            out = padd(F, out, term) if out else term
    return pmod(F, out, f) if out else []


# --------------------------------------------------------------------------
# construction from plain integer coefficients
# --------------------------------------------------------------------------
def poly_from_int_coeffs(F, coeffs: Sequence[int]) -> List[int]:
    """Build a polynomial from integer coefficients given as elements of the
    PRIME field (0 <= c < p). Used for the published CP1 polynomials, whose
    printed coefficients all lie in the prime field and therefore carry no
    primitive-element convention.
    """
    out = [F.from_int(c % F.p) for c in coeffs]
    return trim(out)


def format_poly(F, a: Sequence[int]) -> str:
    """Human-readable polynomial string in STANDARD DESCENDING power order,
    e.g. x^4+2x^3+x^2+2x+1 -- the order used in the source papers, so that
    comparisons against printed polynomials are string-exact.
    """
    if not a:
        return "0"
    parts = []
    for i in range(len(a) - 1, -1, -1):
        c = a[i]
        if c == 0:
            continue
        cs = str(c)
        if i == 0:
            parts.append(cs)
        elif i == 1:
            parts.append("x" if c == 1 else f"{cs}x")
        else:
            parts.append(f"x^{i}" if c == 1 else f"{cs}x^{i}")
    return "+".join(parts)
