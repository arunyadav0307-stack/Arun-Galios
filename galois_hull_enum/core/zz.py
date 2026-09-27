"""Dense polynomials over Z[z] -- exact integer arithmetic, no floats.

INDEPENDENT of sympy: sympy is used only in scripts/run_phase2.py as a
second, independent symbolic check (see task 10 of the P2 brief).

Representation: a tuple of Python ints, index = degree, LOW to HIGH, with no
trailing zero coefficient. The zero polynomial is the empty tuple ().

Nothing here is imported from or adapted from the throwaway blueprinting
prototype, which lives outside this repository.
"""

from __future__ import annotations

from typing import List, Sequence, Tuple

__all__ = [
    "Z", "ONE", "zadd", "zneg", "zsub", "zscal", "zmul", "zpow",
    "zdeg", "zcoeff", "ztrim", "zval", "zformat", "zlist",
    "mat_zero", "mat_eye", "mat_add", "mat_sub", "mat_scal", "mat_mul",
    "mat_trace", "mat_hadamard", "mat_entrywise",
    "zdet", "zmat_charpoly_via_traces",
]

# A Z[z] polynomial is a tuple of ints, low degree first.
Z = Tuple[int, ...]
ONE: Z = (1,)


# --------------------------------------------------------------------------
# polynomial arithmetic
# --------------------------------------------------------------------------
def ztrim(a: Sequence[int]) -> Z:
    out = list(a)
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def zadd(a: Sequence[int], b: Sequence[int]) -> Z:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        if i < len(a):
            out[i] += a[i]
        if i < len(b):
            out[i] += b[i]
    return ztrim(out)


def zneg(a: Sequence[int]) -> Z:
    return tuple(-c for c in a)


def zsub(a: Sequence[int], b: Sequence[int]) -> Z:
    return zadd(a, zneg(b))


def zscal(c: int, a: Sequence[int]) -> Z:
    if c == 0:
        return ()
    return tuple(c * x for x in a)


def zmul(a: Sequence[int], b: Sequence[int]) -> Z:
    if not a or not b:
        return ()
    out = [0] * (len(a) + len(b) - 1)
    for i, ca in enumerate(a):
        if ca:
            for j, cb in enumerate(b):
                if cb:
                    out[i + j] += ca * cb
    return ztrim(out)


def zpow(a: Sequence[int], e: int) -> Z:
    if e < 0:
        raise ValueError("zpow needs a non-negative exponent")
    result = (1,)
    base = tuple(a)
    while e > 0:
        if e & 1:
            result = zmul(result, base)
        e >>= 1
        if e:
            base = zmul(base, base)
    return result


def zdeg(a: Sequence[int]) -> int:
    """Degree; the zero polynomial has degree -1."""
    return len(a) - 1


def zcoeff(a: Sequence[int], i: int) -> int:
    return a[i] if 0 <= i < len(a) else 0


def zval(a: Sequence[int], x: int) -> int:
    """Evaluate at an integer x (Horner)."""
    acc = 0
    for c in reversed(list(a)):
        acc = acc * x + c
    return acc


def zformat(a: Sequence[int]) -> str:
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
            parts.append("z" if c == 1 else f"{cs}z")
        else:
            parts.append(f"z^{i}" if c == 1 else f"{cs}z^{i}")
    return "+".join(parts)


def zlist(a: Sequence[int]) -> List[int]:
    """Coefficient list [c_0, c_1, ...] (zero-padded, no truncation)."""
    return list(a)


# --------------------------------------------------------------------------
# matrices over Z[z]
# --------------------------------------------------------------------------
Matrix = List[List[Z]]


def mat_zero(rows: int, cols: int) -> Matrix:
    return [[() for _ in range(cols)] for _ in range(rows)]


def mat_eye(n: int) -> Matrix:
    m = mat_zero(n, n)
    for i in range(n):
        m[i][i] = (1,)
    return m


def mat_add(A: Matrix, B: Matrix) -> Matrix:
    return [[zadd(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def mat_sub(A: Matrix, B: Matrix) -> Matrix:
    return [[zsub(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def mat_scal(c: int, A: Matrix) -> Matrix:
    return [[zscal(c, A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def mat_mul(A: Matrix, B: Matrix) -> Matrix:
    n, k, m = len(A), len(B), len(B[0])
    out = mat_zero(n, m)
    for i in range(n):
        for l in range(k):
            ail = A[i][l]
            if not ail:
                continue
            Bl = B[l]
            for j in range(m):
                if Bl[j]:
                    out[i][j] = zadd(out[i][j], zmul(ail, Bl[j]))
    return out


def mat_neg_like(A: Matrix) -> Matrix:
    """The additive inverse of A."""
    return [[zneg(A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def mat_trace(A: Matrix) -> Z:
    n = min(len(A), len(A[0]))
    acc: Z = ()
    for i in range(n):
        acc = zadd(acc, A[i][i])
    return acc


def mat_hadamard(A: Matrix, B: Matrix) -> Matrix:
    return [[zmul(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def mat_entrywise(f, A: Matrix) -> Matrix:
    """Apply a Z[z] -> Z[z] function to every entry."""
    return [[f(A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


# --------------------------------------------------------------------------
# determinants (exact, over Z[z])
# --------------------------------------------------------------------------
def _minor(A: Matrix, skip_r: int, skip_c: int) -> Matrix:
    return [[A[i][j] for j in range(len(A[0])) if j != skip_c]
            for i in range(len(A)) if i != skip_r]


def zdet(A: Matrix) -> Z:
    """Determinant by Laplace expansion. Exact; fine for the small matrices
    used here ((P+1) x (P+1), P <= 6)."""
    n = len(A)
    if n == 0:
        return (1,)
    if n == 1:
        return A[0][0]
    total: Z = ()
    for j in range(n):
        if A[0][j]:
            term = zmul(A[0][j], zdet(_minor(A, 0, j)))
            total = zadd(total, term if j % 2 == 0 else zneg(term))
    return total


def zmat_charpoly_via_traces(A: Matrix) -> List[Z]:
    """Characteristic polynomial coefficients of a square matrix over Z[z],
    via Newton's identities from the power sums tr(A^i), i = 1..n.

    chi_A(lambda) = lambda^n - e_1 lambda^{n-1} + ... + (-1)^n e_n, where the
    e_i are the elementary symmetric functions of the eigenvalues, recovered
    from p_i = tr(A^i) by Newton's identities:
        k e_k = sum_{i=1}^{k} (-1)^{i-1} e_{k-i} p_i .
    Returns [c_n, c_{n-1}, ..., c_0] (coefficient of lambda^n first).
    """
    n = len(A)
    if len(A[0]) != n:
        raise ValueError("characteristic polynomial needs a square matrix")
    p = [mat_trace(_matpow(A, i)) for i in range(1, n + 1)]
    e: List[Z] = [(1,)]
    for k in range(1, n + 1):
        acc: Z = ()
        for i in range(1, k + 1):
            term = zmul(e[k - i], p[i - 1])
            acc = zadd(acc, term if i % 2 == 1 else zneg(term))
        # k e_k = acc  ->  divide by the integer k
        if any(c % k != 0 for c in acc):
            raise ArithmeticError("Newton identity produced a non-integer "
                                  "coefficient; matrix is not over a field of "
                                  "characteristic 0")
        e.append(tuple(c // k for c in acc))
    # chi = lambda^n - e_1 lambda^{n-1} + e_2 lambda^{n-2} - ...
    out: List[Z] = []
    for k in range(n, -1, -1):
        ek = e[n - k] if 0 <= n - k <= n else (1,)
        out.append(ek if (n - k) % 2 == 0 else zneg(ek))
    return out


def _matpow(A: Matrix, e: int) -> Matrix:
    if e == 0:
        return mat_eye(len(A))
    result = mat_eye(len(A))
    base = [row[:] for row in A]
    while e > 0:
        if e & 1:
            result = mat_mul(result, base)
        e >>= 1
        if e:
            base = mat_mul(base, base)
    return result
