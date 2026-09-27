"""Finite fields F_{p^e} for the k-Galois hull enumerator.

INDEPENDENT IMPLEMENTATION. Written from scratch for P1. Nothing is imported,
copied, adapted or derived from the throwaway blueprinting prototype,
which lives outside this repository.

Element encoding
----------------
An element is the integer  a = sum_{i=0}^{e-1} c_i p^i  with 0 <= c_i < p,
i.e. the base-p digits of the encoding are the coefficient vector of the
element in the polynomial basis (1, alpha, ..., alpha^{e-1}), where alpha is a
root of the field's modulus polynomial.

Consequences of this encoding (both are used deliberately):
  * 0 is the zero element and 1 is the one element;
  * an integer c in [0, p) encodes the prime-subfield element c DIRECTLY, so
    the lambdas of the CP1 quartet that live in the prime field (F_9: lambda=2,
    F_27: lambda=alpha^13=2) are convention-independent -- no choice of
    primitive element can change them.

Arithmetic uses log / antilog tables plus Zech logarithms, so mul, inv and pow
are O(1) and add is O(1); no polynomial arithmetic is needed per operation.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

__all__ = ["FiniteField", "prime_factors", "factorize_prime_powers"]


# --------------------------------------------------------------------------
# small integer helpers
# --------------------------------------------------------------------------
def prime_factors(n: int) -> List[int]:
    """Distinct prime divisors of n > 0, ascending."""
    out: List[int] = []
    d = 2
    m = n
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0:
                m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        out.append(m)
    return out


def factorize_prime_powers(n: int) -> List[Tuple[int, int]]:
    """[(prime, exponent), ...] for n > 0."""
    out: List[Tuple[int, int]] = []
    for p in prime_factors(n):
        k = 0
        m = n
        while m % p == 0:
            m //= p
            k += 1
        out.append((p, k))
    return out


def _poly_mul_mod_f(f: Sequence[int], g: Sequence[int], mod: Sequence[int], p: int) -> List[int]:
    """Multiply two polynomials over F_p and reduce modulo the monic `mod`."""
    if not f or not g:
        return [0]
    res = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                if b:
                    res[i + j] = (res[i + j] + a * b) % p
    return _poly_rem(res, mod, p)


def _poly_rem(a: Sequence[int], mod: Sequence[int], p: int) -> List[int]:
    """Remainder of a modulo the MONIC polynomial `mod`, over F_p."""
    a = list(a)
    dm = len(mod) - 1
    inv_lead = pow(mod[-1] % p, -1, p) if p > 1 else 0
    while len(a) - 1 >= dm and any(a):
        while a and a[-1] == 0:
            a.pop()
        if len(a) - 1 < dm:
            break
        coeff = (a[-1] * inv_lead) % p
        shift = len(a) - 1 - dm
        for i in range(dm + 1):
            a[shift + i] = (a[shift + i] - coeff * mod[i]) % p
        while a and a[-1] == 0:
            a.pop()
    return a if a else [0]


# --------------------------------------------------------------------------
# the field
# --------------------------------------------------------------------------
class FiniteField:
    """F_{p^e}, p prime, e >= 1."""

    def __init__(self, p: int, e: int, modpoly: Optional[Sequence[int]] = None) -> None:
        for name, val in (("p", p), ("e", e)):
            if not isinstance(val, int):
                raise TypeError(f"{name} must be int")
        if e < 1:
            raise ValueError("e must be >= 1")
        # p is prime iff its only prime factor is p itself (this also rejects
        # p = 1, whose factor list is empty, and composites like 4, whose
        # factor list [2] has length 1 but is not [4]).
        if prime_factors(p) != [p]:
            raise ValueError(f"p={p} is not prime")

        self.p = p
        self.e = e
        self.q = p ** e
        self.order = self.q  # cardinality

        if modpoly is None:
            self.modpoly = self._find_irreducible_poly()
        else:
            self.modpoly = list(modpoly)
            if len(self.modpoly) != e + 1 or self.modpoly[-1] % p != 1:
                raise ValueError("modpoly must be monic of degree e")

        self._build_tables()
        self.zero = 0
        self.one = 1

    # ---------------- construction ----------------
    def _int_to_polyvec(self, a: int) -> List[int]:
        """Encoding integer -> coefficient vector over F_p, low to high."""
        out: List[int] = []
        for _ in range(self.e):
            out.append(a % self.p)
            a //= self.p
        return out

    def _polyvec_to_int(self, v: Sequence[int]) -> int:
        out = 0
        for c in reversed(list(v)):
            out = out * self.p + (c % self.p)
        return out

    def _find_irreducible_poly(self) -> List[int]:
        """First monic irreducible polynomial of degree e over F_p, in a
        deterministic lexicographic order on the coefficient vector.

        Irreducibility test (standard): x^{p^e} = x mod f  AND
        gcd(x^{p^{e/ell}} - x, f) = 1 for every prime ell | e.
        """
        p, e = self.p, self.e
        # iterate over monic polys: coefficient vector (c_0..c_{e-1}), then 1
        total = p ** e
        for code in range(total):
            coeffs: List[int] = []
            m = code
            for _ in range(e):
                coeffs.append(m % p)
                m //= p
            f = coeffs + [1]
            if self._is_irreducible_fp(f):
                return f
        raise RuntimeError(f"no irreducible polynomial of degree {e} over F_{p}")

    def _is_irreducible_fp(self, f: Sequence[int]) -> bool:
        p, e = self.p, self.e
        deg = len(f) - 1
        if deg != e:
            return False
        # f must not be divisible by x  (i.e. constant term nonzero) for the
        # standard test to apply cleanly; x itself is irreducible and handled.
        if f[0] % p == 0:
            return deg == 1  # only x is irreducible with zero constant term

        # h = x^{p^e} mod f  must equal x
        h = self._powmod_fp([0, 1], p ** e, f)
        if h != [0, 1]:
            return False
        for ell in prime_factors(e):
            g = self._powmod_fp([0, 1], p ** (e // ell), f)
            # gcd(g - x, f)
            gm = list(g) + [0] * max(0, 2 - len(g))
            gm[1] = (gm[1] - 1) % p
            while gm and gm[-1] == 0:
                gm.pop()
            if not gm:
                gm = [0]
            d = self._poly_gcd_fp(gm, list(f))
            if len(d) - 1 > 0:
                return False
        return True

    def _powmod_fp(self, base: Sequence[int], exp: int, mod: Sequence[int]) -> List[int]:
        result = [1]
        b = _poly_rem(list(base), mod, self.p)
        while exp > 0:
            if exp & 1:
                result = _poly_rem(_poly_mul_mod_f(result, b, mod, self.p), mod, self.p)
            b = _poly_rem(_poly_mul_mod_f(b, b, mod, self.p), mod, self.p)
            exp >>= 1
        return result

    def _poly_gcd_fp(self, a: Sequence[int], b: Sequence[int]) -> List[int]:
        p = self.p
        a = _poly_rem(list(a), b, p) if len(a) - 1 >= len(b) - 1 else list(a)
        a = list(a)
        b = list(b)
        while any(b):
            a, b = b, _poly_rem(a, b, p)
        # normalise monic
        if not a:
            return [0]
        inv = pow(a[-1] % p, -1, p)
        return [(c * inv) % p for c in a]

    def _build_tables(self) -> None:
        p, e, q = self.p, self.e, self.q
        nzm1 = q - 1

        # ---- multiplicative generator -------------------------------------
        # Find the first element (in the canonical encoding order 1,2,3,...)
        # whose multiplicative order is exactly q-1. Deterministic; no RNG.
        primes = prime_factors(nzm1) if nzm1 > 1 else []
        gen = None
        for cand in range(1, q):
            if nzm1 == 1:
                gen = cand
                break
            ok = True
            for ell in primes:
                if self._pow_fp(cand, nzm1 // ell) == 1:
                    ok = False
                    break
            if ok:
                gen = cand
                break
        if gen is None:
            raise RuntimeError(f"no generator found for F_{q}")
        self.gen = gen

        # exp_table[i] = gen^i  (i = 0..q-2), log_table[gen^i] = i
        exp_table = [0] * nzm1
        log_table = [-1] * q
        cur = 1
        for i in range(nzm1):
            exp_table[i] = cur
            log_table[cur] = i
            cur = self._mul_raw(cur, gen)
        self.exp_table = exp_table
        self.log_table = log_table

        # Zech logarithms: zech[d] = log(1 + gen^d), or -1 when 1+gen^d == 0
        zech = [-1] * nzm1
        one = 1
        for d in range(nzm1):
            s = self._add_raw(one, exp_table[d])
            zech[d] = -1 if s == 0 else log_table[s]
        self.zech = zech

        # Frobenius tables: frob[j][a] = a^{p^j}, j = 0..e-1
        frob: List[List[int]] = []
        for j in range(e):
            tbl = [0] * q
            pe = p ** j
            for a in range(q):
                if a == 0:
                    tbl[a] = 0
                else:
                    tbl[a] = exp_table[(log_table[a] * pe) % nzm1]
            frob.append(tbl)
        self.frob_tables = frob

        # negation: -a = a * (-1); -1 = gen^{(q-1)/2} for odd p, and -1 = 1
        # in characteristic 2.
        minus_one = 1 if p == 2 else exp_table[nzm1 // 2]
        neg = [0] * q
        for a in range(q):
            neg[a] = self._mul_raw(a, minus_one)
        self.neg_table = neg
        self.minus_one = minus_one

    # ---- raw helpers that do NOT depend on the tables being complete ----
    def _add_raw(self, a: int, b: int) -> int:
        """Digitwise addition mod p (base-p digits are the coefficients)."""
        p = self.p
        res = 0
        place = 1
        for _ in range(self.e):
            s = (a % p) + (b % p)
            if s >= p:
                s -= p
            res += s * place
            a //= p
            b //= p
            place *= p
        return res

    def _mul_raw(self, a: int, b: int) -> int:
        """Polynomial multiplication modulo self.modpoly."""
        va = self._int_to_polyvec(a)
        vb = self._int_to_polyvec(b)
        res = [0] * (2 * self.e - 1)
        p = self.p
        for i, ca in enumerate(va):
            if ca:
                for j, cb in enumerate(vb):
                    if cb:
                        res[i + j] = (res[i + j] + ca * cb) % p
        rem = _poly_rem(res, self.modpoly, p)
        return self._polyvec_to_int(rem)

    def _pow_fp(self, a: int, k: int) -> int:
        """a^k by square-and-multiply, using _mul_raw (tables not yet built)."""
        result = 1
        base = a
        while k > 0:
            if k & 1:
                result = self._mul_raw(result, base)
            base = self._mul_raw(base, base)
            k >>= 1
        return result

    # ---------------- public arithmetic ----------------
    def add(self, a: int, b: int) -> int:
        if a == 0:
            return b
        if b == 0:
            return a
        la = self.log_table[a]
        lb = self.log_table[b]
        z = self.zech[(lb - la) % (self.q - 1)]
        if z < 0:
            return 0
        return self.exp_table[(la + z) % (self.q - 1)]

    def neg(self, a: int) -> int:
        return self.neg_table[a]

    def sub(self, a: int, b: int) -> int:
        return self.add(a, self.neg_table[b])

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return self.exp_table[(self.log_table[a] + self.log_table[b]) % (self.q - 1)]

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("inverse of zero")
        return self.exp_table[(-self.log_table[a]) % (self.q - 1)]

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))

    def powi(self, a: int, k: int) -> int:
        """a^k for any integer k (negative k allowed when a != 0).

        a = g^log, so a^k = g^{k*log}; Python's % returns a non-negative
        residue, which makes the same expression correct for k < 0.
        """
        if a == 0:
            if k == 0:
                return 1
            if k < 0:
                raise ZeroDivisionError("negative power of zero")
            return 0
        if self.q == 2:
            return 1
        return self.exp_table[(self.log_table[a] * k) % (self.q - 1)]

    def frobenius(self, a: int, times: int = 1) -> int:
        """a^{p^times}.  times is reduced mod e (Frobenius has order e)."""
        return self.frob_tables[times % self.e][a]

    def multiplicative_order(self, a: int) -> int:
        if a == 0:
            raise ValueError("0 has no multiplicative order")
        nzm1 = self.q - 1
        o = nzm1
        for ell, k in factorize_prime_powers(nzm1):
            for _ in range(k):
                if o % ell == 0 and self.powi(a, o // ell) == 1:
                    o //= ell
                else:
                    break
        return o

    def element_order_is(self, a: int, r: int) -> bool:
        """True iff a has multiplicative order exactly r."""
        if a == 0:
            return False
        if r == 1:
            return a == 1
        if self.powi(a, r) != 1:
            return False
        for ell in prime_factors(r):
            if self.powi(a, r // ell) == 1:
                return False
        return True

    # ---------------- conversions and iteration ----------------
    def from_int(self, c: int) -> int:
        """The prime-subfield element c (0 <= c < p)."""
        if not (0 <= c < self.p):
            raise ValueError(f"{c} not in F_{self.p}")
        return c

    def vector(self, a: int) -> Tuple[int, ...]:
        """Coefficient vector in the basis (1, alpha, ..., alpha^{e-1})."""
        return tuple(self._int_to_polyvec(a))

    def from_vector(self, v: Sequence[int]) -> int:
        return self._polyvec_to_int(v)

    def elements(self):
        return range(self.q)

    def nonzero_elements(self):
        return range(1, self.q)

    def gen_pow(self, i: int) -> int:
        """gen^i -- the canonical way to name an element by discrete log."""
        return self.exp_table[i % (self.q - 1)]

    def __repr__(self) -> str:
        return f"FiniteField(p={self.p}, e={self.e}, q={self.q}, modpoly={self.modpoly})"
