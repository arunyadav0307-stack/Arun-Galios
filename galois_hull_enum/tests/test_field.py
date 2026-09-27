"""Unit tests for core/field.py -- the F_{p^e} kernel."""

import os
import random
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.field import FiniteField, prime_factors, factorize_prime_powers  # noqa: E402

FIELDS = [(2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 8),
          (3, 1), (3, 2), (3, 3), (3, 4),
          (5, 1), (5, 2), (7, 1), (7, 2)]


@pytest.mark.parametrize("p,e", FIELDS)
def test_cardinality_and_encoding(p, e):
    F = FiniteField(p, e)
    assert F.q == p ** e
    assert F.zero == 0 and F.one == 1
    # exactly q distinct encodings, all in [0, q)
    assert list(F.elements()) == list(range(p ** e))
    assert list(F.nonzero_elements()) == list(range(1, p ** e))


@pytest.mark.parametrize("p,e", FIELDS)
def test_prime_subfield_embedding(p, e):
    """An integer c in [0,p) encodes the prime-field element c directly, so
    the CP1 lambdas that live in the prime field carry no primitive-element
    convention. This is what makes the F_9 and F_27 checks exact."""
    F = FiniteField(p, e)
    for c in range(p):
        a = F.from_int(c)
        assert a == c
        assert F.vector(a) == tuple([c] + [0] * (e - 1))
    # char p: p * 1 == 0
    acc = 0
    for _ in range(p):
        acc = F.add(acc, F.one)
    assert acc == 0


@pytest.mark.parametrize("p,e", FIELDS)
def test_ring_axioms(p, e):
    """Associativity, commutativity, distributivity, inverses, on random
    triples (seeded, so failures are reproducible)."""
    F = FiniteField(p, e)
    rng = random.Random(12345 + p * 100 + e)
    els = list(F.elements())
    for _ in range(300):
        a, b, c = rng.choice(els), rng.choice(els), rng.choice(els)
        assert F.add(F.add(a, b), c) == F.add(a, F.add(b, c))
        assert F.add(a, b) == F.add(b, a)
        assert F.add(a, F.zero) == a
        assert F.add(a, F.neg(a)) == F.zero
        assert F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c))
        assert F.mul(a, b) == F.mul(b, a)
        assert F.mul(a, F.one) == a
        if a:
            assert F.mul(a, F.inv(a)) == F.one
            assert F.div(F.mul(a, b), a) == b
        assert F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c))
        assert F.sub(a, b) == F.add(a, F.neg(b))


@pytest.mark.parametrize("p,e", FIELDS)
def test_multiplication_tables_agree_with_raw(p, e):
    """Cross-check the log/antilog fast path against slow polynomial
    multiplication modulo the modulus polynomial."""
    F = FiniteField(p, e)
    rng = random.Random(777)
    els = list(F.elements())
    for _ in range(200):
        a, b = rng.choice(els), rng.choice(els)
        assert F.mul(a, b) == F._mul_raw(a, b)
        assert F.add(a, b) == F._add_raw(a, b)


@pytest.mark.parametrize("p,e", FIELDS)
def test_frobenius_is_ring_hom_of_order_e(p, e):
    """sigma(x) = x^p is a ring homomorphism; sigma^e = id; sigma^0 = id."""
    F = FiniteField(p, e)
    rng = random.Random(4242)
    els = list(F.elements())
    for _ in range(200):
        a, b = rng.choice(els), rng.choice(els)
        assert F.frobenius(F.add(a, b)) == F.add(F.frobenius(a), F.frobenius(b))
        assert F.frobenius(F.mul(a, b)) == F.mul(F.frobenius(a), F.frobenius(b))
        assert F.frobenius(a, 0) == a
        assert F.frobenius(a, e) == a
        # sigma^j(a) == a^{p^j}
        j = rng.randrange(e)
        assert F.frobenius(a, j) == F.powi(a, p ** j)
    # fixed field of sigma is F_p
    for a in F.elements():
        if F.frobenius(a) == a:
            assert a < p


@pytest.mark.parametrize("p,e", FIELDS)
def test_generator_and_orders(p, e):
    F = FiniteField(p, e)
    assert F.multiplicative_order(F.gen) == F.q - 1
    assert F.element_order_is(F.gen, F.q - 1)
    assert F.element_order_is(F.one, 1)
    assert not F.element_order_is(F.zero, 1)
    # every nonzero element has order dividing q-1
    for a in list(F.nonzero_elements())[:50]:
        assert F.powi(a, F.q - 1) == 1
        assert (F.q - 1) % F.multiplicative_order(a) == 0
    # a^(order) == 1 and no smaller exponent works
    for a in list(F.nonzero_elements())[:25]:
        o = F.multiplicative_order(a)
        assert F.powi(a, o) == 1
        for ell in prime_factors(o):
            assert F.powi(a, o // ell) != 1
    assert F.powi(F.gen, -(F.q - 1)) == 1
    assert F.powi(F.gen, 0) == 1


@pytest.mark.parametrize("p,e", FIELDS)
def test_gen_pow_covers_the_whole_group(p, e):
    """The lambda rule in regression_385.yaml enumerates lambda as g^i for
    i = 0..q-2; that must be a bijection onto F_q^*."""
    F = FiniteField(p, e)
    seen = {F.gen_pow(i) for i in range(F.q - 1)}
    assert seen == set(F.nonzero_elements())
    assert len(seen) == F.q - 1


def test_prime_factors():
    assert prime_factors(1) == []
    assert prime_factors(2) == [2]
    assert prime_factors(12) == [2, 3]
    assert prime_factors(255) == [3, 5, 17]
    assert factorize_prime_powers(72) == [(2, 3), (3, 2)]


def test_finite_field_rejects_bad_input():
    with pytest.raises(ValueError):
        FiniteField(4, 1)
    with pytest.raises(ValueError):
        FiniteField(2, 0)
