# leading.py -- positional (leading-digit) decomposition of the base-10
# prime digit-sum excess.
#
# Work by Munim, Inc., produced with AI. Apache-2.0 (code).
# Accompanying write-up LEADING.md and any data are CC BY 4.0.
#
# Context. Gupta & Sridhar (arXiv:1812.11841) observed that the base-10
# digit sum of primes is odd more often than even. This repository's
# DECADES.md (v0.4.6) found that the running excess
#
#     E(x) = #{p <= x : s(p) odd} - #{p <= x : s(p) even}
#
# grows by a large factor across a decade [10^(m-1), 10^m) when the number
# of digits m is EVEN, but barely changes when m is ODD. The mechanism was
# left open. This module gives an exact positional decomposition of E that
# (a) unifies the decade study with the even-base reduction in EVENBASE.md,
# and (b) reduces the open even/odd law to a precise statement about nine
# leading-digit band correlations T_d(m).
#
# Everything here is elementary and exact; no claim is numerical beyond what
# the accompanying test recomputes with a stdlib sieve and gates against the
# published anchors E(10^m) for m <= 6.

from typing import Dict, List


# --- base-10 digit-sum character ------------------------------------------

def s10(n: int) -> int:
    """Base-10 digit sum of n >= 0."""
    t, acc = n, 0
    while t:
        acc += t % 10
        t //= 10
    return acc


def chi(n: int) -> int:
    """Digit-sum parity character chi(n) = (-1)^{s10(n)} in {+1, -1}."""
    return 1 if s10(n) % 2 == 0 else -1


# --- the two classical building blocks, stated exactly ---------------------

def chi_via_truncations(n: int) -> int:
    """chi(n) = prod_{i>=0} (-1)^{floor(n / 10^i)}.

    Proof. The Legendre-type identity s10(n) = n - 9 * sum_{i>=1} floor(n/10^i)
    gives, modulo 2 (9 is odd), s10(n) == sum_{i>=0} floor(n/10^i) (mod 2),
    since floor(n/10^0) = n. Exponentiating by -1 yields the product form.

    This is exactly the uniform even-base reduction of EVENBASE.md
    specialised to b = 10, so the decade study and the even-base thread are
    the same object. Finite because floor(n/10^i) = 0 once 10^i > n.
    """
    prod, t = 1, n
    # i = 0 term is (-1)^n; successive terms use the running truncation t.
    while True:
        if t % 2:
            prod = -prod
        if t == 0:
            break
        t //= 10
    return prod


def chi_leading(n: int) -> int:
    """chi(n) = (-1)^d * chi(n mod 10^(m-1)) where d is the leading digit and
    m the number of digits of n >= 1.

    Proof. n = d*10^(m-1) + r with 0 <= r < 10^(m-1); digit sums add, so
    s10(n) = d + s10(r), hence chi(n) = (-1)^d chi(r).
    """
    if n == 0:
        return 1
    m = len(str(n))
    p = 10 ** (m - 1)
    d = n // p
    r = n % p
    return (1 if d % 2 == 0 else -1) * chi(r)


# --- prime sieve (stdlib) --------------------------------------------------

def primes_upto(limit: int) -> List[int]:
    """All primes <= limit by a simple sieve of Eratosthenes."""
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    i = 2
    while i * i <= limit:
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        i += 1
    return [i for i in range(2, limit + 1) if sieve[i]]


# --- the excess and its positional decomposition ---------------------------

def excess_upto(x: int, primes: List[int]) -> int:
    """E(x) = #{p <= x : s odd} - #{p <= x : s even} = - sum_{p<=x} chi(p)."""
    acc = 0
    for p in primes:
        if p > x:
            break
        acc -= chi(p)
    return acc


def leading_profile(m: int, primes: List[int]) -> Dict[int, int]:
    """T_d(m) = sum over primes p in [10^(m-1), 10^m) with leading digit d of
    chi(r), where r = p - d*10^(m-1) is the trailing (m-1)-digit part.

    Returns {d: T_d(m)} for d in 1..9. Requires `primes` to cover 10^m.
    """
    lo, hi = 10 ** (m - 1), 10 ** m
    base = 10 ** (m - 1)
    T = {d: 0 for d in range(1, 10)}
    for p in primes:
        if p < lo:
            continue
        if p >= hi:
            break
        d = p // base
        r = p - d * base
        T[d] += chi(r)
    return T


def decade_excess_from_profile(T: Dict[int, int]) -> int:
    """E_decade(m) = E(10^m) - E(10^(m-1)) = sum_{d=1}^9 (-1)^(d+1) T_d(m).

    Derivation. For p in the decade, chi(p) = (-1)^d chi(r); summing
    -chi over the decade and grouping by leading digit d gives
    - sum_d (-1)^d T_d = sum_d (-1)^(d+1) T_d.
    """
    return sum((1 if d % 2 else -1) * T[d] for d, _ in sorted(T.items()))
    # (-1)^(d+1) = +1 for odd d, -1 for even d.


def integer_decade_character_sum(m: int) -> int:
    """sum_{n in [10^(m-1), 10^m)} chi(n), the 'all integers' baseline.

    Equals 0 for every m >= 2: grouping by leading digit d gives
    sum_d (-1)^d * S(m-1) with S(k) = sum_{r=0}^{10^k - 1} chi(r) = 0 for
    k >= 1 (a full base-10 block, where sum_{e=0}^9 (-1)^e = 0 in each
    position). So the entire excess is a *primes* effect.
    """
    lo, hi = 10 ** (m - 1), 10 ** m
    return sum(chi(n) for n in range(lo, hi))
