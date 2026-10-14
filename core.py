"""Reference (pure standard-library) routines for digit-sum parity of primes.

This module has NO third-party dependencies so it can be imported by the test
suite in any environment.  The fast, vectorised version used for the large
computation lives in ``fast.py`` (requires numpy) and is checked against the
functions here on small ranges.

Definitions
-----------
For an integer n >= 0 and an integer base b >= 2, s_b(n) is the sum of the
base-b digits of n.  The *digit-sum parity* is s_b(n) mod 2.

Two facts proved in NOTES.md and exercised by the tests:

1. Odd bases are trivial.  If b is odd then n == s_b(n) (mod 2) for every n,
   so every prime p > 2 has an odd base-b digit sum.  Parity studies are only
   meaningful in even bases.

2. Even-base ranges are exactly balanced.  For even b and any k >= 1, exactly
   half of the integers in [0, b**k) have an even base-b digit sum.

Work by Munim, Inc.  Produced with AI (see README.md).  Apache-2.0.
"""

from __future__ import annotations


def digit_sum(n: int, base: int) -> int:
    """Sum of the base-`base` digits of the non-negative integer `n`."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if base < 2:
        raise ValueError("base must be >= 2")
    s = 0
    while n:
        s += n % base
        n //= base
    return s


def digit_sum_parity(n: int, base: int) -> int:
    """Return s_base(n) mod 2 (0 for even digit sum, 1 for odd)."""
    return digit_sum(n, base) & 1


def primes_upto(limit: int) -> list[int]:
    """All primes <= limit via a plain sieve of Eratosthenes (for tests)."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0] = sieve[1] = 0
    i = 2
    while i * i <= limit:
        if sieve[i]:
            sieve[i * i :: i] = b"\x00" * len(sieve[i * i :: i])
        i += 1
    return [i for i in range(2, limit + 1) if sieve[i]]


def prime_parity_counts(limit: int, base: int) -> tuple[int, int]:
    """(odd_count, even_count) of digit-sum parity over primes <= limit."""
    odd = even = 0
    for p in primes_upto(limit):
        if digit_sum_parity(p, base):
            odd += 1
        else:
            even += 1
    return odd, even


def integer_parity_counts(limit_exclusive: int, base: int) -> tuple[int, int]:
    """(odd_count, even_count) of digit-sum parity over 0 <= n < limit."""
    odd = even = 0
    for n in range(limit_exclusive):
        if digit_sum_parity(n, base):
            odd += 1
        else:
            even += 1
    return odd, even


def excess(odd: int, even: int) -> int:
    """Odd-minus-even excess (positive = odd digit sums lead)."""
    return odd - even
