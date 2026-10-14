"""Reference (pure standard-library) routines for digit-sum parity of primes.

This module has NO third-party dependencies so it can be imported by the test
suite in any environment.  The fast, vectorised version used for the large
computation lives in ``fast.py`` (requires numpy) and is checked against the
functions here on small ranges.

Definitions
-----------
For an integer n >= 0 and an integer base b >= 2, s_b(n) is the sum of the
base-b digits of n.  The *digit-sum parity* is s_b(n) mod 2.

Facts proved in NOTES.md and exercised by the tests:

1. Odd bases are trivial.  If b is odd then n == s_b(n) (mod 2) for every n,
   so every prime p > 2 has an odd base-b digit sum.  Parity studies are only
   meaningful in even bases.

2. Even-base ranges are exactly balanced.  For even b and any k >= 1, exactly
   half of the integers in [0, b**k) have an even base-b digit sum.

3. Even-base excess is bounded by 1 at EVERY cutoff (sharpens Fact 2).  For even
   b and any X, the odd-minus-even digit-sum parity excess over [0, X) lies in
   {-1, 0, +1}, with the closed form in ``integer_parity_excess_closed``.  So
   the integer "boundary term" at an arbitrary cutoff is at most one -- it cannot
   inflate a large bias.

4. Base-10 localization.  Every prime p > 5 ends in a digit in {1,3,7,9} (odd),
   so s_10(p) mod 2 = 1 XOR (s_10(floor(p/10)) mod 2): the whole base-10 prime
   parity bias is carried by the high-order digits (the "beheaded" prime).

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


def integer_parity_excess(limit_exclusive: int, base: int) -> int:
    """Brute-force odd-minus-even digit-sum parity excess over [0, X)."""
    odd, even = integer_parity_counts(limit_exclusive, base)
    return odd - even


def integer_parity_excess_closed(limit_exclusive: int, base: int) -> int:
    """Closed-form odd-minus-even digit-sum parity excess over [0, X), even base.

    Proof (see NOTES.md, Fact 3).  Write X = b*q + t with t = X mod b.  Each
    complete block of b consecutive integers b*m .. b*m + (b-1) has units digit
    running through 0..b-1 while the high part s_b(m) stays fixed, so for EVEN b
    it contributes exactly b/2 odd and b/2 even digit sums -- excess 0.  Only the
    final partial block of t integers (r = 0..t-1, parities s_b(q)+r) can be
    unbalanced.  If t is even the partial block is balanced too; if t is odd it
    has one extra of parity s_b(q).  Hence excess(X) in {-1, 0, +1} for every X.
    """
    if base % 2 != 0:
        raise ValueError("closed form is for even bases only")
    if limit_exclusive <= 0:
        return 0
    q, t = divmod(limit_exclusive, base)
    if t % 2 == 0:
        return 0
    # t odd: the majority parity in the partial block equals s_b(q) mod 2.
    return 1 if (digit_sum(q, base) & 1) else -1
