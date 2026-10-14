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

5. Casting out (b-1)'s, and the degeneracy of odd bases.  Because b == 1
   (mod b-1), for every n we have s_b(n) == n (mod b-1).  Hence the digit-sum
   residue modulo any d dividing (b-1) is fixed by n itself.  For parity this
   means: when b is odd (so 2 | b-1) the digit-sum parity equals n's parity and
   carries no information (this is Fact 1); when b is even (b-1 odd) parity is
   NON-degenerate.  Base 10 has b-1 = 9 (odd), so it is a non-degenerate case --
   exactly the regime in which Mauduit & Rivat (Ann. of Math. 171, 2010) proved
   the digit-sum parity of primes equidistributes.  See RESOLUTION.md.

6. A bounded control for the base-10 prime bias.  Integers (not only primes) whose
   last digit is in {1,3,7,9} -- i.e. coprime to 10 -- have a base-10 digit-sum
   parity excess over [0, X) of at most 8 in absolute value at EVERY X (closed form
   in ``coprime10_parity_excess_closed``).  Primes p > 5 are a subset of these
   integers (Fact 4), so any prime excess that grows like sqrt(n) is a property of
   the primes inside that bounded ambient set, not of the ambient set itself.

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


def parity_degenerate(base: int) -> bool:
    """True if digit-sum parity is trivially fixed by n's parity (Fact 5).

    By casting out (base-1)'s, s_base(n) == n (mod base-1) for every n (because
    base == 1 mod base-1).  When base-1 is even -- i.e. base is odd -- this forces
    s_base(n) == n (mod 2), so every prime p > 2 has digit-sum parity 1 and the
    statistic carries no information.  When base is even (base-1 is odd), parity
    is non-degenerate: the regime of the Mauduit-Rivat equidistribution theorem.
    """
    if base < 2:
        raise ValueError("base must be >= 2")
    return base % 2 == 1


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


_COPRIME10_DIGITS = (1, 3, 7, 9)


def coprime10_parity_excess(limit_exclusive: int) -> int:
    """Brute force: odd-minus-even base-10 digit-sum parity excess over the
    integers n in [0, X) whose last digit is in {1,3,7,9}."""
    exc = 0
    for n in range(max(limit_exclusive, 0)):
        if n % 10 in _COPRIME10_DIGITS:
            exc += 1 if digit_sum_parity(n, 10) else -1
    return exc


def coprime10_parity_excess_closed(limit_exclusive: int) -> int:
    """Closed form of ``coprime10_parity_excess`` (Fact 6).

    Write n = 10*m + d with d in {1,3,7,9}.  Every such d is odd, so
    s_10(n) = s_10(m) + d has parity 1 XOR s_10(m).  Write X = 10*M + r.  Each
    complete block m < M contributes 4 integers, all of parity 1 XOR s_10(m): the
    excess over those blocks is 4 * (#even s_10(m) - #odd s_10(m)) = -4 * E(M),
    where E(M) in {-1,0,1} is the integer excess of Fact 3.  The partial block
    m = M contributes c = #{d in {1,3,7,9} : d < r} integers, each of parity
    1 XOR s_10(M).  Hence |excess| <= 4 + 4 = 8 for every X.
    """
    if limit_exclusive <= 0:
        return 0
    m_top, r = divmod(limit_exclusive, 10)
    full = -4 * integer_parity_excess_closed(m_top, 10)
    c = sum(1 for d in _COPRIME10_DIGITS if d < r)
    partial = c * (1 if digit_sum_parity(m_top, 10) == 0 else -1)
    return full + partial
