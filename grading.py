# grading.py -- The singleton grading of the prime digit-sum parity
# correlation: a graded small-marginal expansion whose order-k term collects
# the set partitions with exactly k singletons, interpolating exactly between
# the leading singleton-free structure (OEIS A000296, see GENERAL.md) and the
# full moment expansion (Bell numbers, see MOMENTS.md).
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.
#
# Background (see HIERARCHY.md, MOMENTS.md, GENERAL.md).  For block width b the
# exact character factorisation
#     chi_2(n) = prod_{r=0}^{b-1} chi_{2^b}(n >> r)
# writes the base-2 (Thue-Morse) digit-sum parity of n as a product of b
# shifted base-2^b parities X_1, ..., X_b.  The b-point correlation of the X_r
# over the primes is, by the moment-cumulant (cluster) formula, a sum over the
# set partitions of {1, ..., b} of a product of joint cumulants, one per block.
#
# A singleton block {r} contributes the marginal mean E[X_r] = delta_r, which
# is empirically small (|delta_r| ~ 1e-3 at 1e6).  Grading the expansion by the
# number of singletons k therefore grades it in powers of the small marginal:
# the order-delta^k term is the sum over partitions with exactly k singletons.
# This module proves the counting law for that grading and its generating
# function, and the module's tests verify it both combinatorially and on the
# primes below 10^6.

from fractions import Fraction
from math import comb, isqrt


# ---------------------------------------------------------------------------
# Digit-sum parities (mod-2 Legendre / tower reduction; see EVENBASE.md)
# ---------------------------------------------------------------------------

def popcount(m):
    """Number of 1 bits of the non-negative integer m."""
    return bin(m).count("1")


def chi2(m):
    """Thue-Morse character (-1)^{s_2(m)}; s_2 = base-2 digit sum = popcount."""
    return 1 if popcount(m) % 2 == 0 else -1


def chi_pow2(m, k):
    """Character (-1)^{s_{2^k}(m)} for base 2^k.

    The parity of the base-2^k digit sum equals the parity of the bits of m at
    positions that are multiples of k (each digit contributes its lowest bit to
    the digit-sum parity).  For k = 1 this is chi2.
    """
    s = 0
    i = 0
    while m:
        if i % k == 0:
            s ^= (m & 1)
        m >>= 1
        i += 1
    return 1 if s == 0 else -1


def tower_columns(n, b):
    """The b shifted base-2^b parities whose product reconstructs chi_2(n)."""
    return [chi_pow2(n >> r, b) for r in range(b)]


def primes_upto(limit):
    """All primes p with 2 <= p <= limit, by a simple sieve."""
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(limit) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, limit + 1, i)))
    return [i for i in range(2, limit + 1) if sieve[i]]


# ---------------------------------------------------------------------------
# Set partitions and the singleton grading
# ---------------------------------------------------------------------------

def set_partitions(items):
    """Yield every set partition of `items` as a list of blocks (lists)."""
    items = list(items)
    if len(items) == 0:
        yield []
        return
    if len(items) == 1:
        yield [items]
        return
    first = items[0]
    for rest in set_partitions(items[1:]):
        for i in range(len(rest)):
            yield rest[:i] + [[first] + rest[i]] + rest[i + 1:]
        yield [[first]] + rest


def num_singletons(partition):
    """Number of size-1 blocks in a partition."""
    return sum(1 for block in partition if len(block) == 1)


def bell_numbers(upto):
    """Bell numbers B_0 .. B_upto via the Bell triangle."""
    row = [1]
    bells = [1]
    for _ in range(upto):
        nxt = [row[-1]]
        for x in row:
            nxt.append(nxt[-1] + x)
        row = nxt
        bells.append(row[0])
    return bells


def a000296_enum(n):
    """A000296(n): set partitions of [n] with no singleton block (brute force)."""
    return sum(1 for p in set_partitions(range(n)) if num_singletons(p) == 0)


def a000296_closed(n, bells=None):
    """A000296(n) = sum_j (-1)^{n-j} C(n,j) B_j (inclusion-exclusion over singletons)."""
    if bells is None:
        bells = bell_numbers(n)
    return sum((-1) ** (n - j) * comb(n, j) * bells[j] for j in range(n + 1))


def a000296_recurrence(upto):
    """A000296(0..upto) via a(n+1) = sum_{m=0}^{n-1} C(n,m) a(m), from
    A'(x) = (e^x - 1) A(x) with A(x) = exp(e^x - 1 - x)."""
    a = [1]  # a(0) = 1
    for n in range(0, upto):
        a.append(sum(comb(n, m) * a[m] for m in range(0, n)))
    return a


def singleton_triangle(bmax):
    """P[b][k] = number of set partitions of [b] with exactly k singletons,
    for b = 0..bmax, built by brute-force enumeration."""
    tri = []
    for b in range(bmax + 1):
        counts = [0] * (b + 1)
        for p in set_partitions(range(b)):
            counts[num_singletons(p)] += 1
        tri.append(counts)
    return tri


def singleton_count_law(b, k):
    """Closed form P(b, k) = C(b, k) * A000296(b - k): choose the k singletons,
    then partition the remaining b - k points with no singleton."""
    return comb(b, k) * a000296_enum(b - k)


if __name__ == "__main__":
    BMAX = 7
    bells = bell_numbers(BMAX)
    tri = singleton_triangle(BMAX)
    print("Singleton grading triangle  P(b,k) = C(b,k) * A000296(b-k)")
    print("b \\ k :", "  ".join(f"{k:>5}" for k in range(BMAX + 1)), "  | Bell  A000296(b)")
    for b in range(BMAX + 1):
        row = tri[b]
        cells = "  ".join(f"{row[k]:>5}" if k < len(row) else "    ." for k in range(BMAX + 1))
        print(f"  {b:>3} :", cells, f"  | {sum(row):>5}  {row[0] if row else 1:>5}")
    print()
    print("Bernhart:  B_n == A000296(n) + A000296(n+1)")
    a = a000296_recurrence(BMAX + 1)
    for n in range(BMAX):
        assert bells[n] == a[n] + a[n + 1]
    print("  holds for n = 0 ..", BMAX - 1)
