# moments.py -- Full moment->cumulant decomposition of the base-2 prime
# digit-sum bias as a four-point function of base-16 digit-sum parities.
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.
#
# Background (see HIERARCHY.md).  For block width b the exact character
# factorisation chi_{2}(n) = prod_{r=0}^{b-1} chi_{2^b}(n >> r) writes the
# Thue-Morse (base-2) digit-sum parity as a product of b shifted base-2^b
# parities.  For b = 4 (base 16),
#
#     chi_2(n) = chi_16(n) * chi_16(n>>1) * chi_16(n>>2) * chi_16(n>>3),
#
# so the standardized base-2 prime excess is a four-point correlation of the
# parities A = chi_16(n), B = chi_16(n>>1), C = chi_16(n>>2), D = chi_16(n>>3).
# HIERARCHY.md stated that for b >= 4 this is no longer a pure top-order
# cumulant.  This module provides the *complete* 15-term bookkeeping and a
# derivation of the b in {2,3} purity boundary from the set partitions alone.

from fractions import Fraction
from math import factorial


# ---- parity characters ----------------------------------------------------

def chi2(n):
    """Thue-Morse character: (-1)^popcount(n), in {+1,-1}."""
    return -1 if (bin(n).count("1") & 1) else 1


def chi16(m):
    """Base-16 digit-sum parity character (-1)^{s_16(m)}, in {+1,-1}."""
    s = 0
    while m:
        s += m & 0xF
        m >>= 4
    return -1 if (s & 1) else 1


def shifted_parities(n):
    """The four shifted base-16 parities (A, B, C, D) at argument n."""
    return (chi16(n), chi16(n >> 1), chi16(n >> 2), chi16(n >> 3))


def factorisation_holds(n):
    """Exact check: chi_2(n) == A*B*C*D for a single n."""
    a, b, c, d = shifted_parities(n)
    return chi2(n) == a * b * c * d


# ---- primes ---------------------------------------------------------------

def primes_below(limit):
    """Simple sieve of Eratosthenes; returns a list of primes < limit."""
    if limit <= 2:
        return []
    sieve = bytearray([1]) * limit
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(limit) if sieve[i]]


# ---- set partitions, moments, cumulants -----------------------------------

def set_partitions(items):
    """Yield every set partition of `items` as a list of blocks (lists)."""
    items = list(items)
    if not items:
        yield []
        return
    if len(items) == 1:
        yield [items]
        return
    first, rest = items[0], items[1:]
    for smaller in set_partitions(rest):
        for i in range(len(smaller)):
            yield smaller[:i] + [[first] + smaller[i]] + smaller[i + 1:]
        yield [[first]] + smaller


def all_subsets(k):
    """Every subset of {0,...,k-1} as a frozenset (including the empty set)."""
    out = []
    for mask in range(1 << k):
        out.append(frozenset(i for i in range(k) if (mask >> i) & 1))
    return out


def moments_from_counts(counts, n, k=4):
    """Empirical joint moments m_S = E[prod_{i in S} X_i] as exact Fractions.

    `counts[pattern]` is the number of samples with sign pattern `pattern`
    (bit i set <=> X_i == -1); `n` is the total number of samples.
    """
    m = {frozenset(): Fraction(1)}
    for S in all_subsets(k):
        if not S:
            continue
        tot = 0
        for pattern, c in counts.items():
            prod = 1
            for i in S:
                if (pattern >> i) & 1:
                    prod = -prod
            tot += c * prod
        m[S] = Fraction(tot, n)
    return m


def cumulant(S, m):
    """Joint cumulant kappa_S from moments via Moebius inversion."""
    total = Fraction(0)
    for part in set_partitions(list(S)):
        kblocks = len(part)
        coeff = ((-1) ** (kblocks - 1)) * factorial(kblocks - 1)
        prod = Fraction(1)
        for block in part:
            prod *= m[frozenset(block)]
        total += coeff * prod
    return total


def all_cumulants(m, k=4):
    """kappa_S for every non-empty subset S of {0,...,k-1}."""
    return {S: cumulant(S, m) for S in all_subsets(k) if S}


def moment_from_cumulants(S, kap):
    """Reconstruct the moment m_S as the partition sum of cumulant products."""
    total = Fraction(0)
    for part in set_partitions(list(S)):
        prod = Fraction(1)
        for block in part:
            prod *= kap[frozenset(block)]
        total += prod
    return total


def singleton_free_partitions(k):
    """Partitions of {0,...,k-1} whose every block has size >= 2."""
    out = []
    for part in set_partitions(range(k)):
        if all(len(block) >= 2 for block in part):
            out.append(part)
    return out
