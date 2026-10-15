# general.py -- The general-b leading-order structure of the prime digit-sum
# parity correlation: a sum over singleton-free set partitions, counted by
# OEIS A000296, with the purity boundary b in {2, 3} as a corollary.
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.
#
# Background (see HIERARCHY.md, MOMENTS.md).  For block width b the exact
# character factorisation
#     chi_2(n) = prod_{r=0}^{b-1} chi_{2^b}(n >> r)
# writes the base-2 (Thue-Morse) digit-sum parity as a product of b shifted
# base-2^b digit-sum parities A_0, ..., A_{b-1}.  Over the primes the b-point
# correlation E[A_0 ... A_{b-1}] expands (moments in terms of joint cumulants)
# as a sum over ALL set partitions of {0,...,b-1}.  A block of size 1 carries a
# single-character marginal mean, which the prime data makes small (~1e-3), so
# at leading order only singleton-free partitions survive.  Their number is
# A000296(b); it equals 1 only for b in {2,3}, which is exactly the purity
# boundary that MOMENTS.md observed at b <= 4.  This module makes the general
# statement explicit and gives the closed form.
#
# Reference, stdlib only.  The exact tests live in test_general.py.

from fractions import Fraction
from math import comb, factorial


def digit_sum_parity(m, base):
    """Parity of the base-`base` digit sum of m."""
    s = 0
    while m:
        s += m % base
        m //= base
    return s & 1


def chi(m, base):
    """(-1)^{digit-sum of m in base}."""
    return 1 - 2 * digit_sum_parity(m, base)


def tower_columns(n, b):
    """The b shifted base-2^b digit-sum parities whose product is chi_2(n)."""
    base = 1 << b
    return [chi(n >> r, base) for r in range(b)]


def set_partitions(elements):
    """All set partitions of `elements`, each a list of blocks (lists)."""
    elements = list(elements)
    if not elements:
        yield []
        return
    first = elements[0]
    for rest in set_partitions(elements[1:]):
        for i in range(len(rest)):
            yield rest[:i] + [[first] + rest[i]] + rest[i + 1:]
        yield [[first]] + rest


def bell_numbers(upto):
    """B_0 .. B_upto via the Bell triangle."""
    row, bells = [1], [1]
    for _ in range(upto):
        nxt = [row[-1]]
        for x in row:
            nxt.append(nxt[-1] + x)
        row = nxt
        bells.append(row[0])
    return bells


def singleton_free_count_enumerate(n):
    """A000296(n) by brute enumeration (partitions with every block >= 2)."""
    return sum(
        1 for pi in set_partitions(range(n))
        if all(len(block) >= 2 for block in pi)
    )


def singleton_free_count_closed(n):
    """A000296(n) via Bell inclusion-exclusion: sum_k (-1)^k C(n,k) B_{n-k}."""
    bells = bell_numbers(n)
    return sum((-1) ** k * comb(n, k) * bells[n - k] for k in range(n + 1))


def block_size_table(n):
    """Singleton-free partitions of [n] grouped by sorted block-size multiset."""
    table = {}
    for pi in set_partitions(range(n)):
        if all(len(block) >= 2 for block in pi):
            key = tuple(sorted(len(block) for block in pi))
            table[key] = table.get(key, 0) + 1
    return table


def is_pure_top_cumulant(b):
    """True iff the b-point correlation is a single top cumulant at leading
    order, i.e. A000296(b) == 1, i.e. b in {2, 3}."""
    return singleton_free_count_enumerate(b) == 1


def joint_cumulant(block, moment):
    """Joint cumulant of a block from the joint-moment map (Mobius inversion)."""
    tot = Fraction(0)
    for pi in set_partitions(list(block)):
        coef = (-1) ** (len(pi) - 1) * factorial(len(pi) - 1)
        prod = Fraction(1)
        for blk in pi:
            prod *= moment[frozenset(blk)]
        tot += coef * prod
    return tot


if __name__ == "__main__":
    print("b : A000296(b) (enumerate / closed) : pure? : block-size table")
    for b in range(2, 8):
        print(
            b,
            ":", singleton_free_count_enumerate(b),
            "/", singleton_free_count_closed(b),
            ":", is_pure_top_cumulant(b),
            ":", block_size_table(b),
        )
