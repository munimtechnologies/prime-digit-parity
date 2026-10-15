# grading_results.py -- Measured singleton-graded decomposition of the base-2
# prime digit-sum parity excess at 1e8, read off the three base-8 parity
# columns.  Turns the theory in GRADING.md (P(b,k)=C(b,k)*A000296(b-k)) into
# a table of measured magnitudes: how much of E_2(1e8)=70509 each suppression
# order k (number of marginal/singleton factors) actually carries.
# Everything here is exact rational arithmetic on the eight octal-triple counts
# first published in data/base8_cumulant.csv (CORRELATION-RESULTS.md, v0.4.8).
#
# Work by Munim, Inc., produced with AI.  MV2 of Munim, Inc.
# Part of munimtechnologies/prime-digit-parity.
# Code: Apache-2.0.  Text and data: CC BY 4.0.

from fractions import Fraction
from itertools import combinations

# The eight octal-triple pattern counts over the 5,761,455 primes below 1e8.
# Index i = 4*a + 2*b + c, where (a, b, c) are the three base-8 digit-sum
# parity columns of a prime p: a = parity of s_8(p), b = parity of s_8(p>>1),
# c = parity of s_8(p>>2).  A bit value 0 means the column parity is EVEN
# (character chi = +1); bit value 1 means ODD (chi = -1).  The tower reduction
# chi_2(p) = chi_8(p) * chi_8(p>>1) * chi_8(p>>2) makes the product A*B*C the
# base-2 (Thue-Morse) character, so the base-2 excess E_2 = -sum(A*B*C).
COUNTS = [693845, 736247, 734665, 716635, 736921, 719439, 715554, 708149]
N = sum(COUNTS)  # = pi(1e8) = 5761455

# Published anchors this file is gated against (see outcomes / trajectory job).
ANCHOR_N = 5761455
ANCHOR_E2 = 70509    # base-2 excess at 1e8  = -sum(A*B*C)
ANCHOR_E8 = -1329    # base-8 excess at 1e8  = -sum(A)  (A is the chi_8 column)


def _sign(i, which):
    """+1 if the requested column bit of pattern index i is 0, else -1."""
    bit = (i >> (2 - which)) & 1  # which: 0->a (bit2), 1->b (bit1), 2->c (bit0)
    return 1 if bit == 0 else -1


def moment(subset):
    """Exact E[ product of the columns in `subset` ] over the empirical law."""
    total = Fraction(0)
    for i, c in enumerate(COUNTS):
        s = 1
        for w in subset:
            s *= _sign(i, w)
        total += Fraction(s * c, N)
    return total


def cumulant(subset):
    """Connected joint cumulant of the columns in `subset` (size 1, 2 or 3)."""
    s = tuple(subset)
    if len(s) == 1:
        return moment(s)
    if len(s) == 2:
        x, y = s
        return moment((x, y)) - moment((x,)) * moment((y,))
    if len(s) == 3:
        x, y, z = s
        return (moment((x, y, z))
                - moment((x,)) * moment((y, z))
                - moment((y,)) * moment((x, z))
                - moment((z,)) * moment((x, y))
                + 2 * moment((x,)) * moment((y,)) * moment((z,)))
    raise ValueError("subsets of size 1..3 only")


def _set_partitions(items):
    """Yield every set partition of `items` as a list of tuples (blocks)."""
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for smaller in _set_partitions(rest):
        for n, block in enumerate(smaller):
            yield smaller[:n] + [(first,) + block] + smaller[n + 1:]
        yield [(first,)] + smaller


def graded_contributions():
    """E[A*B*C] split by the number of singleton blocks k, as exact Fractions.

    Each set partition of {A,B,C} contributes the product of the connected
    cumulant of its blocks (the moment-cumulant / linked-cluster expansion).
    Grading by the count of size-1 blocks is exactly GRADING.md's P(3,k):
    k=0 -> 1 term (the singleton-free 3-point cumulant), k=1 -> 3 terms,
    k=2 -> 0 terms (P(3,2)=0), k=3 -> 1 term.
    Returned values are contributions to E[A*B*C]; multiply by -N for E_2.
    """
    buckets = {0: Fraction(0), 1: Fraction(0), 2: Fraction(0), 3: Fraction(0)}
    for part in _set_partitions((0, 1, 2)):
        k = sum(1 for block in part if len(block) == 1)
        term = Fraction(1)
        for block in part:
            term *= cumulant(block)
        buckets[k] += term
    return buckets


def main():
    cols = "ABC"
    print(f"N = pi(1e8) = {N}")
    print(f"sum(A) = {int(moment((0,)) * N)}  -> E_8 = {-int(moment((0,)) * N)}")
    print(f"sum(ABC) = {int(moment((0, 1, 2)) * N)}"
          f"  -> E_2 = {-int(moment((0, 1, 2)) * N)}")
    print()
    print("single-column means:")
    for w in range(3):
        print(f"  <{cols[w]}> = {float(moment((w,))):+.6e}")
    print("pair means:")
    for x, y in combinations(range(3), 2):
        print(f"  <{cols[x]}{cols[y]}> = {float(moment((x, y))):+.6e}")
    print(f"triple mean <ABC> = {float(moment((0, 1, 2))):+.6e}")
    print()
    buckets = graded_contributions()
    print("singleton-graded decomposition of E_2 = 70509:")
    for k in (0, 1, 2, 3):
        contrib_e2 = -buckets[k] * N  # E_2 = -sum(ABC)
        share = float(contrib_e2) / ANCHOR_E2 * 100
        print(f"  k={k}: {float(contrib_e2):+12.3f}  ({share:+.5f}% of E_2)")
    total_e2 = -sum(buckets.values()) * N
    print(f"  sum: {float(total_e2):+12.3f}  (= E_2 = {ANCHOR_E2})")


if __name__ == "__main__":
    main()
