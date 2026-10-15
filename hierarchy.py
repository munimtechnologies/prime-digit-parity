"""
hierarchy.py -- the cumulant tower of the base-2 prime digit-sum parity bias.

Part of munimtechnologies/prime-digit-parity. Work by Munim, Inc., produced
with AI. Code: Apache-2.0. The accompanying note HIERARCHY.md is CC BY 4.0.

For q a power of two, q = 2^a, the reduction identity (REDUCE.md, derived from
the positional lemma in HIERARCHY.md section 1) is EXACT:

    s_{2^a}(n) = XOR_{r=0}^{b-1} s_{2^{ab}}( n >> (r*a) )        (mod 2)

Read through the character chi_q(n) = (-1)^{s_q(n)} it becomes a product,

    chi_{2^a}(n) = prod_{r=0}^{b-1} chi_{2^{ab}}( n >> (r*a) ),

so the standardized base-2^a excess is a b-fold joint correlation of shifted
base-2^{ab} parities. Its moment-cumulant expansion isolates a connected b-th
joint cumulant; this module provides the exact (integer / Fraction) tooling that
HIERARCHY.md and test_hierarchy.py use. No value from any sieve is used here.
"""
from fractions import Fraction


# ---- digit-sum parity characters -----------------------------------------

def digit_sum(n, base):
    """Sum of the base-`base` digits of n >= 0."""
    s = 0
    while n:
        s += n % base
        n //= base
    return s


def s_parity(n, q):
    """Parity (0 or 1) of the base-q digit sum of n."""
    return digit_sum(n, q) & 1


def chi(n, q):
    """The parity character chi_q(n) = (-1)^{s_q(n)} in {+1, -1}."""
    return 1 - 2 * s_parity(n, q)


def chi_factors(n, a, b):
    """The b factors chi_{2^{ab}}(n >> r*a); their product equals chi_{2^a}(n)."""
    Q = 1 << (a * b)
    return [chi(n >> (r * a), Q) for r in range(b)]


def reduction_holds(n, a, b):
    """True iff identity (star) holds at n for (a, b). Exact, by construction."""
    prod = 1
    for f in chi_factors(n, a, b):
        prod *= f
    return prod == chi(n, 1 << a)


# ---- exact moment / cumulant algebra on finite +/-1 samples ---------------

def mean(col):
    return Fraction(sum(col), len(col))


def joint_moment(cols):
    """E[ X_0 * X_1 * ... ] over paired samples; `cols` is a list of columns."""
    n = len(cols[0])
    total = 0
    for i in range(n):
        p = 1
        for c in cols:
            p *= c[i]
        total += p
    return Fraction(total, n)


def _factorial(k):
    r = 1
    for i in range(2, k + 1):
        r *= i
    return r


def set_partitions(items):
    """Yield every set partition of `items` (a list) as a list of blocks."""
    items = list(items)
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [[first] + part[i]] + part[i + 1:]
        yield [[first]] + part


def joint_cumulant(cols):
    """
    Joint cumulant kappa(X_0, ..., X_{k-1}) via Leonov-Shiryaev:
        kappa = sum_partitions (-1)^{|pi|-1} (|pi|-1)! prod_blocks E[prod block].
    Exact Fraction.
    """
    k = len(cols)
    total = Fraction(0)
    for part in set_partitions(range(k)):
        m = len(part)
        term = Fraction((-1) ** (m - 1)) * _factorial(m - 1)
        for block in part:
            term *= joint_moment([cols[j] for j in block])
        total += term
    return total


def moment_from_cumulants(cols):
    """Rebuild E[prod] = sum_partitions prod_blocks kappa(block). Exact Fraction."""
    k = len(cols)
    total = Fraction(0)
    for part in set_partitions(range(k)):
        term = Fraction(1)
        for block in part:
            term *= joint_cumulant([cols[j] for j in block])
        total += term
    return total


def connected_plus_pairings(cols):
    """
    Sum over set partitions with NO singleton block. For |cols| in {2,3} this is
    just the top cumulant; for |cols| = 4 it is kappa_4 plus the three pairings.
    This is the leading-order contribution when marginals are small.
    """
    k = len(cols)
    total = Fraction(0)
    for part in set_partitions(range(k)):
        if any(len(block) == 1 for block in part):
            continue
        term = Fraction(1)
        for block in part:
            term *= joint_cumulant([cols[j] for j in block])
        total += term
    return total
