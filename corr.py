"""Third-cumulant decomposition of the base-2 prime digit-sum bias.

Companion to CORRELATION.md. Stdlib only. Exact.

The base-2 parity character factors pointwise as
    chi_2(n) = chi_8(n) * chi_8(n//2) * chi_8(n//4),
so the standardised base-2 prime excess equals minus the triple correlation of
the three shifted base-8 parities over primes. This module computes, from the
eight pattern counts of (s_8(p), s_8(p//2), s_8(p//4)) over a set of primes, the
marginal/pairwise/triple moments and the third joint cumulant kappa_3, which (to
leading order in the small base-8 marginals) carries the whole base-2 bias.

Work by Munim, Inc., produced with AI. Signed: MV2 of Munim, Inc.
Code: Apache-2.0. Text/data: CC BY 4.0.
"""

from typing import Dict, Iterable, List, Tuple


def digit_sum_parity(n: int, base: int) -> int:
    """Parity (0 or 1) of the base-`base` digit sum of n >= 0."""
    if n < 0:
        raise ValueError("n must be non-negative")
    if base < 2:
        raise ValueError("base must be >= 2")
    s = 0
    while n:
        s += n % base
        n //= base
    return s & 1


def chi(n: int, base: int) -> int:
    """The +/-1 digit-sum parity character (-1)^{s_base(n)}."""
    return 1 - 2 * digit_sum_parity(n, base)


def triple_pattern(p: int) -> Tuple[int, int, int]:
    """(s_8(p), s_8(p//2), s_8(p//4)) as 0/1 parities."""
    return (
        digit_sum_parity(p, 8),
        digit_sum_parity(p // 2, 8),
        digit_sum_parity(p // 4, 8),
    )


def pattern_index(sa: int, sb: int, sc: int) -> int:
    return (sa << 2) | (sb << 1) | sc


def accumulate(primes: Iterable[int]) -> List[int]:
    """Eight pattern counts, indexed by pattern_index(s_8(p),s_8(p//2),s_8(p//4))."""
    counts = [0] * 8
    for p in primes:
        sa, sb, sc = triple_pattern(p)
        counts[pattern_index(sa, sb, sc)] += 1
    return counts


def moments_from_counts(counts: List[int]) -> Dict[str, float]:
    """Raw moments a,b,c and E[AB],E[AC],E[BC],E[ABC] from the 8 pattern counts.

    A = chi_8(p) = 1-2*s_a, etc.; expectations are sample means over the primes.
    """
    n = sum(counts)
    if n == 0:
        raise ValueError("no primes supplied")
    sA = sB = sC = sAB = sAC = sBC = sABC = 0
    for idx, c in enumerate(counts):
        sa = (idx >> 2) & 1
        sb = (idx >> 1) & 1
        sc = idx & 1
        A = 1 - 2 * sa
        B = 1 - 2 * sb
        C = 1 - 2 * sc
        sA += c * A
        sB += c * B
        sC += c * C
        sAB += c * A * B
        sAC += c * A * C
        sBC += c * B * C
        sABC += c * A * B * C
    return {
        "N": float(n),
        "a": sA / n,
        "b": sB / n,
        "c": sC / n,
        "E_AB": sAB / n,
        "E_AC": sAC / n,
        "E_BC": sBC / n,
        "E_ABC": sABC / n,
    }


def third_cumulant(m: Dict[str, float]) -> float:
    """kappa_3 = E[ABC] - a E[BC] - b E[AC] - c E[AB] + 2 a b c."""
    return (
        m["E_ABC"]
        - m["a"] * m["E_BC"]
        - m["b"] * m["E_AC"]
        - m["c"] * m["E_AB"]
        + 2.0 * m["a"] * m["b"] * m["c"]
    )


def centred_triple_moment(counts: List[int], m: Dict[str, float]) -> float:
    """E[(A-a)(B-b)(C-c)] directly; equals third_cumulant(m) exactly."""
    n = sum(counts)
    a, b, c = m["a"], m["b"], m["c"]
    total = 0.0
    for idx, cnt in enumerate(counts):
        A = 1 - 2 * ((idx >> 2) & 1)
        B = 1 - 2 * ((idx >> 1) & 1)
        C = 1 - 2 * (idx & 1)
        total += cnt * (A - a) * (B - b) * (C - c)
    return total / n


def _sieve(limit: int) -> List[int]:
    """Primes < limit (simple sieve, for the self-check only)."""
    if limit < 3:
        return []
    bs = bytearray([1]) * limit
    bs[0] = bs[1] = 0
    i = 2
    while i * i < limit:
        if bs[i]:
            bs[i * i :: i] = bytearray(len(bs[i * i :: i]))
        i += 1
    return [i for i in range(limit) if bs[i]]


if __name__ == "__main__":
    # Pointwise identity chi_2 = chi_8 * chi_8(.//2) * chi_8(.//4).
    for n in range(1 << 14):
        assert chi(n, 2) == chi(n, 8) * chi(n // 2, 8) * chi(n // 4, 8)
    # Cumulant == centred third moment, on the primes below 1e5.
    ps = _sieve(100000)
    counts = accumulate(ps)
    mom = moments_from_counts(counts)
    k_formula = third_cumulant(mom)
    k_direct = centred_triple_moment(counts, mom)
    assert abs(k_formula - k_direct) < 1e-12
    print("corr.py self-check OK; primes < 1e5:",
          "E[ABC]=%.6f kappa_3=%.6f a,b,c=%.4f,%.4f,%.4f"
          % (mom["E_ABC"], k_formula, mom["a"], mom["b"], mom["c"]))
