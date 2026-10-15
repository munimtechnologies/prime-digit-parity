"""Decade structure of the base-10 prime digit-sum parity excess.

Work by Munim, Inc., produced with AI. License: Apache-2.0.

Companion code for DECADES.md. It encodes the gated value-cutoff excesses of the
base-10 prime digit-sum parity bias and shows, with ordinary arithmetic, the
empirical law reported there:

    Going from a value cutoff 10^(m-1) to 10^m, the running excess

        E(10^m) = (#primes p <= 10^m with odd decimal digit sum)
                - (#primes p <= 10^m with even decimal digit sum)

    multiplies by roughly one order of magnitude when m is EVEN, but barely
    changes (ratio ~1) when m is ODD. Almost all the growth of the base-10 bias
    happens in the even-indexed decades [10^(2k-1), 10^(2k)); the odd-indexed
    decades [10^(2k), 10^(2k+1)) contribute almost nothing.

Nothing here hits the network or reads a file; the numbers are the gated
sandbox outputs, embedded as literals. See DECADES.md for provenance and limits.
"""

from math import sqrt

# Gated value-cutoff data. excess_base10[m] = E(10^m); n_primes[m] = pi(10^m).
# m = 1..9 come from two independent sandbox sieves that agree at the 10^8 anchor.
EXCESS_BASE10 = {
    1: 2, 2: -1, 3: 6, 4: 75, 5: 112,
    6: 1916, 7: 2967, 8: 62627, 9: 62582,
}
N_PRIMES = {
    1: 4, 2: 25, 3: 168, 4: 1229, 5: 9592,
    6: 78498, 7: 664579, 8: 5761455, 9: 50847534, 10: 455052511, 11: 4118054813,
}
# Published, verified standardized score at the 10^10 value cutoff (base 10).
Z_AT_1E10 = 98.93


def z_score(excess, n):
    """Standardized excess under the fair-coin null: z = excess / sqrt(n)."""
    return excess / sqrt(n)


def decade_ratio(m):
    """E(10^m) / E(10^(m-1)) for m where both excesses are recorded."""
    return EXCESS_BASE10[m] / EXCESS_BASE10[m - 1]


def decade_increment(m):
    """New excess accumulated in the decade [10^(m-1), 10^m)."""
    return EXCESS_BASE10[m] - EXCESS_BASE10[m - 1]


def classify_decade(m):
    """Label a decade 'jump' (ratio > 10) or 'flat' (ratio < 2).

    Returns (label, ratio). The empirical law: 'jump' iff m is even.
    """
    r = decade_ratio(m)
    if r > 10.0:
        return "jump", r
    if r < 2.0:
        return "flat", r
    return "intermediate", r


def law_holds(m_lo=4, m_hi=9):
    """True iff every decade m in [m_lo, m_hi] is 'jump' when m is even and
    'flat' when m is odd -- the even/odd decade law, over recorded data."""
    for m in range(m_lo, m_hi + 1):
        label, _ = classify_decade(m)
        want = "jump" if m % 2 == 0 else "flat"
        if label != want:
            return False
    return True


def predict_flat_z(m_from, z_from):
    """If the decade [10^m_from, 10^(m_from+1)) is flat (excess stays put), the
    standardized score falls purely because sqrt(n) grows:

        z(10^(m+1)) = z(10^m) * sqrt( pi(10^m) / pi(10^(m+1)) ).

    Validated below on the known flat decade 10^8 -> 10^9.
    """
    return z_from * sqrt(N_PRIMES[m_from] / N_PRIMES[m_from + 1])


def leading_digit_sign_sum(base_digit_hi=9):
    """sum_{d=1..9} (-1)^d, the alternating-sign weight the leading digit places
    on a decade via chi_10(p) = (-1)^d * chi_10(p - d*10^(m-1)). Equals -1."""
    return sum((-1) ** d for d in range(1, base_digit_hi + 1))


def digit_sum(n):
    """Decimal digit sum."""
    s = 0
    while n:
        s += n % 10
        n //= 10
    return s


def chi10(n):
    """(-1)^(decimal digit sum of n)."""
    return -1 if digit_sum(n) & 1 else 1


if __name__ == "__main__":
    print("Base-10 prime digit-sum parity excess by decade (value cutoffs 10^m):")
    print(f"{'m':>3} {'E(10^m)':>12} {'ratio':>9} {'incr':>12} {'class':>6}")
    for m in range(4, 10):
        label, r = classify_decade(m)
        print(f"{m:>3} {EXCESS_BASE10[m]:>12} {r:>9.3f} "
              f"{decade_increment(m):>12} {label:>6}")
    print()
    print(f"Even/odd decade law holds for m=4..9: {law_holds()}")
    print(f"sum_(d=1..9) (-1)^d = {leading_digit_sign_sum()}")
    print()
    # Validate the flat-decade z prediction on the known 10^8 -> 10^9 drop.
    pred = predict_flat_z(8, 26.09)
    print(f"Predicted z(10^9) from flat 10^8->10^9: {pred:.3f} "
          f"(observed 8.78)")
    # Pre-registered prediction for a future sieve: 10^10 -> 10^11 is an ODD
    # decade, so if the law holds the excess stays ~2.1e6 and z should fall to:
    print(f"Pre-registered: if [10^10,10^11) is flat, z(10^11) ~ "
          f"{predict_flat_z(10, Z_AT_1E10):.1f} (down from 98.93)")
