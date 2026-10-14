"""Exact fair-coin baseline for the base-10 prime digit-sum parity bias.

Work by Munim, Inc., produced with AI. Licensed under Apache-2.0.

Pure standard library (import math only). These routines pin down the exact
statistic behind Gupta & Sridhar, "A Peculiarity in the Parity of Primes"
(arXiv:1812.11841): among the first 50,000,000 primes they report 25,032,384
with an odd digit sum in base 10. This module turns that raw count into an
odd-minus-even excess and a fair-coin z-score, so the project's headline
number is reproducible and testable, and so the pending multi-cutoff
computation can be read against a fixed, explicit null.

Fair-coin null: treat each item's digit-sum parity as an independent fair
coin. Write each as +1 (odd) or -1 (even); the excess odd-even is then a sum
of n such +/-1 terms, with mean 0 and variance n, hence standard deviation
sqrt(n). The z-score is excess / sqrt(n).
"""

from math import sqrt

# Gupta & Sridhar (arXiv:1812.11841), first 50,000,000 primes, base 10.
GS_N_PRIMES = 50_000_000
GS_ODD_DIGIT_SUM = 25_032_384


def excess(odd, even):
    """Odd-minus-even parity excess."""
    return odd - even


def fair_coin_std(n):
    """Standard deviation of (odd - even) under n independent fair coins."""
    return sqrt(n)


def parity_z(odd, even):
    """z-score of the odd-minus-even excess against the fair-coin null.

    Under the null each of the n = odd + even items is equally likely to have
    odd or even digit sum, so (odd - even) has mean 0 and std sqrt(n). Returns
    0.0 when n == 0.
    """
    n = odd + even
    if n == 0:
        return 0.0
    return (odd - even) / sqrt(n)


def gupta_sridhar():
    """Return (n_primes, odd, even, excess, z) for the published G&S count.

    This is the exact figure the whole study is about: the base-10 digit-sum
    parity of the first 50,000,000 primes.
    """
    odd = GS_ODD_DIGIT_SUM
    even = GS_N_PRIMES - odd
    return GS_N_PRIMES, odd, even, excess(odd, even), parity_z(odd, even)


if __name__ == "__main__":
    n, odd, even, exc, z = gupta_sridhar()
    print(f"first {n:,} primes (base 10)")
    print(f"  odd digit sum : {odd:,}")
    print(f"  even digit sum: {even:,}")
    print(f"  excess (odd-even): {exc:,}")
    print(f"  fair-coin std    : {fair_coin_std(n):.4f}")
    print(f"  z-score          : {z:.4f} sigma")
