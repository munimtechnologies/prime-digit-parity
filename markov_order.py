"""Core functions for the Markov-order control on prime-gap digit-sum parity.

The prime-digit-parity arc established that the base-b digit-sum parity sequence of
consecutive prime gaps carries serial memory (significant lag-1..lag-k
autocorrelation). But a *first-order* Markov process already produces non-zero
autocorrelation at every lag (r_k = r_1^k), so significant multi-lag correlation
does NOT by itself prove memory beyond one step. This module tests the Markov
ORDER of the parity sequence directly, via the classical triple
conditional-independence test (Anderson-Goodman 1957; Billingsley 1961):

    H0 (first-order Markov):  s_{n+1}  is independent of  s_{n-1}  given  s_n.

Rejecting H0 means the serial memory is genuinely higher-order.

Everything here is pure Python (no numpy) so it can be validated by a fast,
dependency-free test suite; the sandbox job provides numpy fast paths that are
cross-checked against these reference implementations.

Produced with AI by MV2, Munim, Inc. Apache-2.0.
"""
from __future__ import annotations

import math


def sieve_primes(n):
    """Return a list of all primes <= n (Sieve of Eratosthenes, pure Python)."""
    if n < 2:
        return []
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    i = 2
    while i * i <= n:
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
        i += 1
    return [i for i in range(2, n + 1) if sieve[i]]


def digit_sum_parity(value, base):
    """Parity (0/1) of the sum of the base-`base` digits of a non-negative int."""
    if base < 2:
        raise ValueError("base must be >= 2")
    v = int(value)
    if v < 0:
        raise ValueError("value must be >= 0")
    s = 0
    while v:
        s += v % base
        v //= base
    return s & 1


def parity_table(max_value, base):
    """Lookup list of digit_sum_parity(v, base) for v in 0..max_value inclusive."""
    return [digit_sum_parity(v, base) for v in range(max_value + 1)]


def lag_k_autocorr(bits, k):
    """Sample Pearson autocorrelation at lag k of a 0/1 sequence (coded as +-1)."""
    n = len(bits)
    if n <= k:
        return 0.0
    u = [1.0 - 2.0 * b for b in bits]
    m = sum(u) / n
    u = [x - m for x in u]
    denom = sum(x * x for x in u)
    if denom == 0.0:
        return 0.0
    num = sum(u[i] * u[i + k] for i in range(n - k))
    return num / denom


def triple_counts(bits):
    """Counts N[a,b,c] of consecutive triples (s_{n-1}, s_n, s_{n+1}).

    Returns an 8-element list indexed by 4*a + 2*b + c.
    """
    counts = [0] * 8
    for n in range(len(bits) - 2):
        counts[4 * bits[n] + 2 * bits[n + 1] + bits[n + 2]] += 1
    return counts


def g_test_triples(counts):
    """Triple conditional-independence test of first-order Markov.

    counts: 8-element sequence indexed 4*a + 2*b + c, with a = s_{n-1},
    b = s_n (the conditioning state), c = s_{n+1}.

    H0: s_{n+1} independent of s_{n-1} given s_n. For each middle state b we
    test independence in the 2x2 table over (a, c); the statistics add across
    the two values of b, with df = 2. Returns a dict with the likelihood-ratio
    G^2, the Pearson X^2, df (=2) and p = exp(-G^2 / 2) (exact chi^2_2 survival).
    """
    N = [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
    for idx, cnt in enumerate(counts):
        a = (idx >> 2) & 1
        b = (idx >> 1) & 1
        c = idx & 1
        N[a][b][c] = int(cnt)
    g2 = 0.0
    x2 = 0.0
    for b in (0, 1):  # condition on the middle state
        total = sum(N[a][b][c] for a in (0, 1) for c in (0, 1))
        if total == 0:
            continue
        for a in (0, 1):
            ra = N[a][b][0] + N[a][b][1]
            for c in (0, 1):
                cc = N[0][b][c] + N[1][b][c]
                obs = N[a][b][c]
                exp = ra * cc / total
                if exp > 0:
                    if obs > 0:
                        g2 += 2.0 * obs * math.log(obs / exp)
                    x2 += (obs - exp) ** 2 / exp
    return {"G2": g2, "X2": x2, "df": 2, "p": math.exp(-g2 / 2.0)}


def gammaincc(a, x):
    """Regularized upper incomplete gamma Q(a, x) (Numerical Recipes gser/gcf)."""
    if x < 0 or a <= 0:
        raise ValueError("bad args to gammaincc")
    if x == 0:
        return 1.0
    if x < a + 1.0:
        # Series for the lower P(a, x); return 1 - P.
        ap = a
        term = 1.0 / a
        s = term
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            s += term
            if abs(term) < abs(s) * 1e-16:
                break
        return 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a))
    # Continued fraction (Lentz) for the upper Q(a, x).
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def chi2_sf(x, df):
    """Survival function (upper tail) of chi-square with df degrees of freedom."""
    return gammaincc(df / 2.0, x / 2.0)
