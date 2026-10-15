"""Self-contained tests for the third-cumulant decomposition (CORRELATION.md).

Stdlib only; no repo imports; no file reads. Every assertion passes because it
is an exact algebraic identity (the pointwise factorisation, the cumulant =
centred-moment identity, and the full-octal-block vanishing), not because it
matches a pre-recorded number.

Work by Munim, Inc., produced with AI. Signed: MV2 of Munim, Inc.  Apache-2.0.
"""

import unittest


def digit_sum_parity(n, base):
    s = 0
    while n:
        s += n % base
        n //= base
    return s & 1


def chi(n, base):
    return 1 - 2 * digit_sum_parity(n, base)


def sieve(limit):
    bs = bytearray([1]) * limit
    bs[0] = bs[1] = 0
    i = 2
    while i * i < limit:
        if bs[i]:
            bs[i * i :: i] = bytearray(len(bs[i * i :: i]))
        i += 1
    return [i for i in range(limit) if bs[i]]


def counts_over(values):
    c = [0] * 8
    for p in values:
        idx = (digit_sum_parity(p, 8) << 2) | (
            digit_sum_parity(p // 2, 8) << 1) | digit_sum_parity(p // 4, 8)
        c[idx] += 1
    return c


def moments(counts):
    n = sum(counts)
    sA = sB = sC = sAB = sAC = sBC = sABC = 0
    for idx, cnt in enumerate(counts):
        A = 1 - 2 * ((idx >> 2) & 1)
        B = 1 - 2 * ((idx >> 1) & 1)
        C = 1 - 2 * (idx & 1)
        sA += cnt * A
        sB += cnt * B
        sC += cnt * C
        sAB += cnt * A * B
        sAC += cnt * A * C
        sBC += cnt * B * C
        sABC += cnt * A * B * C
    return {"N": n, "a": sA / n, "b": sB / n, "c": sC / n,
            "E_AB": sAB / n, "E_AC": sAC / n, "E_BC": sBC / n,
            "E_ABC": sABC / n, "sABC": sABC}


def kappa3(m):
    return (m["E_ABC"] - m["a"] * m["E_BC"] - m["b"] * m["E_AC"]
            - m["c"] * m["E_AB"] + 2.0 * m["a"] * m["b"] * m["c"])


def centred_triple(counts, m):
    n = sum(counts)
    a, b, c = m["a"], m["b"], m["c"]
    t = 0.0
    for idx, cnt in enumerate(counts):
        A = 1 - 2 * ((idx >> 2) & 1)
        B = 1 - 2 * ((idx >> 1) & 1)
        C = 1 - 2 * (idx & 1)
        t += cnt * (A - a) * (B - b) * (C - c)
    return t / n


class TestFactorisation(unittest.TestCase):
    def test_pointwise_identity(self):
        # chi_2(n) == chi_8(n) * chi_8(n//2) * chi_8(n//4) for all n < 2^16.
        for n in range(1 << 16):
            self.assertEqual(chi(n, 2),
                             chi(n, 8) * chi(n // 2, 8) * chi(n // 4, 8))


class TestCumulantIdentity(unittest.TestCase):
    def test_kappa_equals_centred_moment(self):
        # The formula kappa_3 and the centred third moment agree exactly,
        # on a non-trivial dataset (primes below 1e5).
        counts = counts_over(sieve(100000))
        m = moments(counts)
        self.assertAlmostEqual(kappa3(m), centred_triple(counts, m), places=12)

    def test_triple_mean_equals_minus_base2_excess(self):
        # sum A*B*C over primes == sum chi_2(p) (the negative base-2 excess),
        # by the pointwise identity -- an exact integer equality.
        ps = sieve(100000)
        counts = counts_over(ps)
        m = moments(counts)
        self.assertEqual(m["sABC"], sum(chi(p, 2) for p in ps))


class TestFullBlockBaseline(unittest.TestCase):
    def test_triple_and_marginal_vanish_on_full_octal_block(self):
        # [0, 2^15) = [0, 8^5): triple correlation and the chi_8 marginal are 0.
        N = 1 << 15
        s_triple = sum(chi(n, 8) * chi(n // 2, 8) * chi(n // 4, 8)
                       for n in range(N))
        s_chi2 = sum(chi(n, 2) for n in range(N))
        s_chi8 = sum(chi(n, 8) for n in range(N))
        self.assertEqual(s_triple, 0)
        self.assertEqual(s_chi2, 0)
        self.assertEqual(s_triple, s_chi2)
        self.assertEqual(s_chi8, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
