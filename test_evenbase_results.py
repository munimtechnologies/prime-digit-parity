"""Self-contained verification of EVENBASE-RESULTS.md.

Reproduces, with an independent pure-stdlib sieve to 1e6, the even-base prime
digit-sum parity excess E_b and the per-truncation-layer factor means
eps(floor(p/b^i)) for bases 6, 8, 10, and checks the three claims the note
makes about the 1e6 slice:

  * the anchor gate: E_b(1e6) = {6: 1936, 8: -510, 10: 1916};
  * the shipped factor means at 1e6 (data/evenbase_factor_means.csv, x=1e6 rows)
    are reproduced exactly (values embedded inline here so the test reads no file);
  * P-even-2 in miniature: multiplying the measured marginal factor means as if
    the truncation factors were independent predicts an excess of |E| < 1 for
    every base, while the true |E| is in the hundreds -- so the excess lives in
    the joint correlation of the factors, not their marginals.

Also checks the even-base reduction identity (E') chi_b(n)=prod_i eps(floor(n/b^i))
against a direct digit-sum for small n. Stdlib only; no repo imports; no file
reads. Run: python -m unittest -v test_evenbase_results

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Apache-2.0.
"""
import math
import unittest

BASES = (6, 8, 10)
LIMIT = 1_000_000

E_ANCHOR = {6: 1936, 8: -510, 10: 1916}
PI_ANCHOR = 78498

# Factor means eps(floor(p/b^i)) over primes p<=1e6, i>=1, from the anchor-gated
# sandbox run (compete evenbase-factors-v2, x=1e6 slice). Inlined, not read from
# the CSV, so the test is standalone; the CSV ships the same numbers for 1e6/7/8.
EXPECTED_FACTOR_MEANS_1E6 = {
    6: {1: -0.001872730747181349, 2: -0.0012868373106374304,
        3: -0.002740564173815503, 4: -0.003193337420805232,
        5: -0.004270251441693651, 6: -0.024173079794237006,
        7: -0.23023432154019716},
    8: {1: -7.643896348765511e-05, 2: 0.0013761467889908258,
        3: 0.0026657823242050484, 4: -0.005748453819898889,
        5: -0.0205638385832022, 6: -0.29525388302281164},
    10: {1: -0.002777282340051469, 2: -0.0017967963503370587,
         3: -0.00015319800842589047, 4: -0.009693408740892208,
         5: -0.115926044176124},
}


def sieve(limit):
    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i::i] = b"\x00" * len(range(i * i, limit + 1, i))
        i += 1
    return s


def digit_sum_parity_direct(n, b):
    s = 0
    while n:
        s ^= (n % b) & 1
        n //= b
    return s


def chi_reduction(n, b):
    """chi_b(n) via the even-base reduction prod_{i>=0} eps(floor(n/b^i))."""
    q, parity = n, 0
    while q:
        parity ^= q & 1
        q //= b
    return 1 if parity == 0 else -1


def measure(limit):
    s = sieve(limit)
    excess = {b: 0 for b in BASES}
    fsum = {b: {} for b in BASES}
    fcnt = {b: {} for b in BASES}
    n_primes = 0
    for p in range(2, limit + 1):
        if not s[p]:
            continue
        n_primes += 1
        for b in BASES:
            q, i, parity = p, 0, 0
            while q:
                bit = q & 1
                if i >= 1:
                    f = 1 - 2 * bit
                    fsum[b][i] = fsum[b].get(i, 0) + f
                    fcnt[b][i] = fcnt[b].get(i, 0) + 1
                parity ^= bit
                q //= b
                i += 1
            excess[b] += 1 if parity else -1
    means = {b: {i: fsum[b][i] / fcnt[b][i] for i in fsum[b]} for b in BASES}
    return excess, means, n_primes


class TestEvenbaseResults(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.excess, cls.means, cls.n_primes = measure(LIMIT)

    def test_identity_reduction_small_n(self):
        for b in (2, 4, 6, 8, 10, 12, 16):
            for n in range(0, 40000):
                self.assertEqual(
                    chi_reduction(n, b),
                    1 if digit_sum_parity_direct(n, b) == 0 else -1,
                    msg=f"reduction identity fails at n={n}, b={b}")

    def test_anchor_gate(self):
        self.assertEqual(self.n_primes, PI_ANCHOR)
        for b in BASES:
            self.assertEqual(self.excess[b], E_ANCHOR[b],
                             msg=f"E_{b}(1e6) mismatch")

    def test_factor_means_reproduced(self):
        for b in BASES:
            self.assertEqual(sorted(self.means[b]),
                             sorted(EXPECTED_FACTOR_MEANS_1E6[b]),
                             msg=f"layer set for base {b}")
            for i, want in EXPECTED_FACTOR_MEANS_1E6[b].items():
                got = self.means[b][i]
                self.assertTrue(
                    math.isclose(got, want, rel_tol=1e-9, abs_tol=1e-12),
                    msg=f"factor mean base {b} layer {i}: {got} != {want}")

    def test_interior_factors_balanced(self):
        # Low-order (interior) truncation factors are balanced to ~sqrt size:
        # |mean| <~ a few / sqrt(n). The top (leading) layer is NOT balanced --
        # it is the deterministic leading-truncation bias discussed in the note,
        # an order of magnitude (or more) larger than any interior layer.
        sd = 1.0 / math.sqrt(self.n_primes)  # fair-coin sd of a mean over n draws
        for b in BASES:
            layers = sorted(self.means[b])
            interior = layers[:-2]  # exclude the top two (leading) layers
            max_interior = max(abs(self.means[b][i]) for i in interior)
            for i in interior:
                self.assertLess(abs(self.means[b][i]), 20 * sd,
                                msg=f"interior factor base {b} layer {i} "
                                    f"not sqrt-size")
            top = layers[-1]
            self.assertGreater(abs(self.means[b][top]), 10 * max_interior,
                               msg=f"leading layer base {b} not >> interior")

    def test_p_even_2_excess_is_joint_not_marginal(self):
        # Treating the truncation factors as INDEPENDENT: predicted mean of
        # chi_b(p) = eps(p) * prod_i E[eps_i] = (-1) * prod_i factor_mean[i]
        # (eps(p) = -1 for every odd prime). Predicted excess = -n * mean.
        for b in BASES:
            prod = 1.0
            for i in self.means[b]:
                prod *= self.means[b][i]
            predicted_mean_chi = -prod
            predicted_excess = -self.n_primes * predicted_mean_chi
            self.assertLess(abs(predicted_excess), 1.0,
                            msg=f"independence should predict ~0 excess, base {b}")
            self.assertGreaterEqual(abs(self.excess[b]), 400,
                                    msg=f"true excess should be large, base {b}")
            # The gap between prediction and truth is the joint correlation.
            self.assertGreater(abs(self.excess[b]) - abs(predicted_excess), 400)


if __name__ == "__main__":
    unittest.main(verbosity=2)
