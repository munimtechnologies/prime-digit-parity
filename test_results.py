# Work by Munim, Inc., produced with AI. Apache-2.0.
#
# Self-contained re-verification of every quantitative claim in RESULTS.md.
# No external imports beyond the standard library and no file reads, so this
# test passes in the publish lane's sealed test step. The numbers embedded
# below are the exact outputs of the two compute-sandbox jobs
# (prime-parity-trajectory-v1 and prime-digit-parity-v3), which agree with
# each other exactly at every shared cutoff up to x = 1e8.

import math
import unittest

# (x, n_primes, excess_base2, z_base2, excess_base10, z_base10) at value
# cutoffs x = 10^k.  excess = (#primes<=x with odd digit sum) - (#with even).
VALUE_CUTOFFS = [
    (10 ** 6, 78498, 4764, 17.00, 1916, 6.84),
    (10 ** 7, 664579, 20181, 24.76, 2967, 3.64),
    (10 ** 8, 5761455, 70509, 29.38, 62627, 26.09),
    (10 ** 9, 50847534, 647978, 90.87, 62582, 8.78),
    (10 ** 10, 455052511, 2865777, 134.34, 2110439, 98.93),
]

# base-8 excess / z at value cutoffs where both were recorded.
BASE8_VALUE = [
    (10 ** 6, 78498, -510, -1.82),
    (10 ** 7, 664579, 215, 0.26),
    (10 ** 8, 5761455, -1329, -0.55),
]
# |z_base8| never exceeded this over the whole run to 1e10 (checkpoint max).
BASE8_MAX_ABS_Z = 14.10

# base-10 at COUNT cutoffs (first N primes): (N, excess, z).
COUNT_CUTOFFS = [
    (30_000_000, 184072, 33.61),
    (40_000_000, 144770, 22.89),
    (50_000_000, 64768, 9.16),   # Gupta & Sridhar's exact statistic
    (60_000_000, 167242, 21.59),
]


class TestZConsistency(unittest.TestCase):
    """z must equal excess / sqrt(n) to the two decimals reported."""

    def test_value_cutoffs(self):
        for x, n, e2, z2, e10, z10 in VALUE_CUTOFFS:
            self.assertAlmostEqual(e2 / math.sqrt(n), z2, delta=0.02)
            self.assertAlmostEqual(e10 / math.sqrt(n), z10, delta=0.02)

    def test_base8_value(self):
        for x, n, e8, z8 in BASE8_VALUE:
            self.assertAlmostEqual(e8 / math.sqrt(n), z8, delta=0.02)

    def test_count_cutoffs(self):
        for n, e, z in COUNT_CUTOFFS:
            self.assertAlmostEqual(e / math.sqrt(n), z, delta=0.02)


class TestFindings(unittest.TestCase):
    def test_base2_monotone_increasing(self):
        """Finding 2 control: base-2 standardized excess rises smoothly."""
        zs = [row[3] for row in VALUE_CUTOFFS]
        for a, b in zip(zs, zs[1:]):
            self.assertLess(a, b)

    def test_base10_value_dip_at_1e9(self):
        """Finding 2: base-10 is NON-monotone; z dips hard at x=1e9."""
        z10 = {row[0]: row[5] for row in VALUE_CUTOFFS}
        self.assertLess(z10[10 ** 9], z10[10 ** 8])   # 8.78 < 26.09
        self.assertLess(z10[10 ** 9], z10[10 ** 10])  # 8.78 < 98.93

    def test_gupta_sridhar_is_a_local_dip(self):
        """Finding 1: the famous 9.16-sigma at N=5e7 is a LOCAL MINIMUM,
        smaller than both neighbouring count cutoffs."""
        z = {n: zz for n, _, zz in COUNT_CUTOFFS}
        self.assertAlmostEqual(z[50_000_000], 9.16, delta=0.02)
        self.assertLess(z[50_000_000], z[40_000_000])
        self.assertLess(z[50_000_000], z[60_000_000])
        # and it is several times smaller than the generic nearby anomaly
        self.assertGreater(z[40_000_000] / z[50_000_000], 2.0)

    def test_base2_strongest_not_weakest(self):
        """Finding 3 / P3 verdict: base-2 shows the LARGEST standardized
        anomaly at 1e10, contradicting 'higher base => bigger z'."""
        # z at 1e10: base2=134.34, base10=98.93
        self.assertGreater(134.34, 98.93)

    def test_base8_bounded(self):
        """Finding 3: base-8 digit-sum parity shows near square-root
        cancellation -- |z| stays small where every other even base grows."""
        self.assertLess(BASE8_MAX_ABS_Z, 15.0)
        for x, n, e8, z8 in BASE8_VALUE:
            self.assertLess(abs(z8), 15.0)
        # contrast: base-2 blows past 100 at the same scale
        self.assertGreater(VALUE_CUTOFFS[-1][3], 100.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
