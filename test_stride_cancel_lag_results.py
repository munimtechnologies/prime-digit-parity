"""Pure-stdlib checks of the committed stride-cancel-lag result.

Verifies data/stride_cancel_lag.json against the pre-registered hypotheses in
notes/stride-cancel-lag-prereg.md (P-L2, P-L3, P-PERSIST) with the frozen
thresholds (stream loud |z|>=30, XOR quiet |z|<=30, base-2 persist |z|>=10).
No numpy, no network; data located relative to this file.
Produced with AI by MV2, Munim, Inc.
"""
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "data", "stride_cancel_lag.json"), encoding="utf-8") as f:
        return json.load(f)


def absz(d, sieve, stream, lag):
    return abs(d["per_sieve"][sieve][stream]["lags"][str(lag)]["z"])


class TestStrideCancelLagResults(unittest.TestCase):
    def setUp(self):
        self.d = load()
        self.t = self.d["grading"]["thresholds"]

    def test_self_tests_passed(self):
        self.assertTrue(self.d["self_tests"]["all_passed"])

    def test_thresholds_frozen(self):
        self.assertEqual(self.t["stream_loud"], 30)
        self.assertEqual(self.t["xor_quiet"], 30)
        self.assertEqual(self.t["base2_persist"], 10)

    def test_pl2_sumsq_cancellation_persists(self):
        # P-L2: at lag-2 a stream stays loud while the base-2 XOR stays quiet.
        g = self.d["grading"]["P_L2_sumsq"]
        self.assertTrue(g["cancellation_persists"])
        self.assertFalse(g["falsified_lag_specific"])
        a = absz(self.d, "sums_of_two_squares", "a_base4_even", 2)
        b = absz(self.d, "sums_of_two_squares", "b_base4_odd", 2)
        x = absz(self.d, "sums_of_two_squares", "x_base2", 2)
        self.assertGreaterEqual(max(a, b), self.t["stream_loud"])
        self.assertLessEqual(x, self.t["xor_quiet"])

    def test_pl3_sumsq_falsified_lag_specific(self):
        # P-L3: at lag-3 the base-2 XOR is no longer quiet -> cancellation breaks.
        g = self.d["grading"]["P_L3_sumsq"]
        self.assertTrue(g["falsified_lag_specific"])
        self.assertFalse(g["cancellation_persists"])
        self.assertGreater(
            absz(self.d, "sums_of_two_squares", "x_base2", 3), self.t["xor_quiet"]
        )

    def test_persist_primes_lag2(self):
        self.assertTrue(self.d["grading"]["P_PERSIST_primes_l2"])
        self.assertGreaterEqual(
            absz(self.d, "primes", "x_base2", 2), self.t["base2_persist"]
        )

    def test_persist_squarefree_lag2_fails(self):
        # Squarefree base-2 does NOT persist at lag-2 (below the 10 threshold).
        self.assertFalse(self.d["grading"]["P_PERSIST_squarefree_l2"])
        self.assertLess(
            absz(self.d, "squarefree", "x_base2", 2), self.t["base2_persist"]
        )

    def test_squarefree_base2_oscillation(self):
        # Loud at lag-1, quiet at lag-2, loud again at lag-3.
        l1 = absz(self.d, "squarefree", "x_base2", 1)
        l2 = absz(self.d, "squarefree", "x_base2", 2)
        l3 = absz(self.d, "squarefree", "x_base2", 3)
        self.assertGreater(l1, 100.0)
        self.assertLess(l2, 10.0)
        self.assertGreater(l3, 100.0)

    def test_grading_absz_consistency(self):
        g2 = self.d["grading"]["P_L2_sumsq"]["abs_z"]
        self.assertAlmostEqual(
            g2["a"], absz(self.d, "sums_of_two_squares", "a_base4_even", 2), places=6
        )
        self.assertAlmostEqual(
            g2["x_base2"], absz(self.d, "sums_of_two_squares", "x_base2", 2), places=6
        )


if __name__ == "__main__":
    unittest.main()
