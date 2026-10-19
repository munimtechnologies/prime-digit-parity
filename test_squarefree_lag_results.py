"""Pure-stdlib checks of the committed squarefree lag-period result.

Verifies data/squarefree_lag_period.json against the pre-registered predictions in
notes/squarefree-lag-period-prereg.md (P-OSC-1..P-OSC-4). No numpy, no network; data is
located relative to this file. Every z is independently RE-DERIVED here from the committed
(S, E_null, var_null) and checked against the job's reported value. The main hypothesis
P-OSC-1 was FALSIFIED and P-OSC-3 was FALSIFIED; these tests assert the falsifications as
well as the confirmations, so the honest negative result is locked to the data.
Run with: python -m unittest -v
Produced with AI by MV2, Munim, Inc.
"""
import json
import math
import os
import statistics
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ODD = ["1", "3", "5", "7"]
EVEN = ["2", "4", "6", "8"]
LOUD = 30.0


def load():
    with open(os.path.join(HERE, "data", "squarefree_lag_period.json"), encoding="utf-8") as f:
        return json.load(f)


def z_of(entry):
    return (entry["S"] - entry["E_null"]) / math.sqrt(entry["var_null"])


class TestSquarefreeLagResults(unittest.TestCase):
    def setUp(self):
        self.d = load()
        self.sf = self.d["per_sieve"]["squarefree"]["x_base2"]["lags"]

    def test_self_tests_passed(self):
        self.assertTrue(self.d["self_tests"]["all_passed"])

    def test_z_rederivation(self):
        # Re-derive z = (S - E_null)/sqrt(var_null) for the squarefree base-2 stream.
        for lag, entry in self.sf.items():
            self.assertTrue(
                math.isclose(z_of(entry), entry["z"], rel_tol=1e-9, abs_tol=1e-9),
                msg=f"lag {lag}: re-derived {z_of(entry)} vs reported {entry['z']}",
            )

    def test_posc2_odd_lags_loud(self):
        # P-OSC-2 HOLDS: every odd lag |z| >= 30.
        for lag in ODD:
            self.assertGreaterEqual(abs(self.sf[lag]["z"]), LOUD)
        self.assertTrue(self.d["grading"]["P_OSC_squarefree"]["P_OSC_2_odd_all_loud"])

    def test_posc1_main_hypothesis_falsified(self):
        # P-OSC-1 FALSIFIED: median(odd|z|) is NOT >= 5 * median(even|z|).
        med_odd = statistics.median(abs(self.sf[l]["z"]) for l in ODD)
        med_even = statistics.median(abs(self.sf[l]["z"]) for l in EVEN)
        self.assertLess(med_odd, 5.0 * med_even)
        self.assertFalse(self.d["grading"]["P_OSC_squarefree"]["P_OSC_1_ratio_ge_5"])

    def test_posc3_even_lags_not_all_quiet_falsified(self):
        # P-OSC-3 FALSIFIED: NOT every even lag is quiet -- lags 4, 6, 8 are loud.
        loud_even = [l for l in EVEN if abs(self.sf[l]["z"]) > LOUD]
        self.assertTrue(set(loud_even) >= {"4", "6", "8"})
        self.assertFalse(self.d["grading"]["P_OSC_squarefree"]["P_OSC_3_even_all_quiet"])

    def test_lag2_is_uniquely_quiet(self):
        # The real structure: lag 2 is the ONLY quiet lag; all of 1,3,4,5,6,7,8 are loud.
        self.assertLess(abs(self.sf["2"]["z"]), LOUD)
        for lag in ["1", "3", "4", "5", "6", "7", "8"]:
            self.assertGreater(abs(self.sf[lag]["z"]), LOUD)

    def test_posc4_xor_specific(self):
        # P-OSC-4 (exploratory): the base-4 streams are not all-quiet on even lags.
        a = self.d["per_sieve"]["squarefree"]["a_base4_even"]["lags"]
        self.assertTrue(any(abs(a[l]["z"]) > LOUD for l in EVEN))
        self.assertTrue(self.d["grading"]["P_OSC_4_xor_specific"]["prediction_holds"])

    def test_counts_present(self):
        c = self.d["counts"]["squarefree"]
        self.assertEqual(c["max_gap"], 10)
        self.assertEqual(c["count"], 60792694)


if __name__ == "__main__":
    unittest.main()
