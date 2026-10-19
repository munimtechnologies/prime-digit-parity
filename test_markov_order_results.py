"""Pure-stdlib checks of the committed Markov-order result.

Verifies data/markov_order.json against the five pre-registered predictions in
notes/markov-order-prereg.md (P-M1..P-M5). No numpy, no network; data is located
relative to this file. Run with: python -m unittest -v
Produced with AI by MV2, Munim, Inc.
"""
import json
import math
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))


def load():
    with open(os.path.join(HERE, "data", "markov_order.json"), encoding="utf-8") as f:
        return json.load(f)


class TestMarkovOrderResults(unittest.TestCase):
    def setUp(self):
        self.d = load()

    def test_self_tests_passed(self):
        st = self.d["self_tests"]
        self.assertTrue(st["all_passed"])
        # A synthetic first-order chain must NOT be flagged; a second-order one MUST.
        self.assertLess(st["markov1_G2"], 40.0)
        self.assertGreater(st["markov2_G2"], 100.0)

    def test_pm1_main_higher_order_base2(self):
        # P-M1: global base-2 test rejects first-order Markov (G2 >= 30).
        pb = self.d["per_base"]["2"]["global"]
        self.assertGreaterEqual(pb["G2"], 30.0)
        self.assertLessEqual(pb["p"], 3.1e-7)

    def test_pm2_blocked_robust_base2(self):
        # P-M2: blocked base-2 reading also rejects (not a drift artifact).
        self.assertLessEqual(self.d["per_base"]["2"]["blocked"]["p"], 1e-3)

    def test_pm3_base2_dominance(self):
        # P-M3: base 2 has the largest global G2 of the six bases.
        g = self.d["global_G2_by_base"]
        self.assertEqual(self.d["base_with_largest_global_G2"], "2")
        for b, v in g.items():
            if b != "2":
                self.assertGreater(g["2"], v)

    def test_pm4_first_order_relation_violated_base2(self):
        # P-M4: base-2 |r2 - r1^2| / |r1^2| > 0.10 (first-order would give r2 = r1^2).
        pb = self.d["per_base"]["2"]
        ratio = abs(pb["r2"] - pb["r1_squared"]) / abs(pb["r1_squared"])
        self.assertGreater(ratio, 0.10)

    def test_pm5_sanity(self):
        # P-M5: gates pass, G2 >= 0, 0 < p <= 1, triples sum to n_gaps - 2.
        n_gaps = self.d["n_gaps"]
        for b, pb in self.d["per_base"].items():
            self.assertGreaterEqual(pb["global"]["G2"], 0.0)
            self.assertTrue(0.0 <= pb["global"]["p"] <= 1.0)
            self.assertEqual(sum(pb["triple_counts"]), pb["n_triples"])
            self.assertEqual(pb["n_triples"], n_gaps - 2)

    def test_all_bases_reject_first_order(self):
        # Beyond P-M1: higher-order memory appears in every base tested.
        for b, pb in self.d["per_base"].items():
            self.assertGreaterEqual(pb["global"]["G2"], 30.0)

    def test_df2_survival_closed_form(self):
        # Chi-square df=2 survival is p = exp(-G2/2); underflows to 0 for large G2.
        for b, pb in self.d["per_base"].items():
            g2 = pb["global"]["G2"]
            if g2 < 1400:
                expected = math.exp(-g2 / 2.0)
                self.assertTrue(
                    math.isclose(pb["global"]["p"], expected, rel_tol=1e-6, abs_tol=1e-300)
                )
            else:
                self.assertEqual(pb["global"]["p"], 0.0)


if __name__ == "__main__":
    unittest.main()
