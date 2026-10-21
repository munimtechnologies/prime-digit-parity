# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc.
#
# In-sample non-monotonicity audit of the density-scaling-law coefficient k(p)
# for the base-2 gap-parity "lag-2 unique-quiet" signature (prime-digit-parity).
#
# Produced with AI by MV2 of Munim, Inc. Frontier/reputation work; $0 earned,
# zero external engagement. Fully re-derivable offline with the standard library
# and no network access.
#
# WHAT THIS CHECKS (pure stdlib, no sealed job, no magnitudes-that-drift):
#   The published density-scaling law says, for the admitted-class relaxation at
#   a single prime p (cap 2),
#         r(p) - r_control  ~  C * d(p),   d(p) = 1/p^2 - 1/p^3 = (p-1)/p^3,
#   with per-prime coefficient k(p) = (r(p) - r_control) / d(p) reported in
#   data/axis_d5_density_scaling_law.json for the twelve primes 2..37.
#
#   This test re-derives every k(p) from the published ratios and the EXACT
#   density (p-1)/p^3, then pins down a feature census-in-sample (149) did not
#   test: k(p) is NON-MONOTONE across the measured primes, with two reversals
#   that survive the data's precision limits:
#       (bottom)  k(2)  = 9.017  <  k(3)  = 10.235
#       (top)     k(31) = 7.841  <  k(37) = 8.419     [inside the full-precision
#                                                      D5-exact block 29,31,37]
#   No single-term monotone correction model k(p) = C + a*g(p) with g monotone
#   in p (e.g. g = 1/p or g = d(p)) can produce either reversal. All assertions
#   are fixed published constants with pre-registered tolerances; none depend on
#   N, so a different run length cannot make them pass or fail by accident.

import json
import math
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- Published constants (verbatim from data/axis_d5_density_scaling_law.json).
# r_control is bit-identical across axes D5/D6/D7.
R_CONTROL = 0.014860313424210063

# grading.ratio_by_relaxed_prime : primes 2,3,5,7 are stored at 3 decimal places
# (a coarse public table); 11..37 carry full double precision; 29,31,37 are the
# "D5 data (exact)" block.
RATIO = {
    2: 1.142,
    3: 0.773,
    5: 0.342,
    7: 0.169,
    11: 0.07770937284253937,
    13: 0.06033990564490313,
    17: 0.042025554336370244,
    19: 0.03664851374134167,
    23: 0.029122410406888552,
    29: 0.02388594394194274,
    31: 0.022755927382039735,
    37: 0.020843655330602732,
}

# grading.density_scaling_law.coefficient_k_by_prime (the column we re-derive).
K_PUB = {
    2: 9.01711749260632,
    3: 10.234885768773166,
    5: 10.223115205493436,
    7: 8.81165208258266,
    11: 8.36520980857963,
    13: 8.326555342405221,
    17: 8.341426787590185,
    19: 8.302514776400317,
    23: 7.887587908556781,
    29: 7.8616465248922225,
    31: 7.840607847256792,
    37: 8.418728266236329,
}

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
SMALL = [2, 3, 5, 7]            # primes whose ratios are the coarse 3-d.p. table
LARGE = [11, 13, 17, 19, 23, 29, 31, 37]
EXACT_BLOCK = [29, 31, 37]      # full-precision "D5 data (exact)"


def density(p):
    """Exact admitted-class density d(p) = 1/p^2 - 1/p^3 = (p-1)/p^3."""
    return (p - 1) / float(p ** 3)


def k_derived(p):
    return (RATIO[p] - R_CONTROL) / density(p)


def rounding_unc_k(p, decimals=3):
    """Half-ULP of `decimals`-place rounding of r(p), propagated to k(p)."""
    half = 0.5 * 10 ** (-decimals)
    return half / density(p)


def median(values):
    s = sorted(values)
    n = len(s)
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


class ExactDensity(unittest.TestCase):
    def test_density_is_pminus1_over_pcubed(self):
        # The density term is exact rational (p-1)/p^3; confirm the closed form
        # reproduces the three verbatim example densities from the file.
        self.assertEqual(density(2), 0.125)
        self.assertAlmostEqual(density(7), 0.017492711370262388, delta=1e-15)
        self.assertAlmostEqual(density(37), 0.0007107180226245237, delta=1e-15)


class CoefficientRederivation(unittest.TestCase):
    def test_k_matches_published_column(self):
        # Internal-consistency audit: the published (ratio, density, k) triple
        # must satisfy k = (r - r_control)/d for every one of the 12 primes.
        for p in PRIMES:
            self.assertAlmostEqual(
                k_derived(p), K_PUB[p], delta=1e-9,
                msg="k(%d) re-derivation disagrees with published column" % p)


class SmallVsLargeSeparation(unittest.TestCase):
    def test_small_primes_do_not_overlap_large(self):
        # Every prime <= 7 has a strictly larger coefficient than every prime
        # >= 11: min over {2,3,5,7} exceeds max over {11..37}. Clean gap.
        lo_small = min(K_PUB[p] for p in SMALL)
        hi_large = max(K_PUB[p] for p in LARGE)
        self.assertGreater(lo_small, hi_large)
        self.assertGreater(lo_small - hi_large, 0.3)  # margin ~0.39


class BottomReversal(unittest.TestCase):
    def test_k2_below_k3_robust_to_rounding(self):
        # The smallest prime does NOT have the most extreme k: k(2) < k(3).
        # The gap dwarfs the combined 3-d.p. rounding uncertainty on both, so
        # it is not a precision artifact of the coarse small-prime table.
        gap = K_PUB[3] - K_PUB[2]
        self.assertGreater(gap, 0.0)
        unc = rounding_unc_k(2) + rounding_unc_k(3)
        self.assertGreater(gap, 10.0 * unc)


class TopReversalInExactBlock(unittest.TestCase):
    def test_small_ratios_are_three_decimals(self):
        for p in SMALL:
            self.assertEqual(round(RATIO[p], 3), RATIO[p])

    def test_exact_block_ratios_are_full_precision(self):
        for p in EXACT_BLOCK:
            self.assertNotEqual(round(RATIO[p], 3), RATIO[p])

    def test_sequence_29_31_37_is_down_then_up(self):
        # Within the full-precision D5-exact block the coefficient decreases
        # 29 -> 31 and then jumps UP 31 -> 37. Not attributable to 3-d.p.
        # rounding (these ratios are exact) nor to plausible sampling scale.
        self.assertGreater(K_PUB[29], K_PUB[31])          # 29 -> 31 down
        self.assertLess(K_PUB[31], K_PUB[37])             # 31 -> 37 up
        self.assertGreater(K_PUB[37] - K_PUB[31], 0.3)    # +0.578, ~10x noise scale


class MidPrimeFlatness(unittest.TestCase):
    def test_mid_band_near_constant(self):
        # Locally (11..19) the coefficient IS near-constant: relative spread <1%.
        mid = [K_PUB[p] for p in (11, 13, 17, 19)]
        rel = (max(mid) - min(mid)) / (sum(mid) / len(mid))
        self.assertLess(rel, 0.01)


class NotUniversal(unittest.TestCase):
    def test_full_range_spread_exceeds_quarter_of_median(self):
        ks = [K_PUB[p] for p in PRIMES]
        spread = max(ks) - min(ks)
        self.assertGreater(spread / median(ks), 0.25)  # ~0.287


class MonotoneCorrectionRejected(unittest.TestCase):
    def test_at_least_two_increases_with_prime(self):
        # Any k(p) = C + a*g(p) with a>0 and g monotone decreasing in p (1/p,
        # d(p), p^-s) is monotone decreasing in p; a>0 with g increasing is
        # monotone increasing. Either way k would be monotone. The data has
        # increases at BOTH ends (2->3 and 31->37), so no single-term monotone
        # correction fits. Count sign changes in consecutive differences.
        diffs = [K_PUB[PRIMES[i + 1]] - K_PUB[PRIMES[i]]
                 for i in range(len(PRIMES) - 1)]
        increases = [d for d in diffs if d > 0]
        self.assertGreaterEqual(len(increases), 2)
        # and the two robust ones specifically:
        self.assertGreater(K_PUB[3] - K_PUB[2], 0.0)
        self.assertGreater(K_PUB[37] - K_PUB[31], 0.0)


class CrossCheckCommittedFile(unittest.TestCase):
    """If the committed data file is present in the run, cross-check that the
    embedded constants are bit-faithful to it (catches any transcription error
    before anything ships). Skips cleanly when the file is not in the checkout."""

    def _load(self):
        path = os.path.join(HERE, "data", "axis_d5_density_scaling_law.json")
        if not os.path.exists(path):
            self.skipTest("committed data file not present in this run")
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)

    def test_constants_match_file(self):
        doc = self._load()
        g = doc["grading"]
        dsl = g["density_scaling_law"]
        self.assertTrue(
            math.isclose(dsl["r_control_measured"], R_CONTROL,
                         rel_tol=0.0, abs_tol=1e-15))
        ratio_map = g["ratio_by_relaxed_prime"]
        k_map = dsl["coefficient_k_by_prime"]
        dens_map = dsl["density_by_prime"]
        for p in PRIMES:
            s = str(p)
            self.assertTrue(math.isclose(ratio_map[s], RATIO[p],
                                         rel_tol=0.0, abs_tol=1e-12),
                            msg="ratio(%d) differs from file" % p)
            self.assertTrue(math.isclose(k_map[s], K_PUB[p],
                                         rel_tol=0.0, abs_tol=1e-9),
                            msg="k(%d) differs from file" % p)
            self.assertTrue(math.isclose(dens_map[s], density(p),
                                         rel_tol=1e-12, abs_tol=0.0),
                            msg="density(%d) differs from (p-1)/p^3" % p)


if __name__ == "__main__":
    unittest.main(verbosity=2)
