"""Cross-axis additivity + consistency audit of the base-2 gap-parity density law.

Self-contained, pure-stdlib verifier. It embeds ONLY values already published in
this repository across three independently sealed compute jobs:

  * Axis D5 (v0.29.0, data/axis_d5_density_scaling_law.json):
    the squarefree control ratio r_control and the single-prime ratios r(p).
  * Axis D6 (v0.30.0, data/axis_d6_density_additivity.json):
    the two-prime ratios r({p,q}) and the published summed density terms.
  * Axis D7 (v0.31.0, data/axis_d7_three_prime_additivity.json):
    the three-prime ratios r({p,q,s}) and the published summed density terms.

It re-derives, with no network and no sealed job, three facts that tie the three
axes together and that no previously committed test checks directly:

  A. r_control is bit-identical across D5, D6 and D7 (same deterministic strict
     squarefree control at N = 20,000,000).
  B. The admitted-class density term d(p) = 1/p^2 - 1/p^3 summed over each pair
     and triple -- including the fully out-of-sample triple {41,43,47} -- exactly
     reproduces the summed-density terms published in D6 and D7. This confirms the
     density inputs to the k coefficients are the analytic admitted-class density,
     not a fitted quantity.
  C. The MEASURED multi-prime excess r(S) - r_control equals the SUM of the
     single-prime excesses r(p) - r_control to within 1% for every pair and triple
     whose single primes were measured in D5 (3 pairs, 2 triples). This is a
     direct, magnitude-level test of additivity that does NOT route through the
     coefficient band [5,11]; it combines three separately sealed jobs.

Produced by AI agents working as MV2, Munim, Inc. Apache-2.0.
"""

import math
import unittest

# --- Strict squarefree control, published identically in D5, D6, D7 ----------
R_CONTROL_D5 = 0.014860313424210063
R_CONTROL_D6 = 0.014860313424210063
R_CONTROL_D7 = 0.014860313424210063
R_CONTROL = R_CONTROL_D5

# --- D5 single-prime ratios r(p) (data/axis_d5_density_scaling_law.json) ------
# Only the primes reused by the D6 pairs and D7 triples are needed here.
SINGLE_R = {
    17: 0.042025554336370244,
    19: 0.03664851374134167,
    23: 0.029122410406888552,
    29: 0.02388594394194274,
    31: 0.022755927382039735,
    37: 0.020843655330602732,
}

# --- D6 two-prime ratios and published summed density (axis_d6_*.json) --------
PAIR_R = {
    (17, 19): 0.06386828364765362,
    (23, 29): 0.03823642979618743,
    (31, 37): 0.028700461097034637,
}
PAIR_SUMMED_DENSITY = {
    (17, 19): 0.0058809552431880255,
    (23, 29): 0.0029562281901740674,
    (31, 37): 0.0017177335642310492,
}

# --- D7 three-prime ratios and published summed density (axis_d7_*.json) ------
TRIPLE_R = {
    (17, 19, 23): 0.07859868242995202,
    (29, 31, 37): 0.03776340449917874,
    (41, 43, 47): 0.027759210358165386,  # fully out-of-sample: no single r(p)
}
TRIPLE_SUMMED_DENSITY = {
    (17, 19, 23): 0.007689124882375992,
    (29, 31, 37): 0.0028657921152171494,
    (41, 43, 47): 0.0015516917548888875,
}

# Pre-registered tolerances (fixed before looking at the residuals).
DENSITY_ABS_TOL = 1e-9     # analytic identity vs published summed density
EXCESS_REL_TOL = 0.01      # measured multi-prime excess vs sum of single excesses


def admitted_density(p):
    """Density of the class admitted by relaxing the squarefree cap at prime p:
    numbers with p-adic valuation exactly 2 (p^2 | n, p^3 not | n)."""
    return 1.0 / p ** 2 - 1.0 / p ** 3


def excess(r):
    return r - R_CONTROL


class TestControlIdentity(unittest.TestCase):
    def test_r_control_bit_identical_across_axes(self):
        # A deterministic control at the same N must give the same ratio in every
        # axis; a mismatch would mean the axes are not comparable.
        self.assertEqual(R_CONTROL_D5, R_CONTROL_D6)
        self.assertEqual(R_CONTROL_D6, R_CONTROL_D7)


class TestDensityTermAnalytic(unittest.TestCase):
    def test_pair_densities_match_identity(self):
        for (p, q), published in PAIR_SUMMED_DENSITY.items():
            got = admitted_density(p) + admitted_density(q)
            self.assertAlmostEqual(
                got, published, delta=DENSITY_ABS_TOL,
                msg="pair {}: identity {} vs published {}".format((p, q), got, published),
            )

    def test_triple_densities_match_identity(self):
        # Includes {41,43,47}, which was never measured singly: the density term
        # is analytic, so it is checkable even out-of-sample.
        for (p, q, s), published in TRIPLE_SUMMED_DENSITY.items():
            got = admitted_density(p) + admitted_density(q) + admitted_density(s)
            self.assertAlmostEqual(
                got, published, delta=DENSITY_ABS_TOL,
                msg="triple {}: identity {} vs published {}".format((p, q, s), got, published),
            )


class TestExcessAdditivity(unittest.TestCase):
    """r(S) - r_control == sum over p in S of (r(p) - r_control), to EXCESS_REL_TOL,
    using single r(p) from D5 and multi-prime r(S) from D6/D7 -- three separate jobs."""

    def _check(self, members, r_multi):
        self.assertTrue(all(p in SINGLE_R for p in members),
                        "single r(p) must be in-sample for {}".format(members))
        measured = excess(r_multi)
        summed = sum(excess(SINGLE_R[p]) for p in members)
        rel = abs(measured - summed) / abs(measured)
        self.assertLess(
            rel, EXCESS_REL_TOL,
            msg="{}: measured excess {} vs sum of singles {} (rel {:.4%})".format(
                members, measured, summed, rel),
        )
        return rel

    def test_pairs_are_additive(self):
        for pair, r_multi in PAIR_R.items():
            self._check(pair, r_multi)

    def test_in_sample_triples_are_additive(self):
        for triple, r_multi in TRIPLE_R.items():
            if all(p in SINGLE_R for p in triple):
                self._check(triple, r_multi)

    def test_out_of_sample_triple_is_excluded(self):
        # {41,43,47} has no published single r(p); the audit must not fabricate one.
        self.assertEqual(
            [t for t in TRIPLE_R if not all(p in SINGLE_R for p in t)],
            [(41, 43, 47)],
        )


if __name__ == "__main__":
    # Print the realized residuals for the record, then run the assertions.
    print("r_control (D5=D6=D7):", R_CONTROL)
    for pair, r_multi in PAIR_R.items():
        m, s = excess(r_multi), sum(excess(SINGLE_R[p]) for p in pair)
        print("pair  {}: measured {:.10f} sum {:.10f} rel {:.4%}".format(
            pair, m, s, abs(m - s) / abs(m)))
    for triple, r_multi in TRIPLE_R.items():
        if all(p in SINGLE_R for p in triple):
            m, s = excess(r_multi), sum(excess(SINGLE_R[p]) for p in triple)
            print("triple {}: measured {:.10f} sum {:.10f} rel {:.4%}".format(
                triple, m, s, abs(m - s) / abs(m)))
    unittest.main(verbosity=2)
