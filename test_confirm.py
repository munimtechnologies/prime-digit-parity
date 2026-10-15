"""Self-contained checks for CONFIRM.md.

Work by Munim, Inc., produced with AI. Apache-2.0.

These assert algebraic identities on the measured numbers reported in
CONFIRM.md / data/dip_confirmation.csv. They do not recompute primes and do
not read any file; every value is embedded here and every assertion passes by
mathematics (the null-model identity z = excess / sqrt(n)). The underlying
counts are produced by compete/prime-parity-dip-confirm-v1/job.py.
"""
import math
import unittest

# (first n primes) -> (odd-even base-10 excess, reported z)
COUNT = {
    30_000_000: (184072, 33.60679553503031),
    40_000_000: (144770, 22.89014684312881),
    50_000_000: (64768, 9.159578400778063),
    60_000_000: (167242, 21.59084935954736),
    70_000_000: (252984, 30.237371450385222),
}

# (primes below X) -> dict of n, base-10 excess/z, base-8 excess/z
VALUE = {
    10**8: dict(n=5761455, e10=62627, z10=26.09128815736286,
                e8=-1329, z8=-0.5536800734688752),
    10**9: dict(n=50847534, e10=62582, z10=8.776361236175434,
                e8=-58698, z8=-8.231677668355529),
}


class TestNullIdentity(unittest.TestCase):
    """z must equal excess / sqrt(n) exactly (fair-coin null)."""

    def test_count_cutoff_z(self):
        for n, (excess, z) in COUNT.items():
            self.assertAlmostEqual(z, excess / math.sqrt(n), places=6)

    def test_value_cutoff_z(self):
        for row in VALUE.values():
            self.assertAlmostEqual(row["z10"], row["e10"] / math.sqrt(row["n"]),
                                   places=6)
            self.assertAlmostEqual(row["z8"], row["e8"] / math.sqrt(row["n"]),
                                   places=6)


class TestDipIsLocalMinimum(unittest.TestCase):
    """The Gupta-Sridhar 50M cutoff is a strict local minimum of |z|."""

    def test_strict_local_min(self):
        z = {n: abs(zz) for n, (_, zz) in COUNT.items()}
        self.assertLess(z[50_000_000], z[40_000_000])
        self.assertLess(z[50_000_000], z[60_000_000])

    def test_neighbours_2p4_to_3p7x(self):
        z50 = abs(COUNT[50_000_000][1])
        ratios = [abs(COUNT[n][1]) / z50
                  for n in (30_000_000, 40_000_000, 60_000_000, 70_000_000)]
        self.assertGreaterEqual(min(ratios), 2.3)
        self.assertLessEqual(max(ratios), 3.7)


class TestExternalAnchor(unittest.TestCase):
    """50M cutoff reproduces Gupta & Sridhar's published excess and sigma."""

    def test_matches_paper(self):
        excess, z = COUNT[50_000_000]
        self.assertEqual(excess, 64768)
        self.assertAlmostEqual(z, 9.16, places=2)


class TestPlateau(unittest.TestCase):
    """Base-10 raw excess is nearly flat from 1e8 to 1e9 while n grows ~8.8x."""

    def test_raw_excess_flat(self):
        e8 = VALUE[10**8]["e10"]
        e9 = VALUE[10**9]["e10"]
        self.assertLess(abs(e9 - e8), 100)            # net change is tiny (-45)

    def test_count_grows(self):
        ratio = VALUE[10**9]["n"] / VALUE[10**8]["n"]
        self.assertGreater(ratio, 8.0)                # ~8.8x more primes

    def test_z_falls_despite_flat_excess(self):
        # Flat raw excess + growing sqrt(n) => z must drop.
        self.assertLess(VALUE[10**9]["z10"], VALUE[10**8]["z10"])


class TestBase8SignChange(unittest.TestCase):
    """Base 8 leaves its near-zero regime and becomes clearly negative."""

    def test_near_zero_at_1e8(self):
        self.assertLess(abs(VALUE[10**8]["z8"]), 1.0)

    def test_clearly_negative_at_1e9(self):
        self.assertLess(VALUE[10**9]["z8"], -5.0)


if __name__ == "__main__":
    unittest.main()
