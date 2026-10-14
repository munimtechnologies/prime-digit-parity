"""Tests for stats.py (standard library only)."""

import unittest

import stats


class TestStats(unittest.TestCase):
    def test_excess(self):
        self.assertEqual(stats.excess(10, 4), 6)
        self.assertEqual(stats.excess(4, 10), -6)
        self.assertEqual(stats.excess(7, 7), 0)

    def test_fair_coin_std(self):
        self.assertAlmostEqual(stats.fair_coin_std(50_000_000),
                               7071.067811865476, places=4)
        self.assertAlmostEqual(stats.fair_coin_std(100), 10.0, places=9)

    def test_parity_z_balanced(self):
        self.assertEqual(stats.parity_z(5, 5), 0.0)
        self.assertEqual(stats.parity_z(0, 0), 0.0)

    def test_parity_z_sign(self):
        self.assertGreater(stats.parity_z(60, 40), 0)
        self.assertLess(stats.parity_z(40, 60), 0)

    def test_parity_z_unit(self):
        # 2 odd, 0 even: excess 2, n 2, z = 2/sqrt(2) = sqrt(2)
        self.assertAlmostEqual(stats.parity_z(2, 0), 2 ** 0.5, places=9)

    def test_gupta_sridhar_exact(self):
        n, odd, even, exc, z = stats.gupta_sridhar()
        self.assertEqual(n, 50_000_000)
        self.assertEqual(odd, 25_032_384)
        self.assertEqual(even, 24_967_616)
        self.assertEqual(exc, 64_768)
        self.assertAlmostEqual(z, 9.16, places=2)


if __name__ == "__main__":
    unittest.main()
