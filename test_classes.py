# Copyright 2026 Munim, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Part of munimtechnologies/prime-digit-parity. Work by Munim, Inc.,
# produced with AI. Stdlib only (unittest).
"""Tests that the closed forms in classes.py match brute force, and that
the sharp coprime-excess bound |E_b^coprime(X)| <= phi(b) holds (and is
attained for base 10, improving the earlier |E| <= 8 of Fact 6 to the
tight constant phi(10) = 4)."""
import unittest

from classes import (
    digit_sum,
    integer_excess_closed,
    class_excess_closed,
    coprime_excess_closed,
    integer_excess_bruteforce,
    class_excess_bruteforce,
    coprime_excess_bruteforce,
    euler_phi,
    reduced_residues,
)

EVEN_BASES = [2, 4, 6, 8, 10, 12, 16]
XMAX = 1200


class TestIntegerExcess(unittest.TestCase):
    def test_closed_matches_bruteforce(self):
        for b in EVEN_BASES:
            for X in range(0, XMAX):
                self.assertEqual(
                    integer_excess_closed(b, X),
                    integer_excess_bruteforce(b, X),
                    msg="integer excess b=%d X=%d" % (b, X),
                )

    def test_bounded_by_one(self):
        for b in EVEN_BASES:
            for X in range(0, XMAX):
                self.assertIn(integer_excess_closed(b, X), (-1, 0, 1))


class TestClassExcess(unittest.TestCase):
    def test_closed_matches_bruteforce(self):
        for b in EVEN_BASES:
            for a in range(b):
                for X in range(0, 400):
                    self.assertEqual(
                        class_excess_closed(b, a, X),
                        class_excess_bruteforce(b, a, X),
                        msg="class excess b=%d a=%d X=%d" % (b, a, X),
                    )

    def test_each_class_bounded_by_one(self):
        for b in EVEN_BASES:
            for a in range(b):
                for X in range(0, 400):
                    self.assertIn(class_excess_closed(b, a, X), (-1, 0, 1))


class TestCoprimeExcess(unittest.TestCase):
    def test_closed_matches_bruteforce(self):
        for b in EVEN_BASES:
            for X in range(0, XMAX):
                self.assertEqual(
                    coprime_excess_closed(b, X),
                    coprime_excess_bruteforce(b, X),
                    msg="coprime excess b=%d X=%d" % (b, X),
                )

    def test_sharp_bound_phi(self):
        # Main theorem: |E_b^coprime(X)| <= phi(b) at every cutoff.
        for b in EVEN_BASES:
            phi = euler_phi(b)
            for X in range(0, XMAX):
                self.assertLessEqual(abs(coprime_excess_closed(b, X)), phi)

    def test_bound_attained(self):
        # phi(b) is the tight constant: attained at X = b*t with t % b odd.
        for b in EVEN_BASES:
            phi = euler_phi(b)
            attained = max(
                abs(coprime_excess_closed(b, X)) for X in range(0, XMAX)
            )
            self.assertEqual(attained, phi, msg="base %d" % b)

    def test_base10_specifics(self):
        # Fact 6 said |E| <= 8; the tight constant is phi(10) = 4.
        self.assertEqual(euler_phi(10), 4)
        self.assertEqual(coprime_excess_closed(10, 110), 4)
        self.assertLessEqual(
            max(abs(coprime_excess_closed(10, X)) for X in range(0, XMAX)), 4
        )


class TestReducedResiduesAllOdd(unittest.TestCase):
    def test_even_base_reduced_residues_are_odd(self):
        # Corollary generalizing Fact 4: for even b every reduced residue
        # (hence every prime > b) ends in an ODD digit.
        for b in EVEN_BASES:
            for a in reduced_residues(b):
                self.assertEqual(a % 2, 1, msg="base %d residue %d" % (b, a))


class TestDigitSum(unittest.TestCase):
    def test_small(self):
        self.assertEqual(digit_sum(0, 10), 0)
        self.assertEqual(digit_sum(19, 10), 10)
        self.assertEqual(digit_sum(255, 16), 30)


if __name__ == "__main__":
    unittest.main()
