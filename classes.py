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
# produced with AI (autonomous research agents). Stdlib only.
"""Exact closed forms for digit-sum parity excess on residue classes.

Setup. For an integer base b >= 2 let s_b(n) be the sum of the base-b
digits of n and let f_b(n) = (-1)**s_b(n), which is +1 when the digit sum
is even and -1 when it is odd. For a cutoff X >= 0 the *digit-sum parity
excess* over [0, X) is

    D_b(X) = sum_{0 <= n < X} f_b(n)
           = (# n < X with even digit sum) - (# n < X with odd digit sum).

For EVEN b this module gives exact closed forms for D_b(X), for the excess
restricted to a single last-digit class n == a (mod b), and for the excess
over the integers coprime to b. It proves (by construction + brute-force
test) the sharp bound |E_b^coprime(X)| <= phi(b) at every cutoff.

Why even b. For odd b, s_b(n) == n (mod 2) because every power of b is odd,
so the parity is degenerate (it just tracks n). For even b every power
b^k with k >= 1 is even, so the parity is a genuine Thue-Morse-type signal.
"""


def digit_sum(n, b):
    """Sum of the base-b digits of n >= 0 (b >= 2)."""
    if b < 2:
        raise ValueError("base must be >= 2")
    if n < 0:
        raise ValueError("n must be >= 0")
    total = 0
    while n:
        total += n % b
        n //= b
    return total


def integer_excess_closed(b, X):
    """Exact D_b(X) for EVEN base b, in {-1, 0, +1}.

    Derivation. Split [0, X) into blocks [b*j, b*(j+1)). Over a full block
    the last digit runs 0..b-1; since b is even, sum_{r=0}^{b-1} (-1)**r = 0,
    so every completed block contributes 0. With X = b*Q + R the only
    surviving part is the partial block j = Q, r = 0..R-1, whose sum is
    (-1)**s_b(Q) * (R mod 2). Hence

        D_b(X) = (-1)**s_b(X // b) * ((X % b) % 2)  in  {-1, 0, +1}.
    """
    if b % 2 != 0:
        raise ValueError("integer_excess_closed is for even bases only")
    if X <= 0:
        return 0
    sign = 1 if digit_sum(X // b, b) % 2 == 0 else -1
    return sign * ((X % b) % 2)


def _count_in_class(b, a, X):
    """# of m >= 0 with b*m + a < X, for a single digit 0 <= a < b."""
    if a < 0 or a >= b:
        raise ValueError("a must be a single base-b digit (0 <= a < b)")
    if X <= a:
        return 0
    return (X - a + b - 1) // b


def class_excess_closed(b, a, X):
    """Exact excess over {n < X : n == a (mod b)} for EVEN b, single digit a.

    n = b*m + a gives s_b(n) = s_b(m) + a (no carry, since a < b), so
    f_b(n) = (-1)**a * f_b(m) and the class excess is

        E_b^(a)(X) = (-1)**a * D_b(M_a),   M_a = #{m : b*m + a < X}.

    Therefore |E_b^(a)(X)| <= 1 for every X.
    """
    if b % 2 != 0:
        raise ValueError("class_excess_closed is for even bases only")
    m = _count_in_class(b, a, X)
    sign = 1 if a % 2 == 0 else -1
    return sign * integer_excess_closed(b, m)


def euler_phi(n):
    """Euler's totient phi(n) for n >= 1."""
    if n < 1:
        raise ValueError("n must be >= 1")
    result = n
    p = 2
    m = n
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def reduced_residues(b):
    """Digits a in [1, b-1] with gcd(a, b) = 1."""
    from math import gcd
    return [a for a in range(1, b) if gcd(a, b) == 1]


def coprime_excess_closed(b, X):
    """Exact digit-sum parity excess over {n < X : gcd(n, b) = 1}, EVEN b.

    Summing class_excess_closed over the reduced residues a gives

        E_b^coprime(X) = sum_{gcd(a,b)=1} (-1)**a * D_b(M_a),

    a sum of phi(b) terms each in {-1, 0, +1}, so

        |E_b^coprime(X)| <= phi(b)   for every X.

    For even b every reduced residue a is odd, so the bound is attained:
    at X = b*t with t % b odd all M_a equal t and every term is +1 or -1
    together, giving |E_b^coprime(b*t)| = phi(b). (E.g. b=10, X=110: +4.)
    """
    if b % 2 != 0:
        raise ValueError("coprime_excess_closed is for even bases only")
    return sum(class_excess_closed(b, a, X) for a in reduced_residues(b))


# ---- Brute-force references (any base), used only to test the closed forms ----

def integer_excess_bruteforce(b, X):
    return sum(1 if digit_sum(n, b) % 2 == 0 else -1 for n in range(max(0, X)))


def class_excess_bruteforce(b, a, X):
    return sum(1 if digit_sum(n, b) % 2 == 0 else -1
               for n in range(max(0, X)) if n % b == a)


def coprime_excess_bruteforce(b, X):
    from math import gcd
    return sum(1 if digit_sum(n, b) % 2 == 0 else -1
               for n in range(max(0, X)) if gcd(n, b) == 1)
