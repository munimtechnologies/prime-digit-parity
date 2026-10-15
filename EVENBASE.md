# A uniform even-base reduction for digit-sum parity

*Work by Munim, Inc., produced with AI. Author: MV2 of Munim, Inc. Code
Apache-2.0; this text and the data table CC BY 4.0.*

This note extends the project's "tower" of parity representations, which so far
covered only the **powers of two** (`REDUCE.md`, `CORRELATION.md`,
`HIERARCHY.md`), to **every even base**, including the non-powers-of-two 6, 10
and 12. It gives one elementary identity that writes the digit-sum parity of any
even base as a product of single-bit characters of the number's successive
base-`b` truncations, and uses it to answer — in the negative — one natural
guess about why the cross-base prime biases are ordered the way `RESULTS.md`
found them.

## 1. The identity

Let `s_b(n)` be the parity (mod 2) of the sum of the base-`b` digits of
`n >= 0`, and let `chi_b(n) = (-1)^{s_b(n)} in {+1,-1}`. Write `eps(x) = (-1)^x`
for the least-significant-bit character.

**Theorem (uniform even-base reduction).** For every **even** base `b` and every
`n >= 0`,

```
s_b(n)  ==  sum_{i>=0} floor(n / b^i)    (mod 2),                      (E)
```

equivalently

```
chi_b(n)  =  prod_{i>=0} eps( floor(n / b^i) ).                        (E')
```

**Proof.** Legendre's classical formula for the base-`b` digit sum `S_b(n)` is

```
S_b(n) = n - (b-1) * sum_{i>=1} floor(n / b^i).
```

Reduce mod 2. Because `b` is even, `b-1` is odd, so `(b-1)*X == X (mod 2)` for
every integer `X`. Hence

```
S_b(n) == n + sum_{i>=1} floor(n/b^i) == sum_{i>=0} floor(n/b^i)   (mod 2),
```

and `s_b(n) = S_b(n) mod 2`, which is (E). Exponentiating `(-1)^{(.)}` turns the
sum into the product (E'). A second, self-contained proof: for even `b`,
`(x mod b) mod 2 = x mod 2`, so the parity of the `i`-th base-`b` digit
`d_i = floor(n/b^i) mod b` equals `floor(n/b^i) mod 2`; summing over `i` gives
(E) directly. ∎

The product in (E') is finite: the factor for index `i` is `+1` once
`b^i > n`, so there are exactly `floor(log_b n) + 1` non-trivial factors. The
`i = 0` factor is `eps(n)`, which is `-1` for **every odd `n`**, hence for every
prime `p > 2` in every even base.

Both forms are machine-checked exhaustively in `test_evenbase.py` for
`b in {2,4,6,8,10,12,16}` over the first tens of thousands of `n`, together with
the Legendre consistency step and the prime corollary of §3.

## 2. Relation to the dyadic tower (what is and isn't new)

The identity (E) is elementary — it is the mod-2 shadow of Legendre's formula,
so I do not claim the identity itself as new. What is new to this project, and
as far as prior-art search shows not stated elsewhere for prime digit-sum
parity, is the **uniform framing across all even bases** and its consequences:

* For `b = 2^a` the `i`-th truncation's lsb is exactly bit `a*i` of `n`
  (proved in `test_evenbase.py::test_dyadic_factors_are_single_bits`), so (E')
  reads `chi_{2^a}(n) = XOR of bits {0, a, 2a, ...}` — exactly the positional
  lemma of `HIERARCHY.md §1`. The dyadic tower is the special case in which the
  truncation factors are **disjoint single bits** of `n`.
* For a non-dyadic even base `b = 2^a * m` with `m > 1` odd (6 = 2·3,
  10 = 2·5, 12 = 4·3), the truncations `floor(n/b^i)` divide by the odd part `m`
  as well, so `eps(floor(n/b^i))` is **not** a single bit of `n`. These bases
  are therefore a structurally different family, outside the dyadic tower — yet
  (E') still represents them as a product of lsb characters of the same shape.

So (E') unifies all even bases into one object — *a product of lsb characters of
base-`b` truncations* — and isolates the odd part `m` as the single feature that
separates the dyadic tower (studied in the earlier notes) from bases 6, 10, 12.

## 3. Consequence for primes

Because the `i = 0` factor is `eps(p) = -1` for every odd prime, (E') gives,
for `p > 2`,

```
chi_b(p) = - prod_{i>=1} eps( floor(p / b^i) ).
```

With the excess `E_b(x) = #{p<=x: s_b(p) odd} - #{p<=x: s_b(p) even}
= - sum_{p<=x} chi_b(p)`, this says the entire even-base prime bias is carried by
the **joint behaviour of the truncation-lsb factors** `eps(floor(p/b^i))`,
`i >= 1` — the parities of `p` with its last 1, 2, 3, ... base-`b` digits
stripped off. (Verified as an exact sign identity for odd `n` in
`test_evenbase.py::test_prime_corollary_sign`.) This is the exact non-dyadic
analogue of the base-2-from-base-8 picture in `REDUCE.md`: no single factor need
be biased for the product to be.

## 4. A natural guess, refuted by the shipped data

How many fluctuating factors does each base have? To `x = 10^8` the count of
non-trivial factors (`i >= 1`) is `floor(log_b 10^8)`:

| base b | odd part m | # factors (i>=1) to 1e8 | excess E_b at 1e8 |
|-------:|-----------:|------------------------:|------------------:|
| 2  | 1 | 26 | 70509 |
| 4  | 1 | 13 | 65633 |
| 6  | 3 | 10 | 57813 |
| 8  | 1 |  8 | **-1329** |
| 10 | 5 |  8 | **62627** |
| 12 | 3 |  7 | 32235 |
| 16 | 1 |  6 | 104603 |

(Excesses are the project's verified `trajectory` sandbox run, two independent
sieves agreeing to 1e8; see `RESULTS.md` / `data/`.)

A natural guess is that the magnitude of the bias tracks the number of factors.
The data **refutes it cleanly**: base 8 and base 10 have the *same* eight
non-trivial factors to 1e8, yet base 8 is essentially unbiased (`E = -1329`)
while base 10 is strongly biased (`E = 62627`); and base 16, with the *fewest*
factors (6), has the *largest* excess (104603). The number of factors does not
determine, or even order, the bias. What differs between base 8 and base 10 is
the odd part `m` (1 vs 5), i.e. the internal structure and correlation of the
factors, not their count. This is consistent with the project's standing
finding that the signal lives in the *joint correlation* of the parities, not in
any one of them, and now shows it holds across the dyadic/non-dyadic divide.

## 5. Pre-registered prediction (small future sandbox job)

`compete/evenbase-factors-v1` measures, for bases 6, 8 and 10, the marginal
balance of each truncation-lsb factor `eps(floor(p/b^i))` and the direct
excesses `E_b` at 1e7 and 1e8. Predictions, to be scored honestly next run:

* **P-even-0 (gate):** direct excesses reproduce the anchor
  `E_6=57813, E_8=-1329, E_10=62627` at 1e8 (and `E_6=3181, E_8=215, E_10=2967`
  at 1e7). A mismatch means a sieve bug — fix before trusting anything else.
* **P-even-1:** each individual factor `eps(floor(p/b^i))` is close to balanced
  over primes (|mean| small, square-root size), for all three bases — so, as in
  the dyadic case, no single factor explains base 10's large excess.
* **P-even-2:** base 8 and base 10 agree on P-even-1 (balanced marginals) yet
  differ by ~60000 in total excess, locating the entire difference in the
  **joint correlation** of the factors, which is governed by the odd part `m`.

A confirmation of P-even-1 with P-even-2 would make precise, for non-dyadic
bases, the same mechanism the dyadic tower already isolates: the bias is a
connected correlation of individually balanced parities.

## 6. Prior art

* **Legendre** — the digit-sum formula `S_b(n) = n - (b-1) Σ_{i>=1} floor(n/b^i)`
  is classical; (E) is its mod-2 reduction for even `b`, stated here as a
  representation, not claimed as new.
* **Mauduit–Rivat (2010, Annals)** — base-`b` digit-sum of primes equidistributes
  mod 2 for `gcd(2, b-1) = 1` (all even `b`), with a power-saving error. An
  asymptotic density statement only: no finite-size law, no cross-base ordering,
  no truncation-factor representation.
* **Gupta–Sridhar (arXiv:1812.11841)** — one base-10 cutoff; no even-base family.
* **This project's earlier notes** — `REDUCE.md`, `CORRELATION.md`,
  `HIERARCHY.md` cover only powers of two. No prior statement, here or
  elsewhere found, of the uniform even-base reduction, its use to place bases
  6/10/12 in the same framework as the dyadic tower, or the factor-count
  refutation using base 8 vs base 10.

**Honest limits.** (E) is elementary. The contribution is organisational and
comparative, plus a clean negative result; it does *not* explain *why* the
joint correlation is large along primes (the project's standing open problem),
and the factor-balance predictions of §5 await the sandbox.
