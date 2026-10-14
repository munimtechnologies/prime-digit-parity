# On the digit-sum parity bias of primes: value cutoffs vs count cutoffs

*Munim, Inc. — research note. Produced with AI. Apache-2.0 (code) / CC BY 4.0 (text).*

## Background and prior art

For an integer `n` and base `b`, let `s_b(n)` be the sum of its base-`b` digits.
The parity `s_b(n) mod 2` of primes has a well-studied history:

- **Base 2** ("odious" = odd digit sum, "evil" = even). Mauduit & Rivat (*Annals
  of Math.* 171, 2010) proved the base-2 digit sums of primes are equidistributed
  mod 2 (density 1/2 each), resolving a 1968 conjecture of Gelfond. Shevelev
  (arXiv:0706.0786, 0707.1761) conjectured the odious primes stay in the majority
  for all `n > 6`; this race (OEIS A130911) is verified to **10^19**
  (B. Chaffin, Jun 2026; T. D. Noe to 10^13 earlier). So base 2 is thoroughly mapped.
- **Base 10**. Gupta & Sridhar, *A Peculiarity in the Parity of Primes*
  (arXiv:1812.11841), report that over the **first 50,000,000 primes** the digit
  sum is "significantly more likely to be odd than even," with the odd count
  about *an order of magnitude farther from the mean than expected* — a bias they
  say they have no explanation for, persisting across bases (with the leading
  parity possibly changing).
- Navarrete-Pedraza (OEIS A380192, 2025) shows prime digit sums are approximately
  normal and uniform mod `n`, consistent with Mauduit–Rivat.

**What appears to be missing**, and what this note targets: a clean separation of
the *sampling convention* from the *arithmetic*. Gupta & Sridhar measured a
**prime-count cutoff** ("first N primes"), not a **value cutoff** (`primes ≤ x`).
Those are different experiments, and for digit-sum parity the difference is not
cosmetic. We did not find a prior study that isolates this.

## Two facts we can prove directly (no computation)

**Fact 1 — odd bases are trivial.** If `b` is odd then `b ≡ 1 (mod 2)`, so for
`n = Σ d_i b^i` we have `n ≡ Σ d_i = s_b(n) (mod 2)`. Hence `s_b(p)` is *odd* for
every prime `p > 2`, in every odd base. Digit-sum parity is only a non-trivial
statistic in **even** bases. (A parity "bias" reported in an odd base is just the
statement that primes are odd.)

**Fact 2 — even-base ranges are exactly balanced.** For even `b` and any `k ≥ 1`,
the digits `0..b-1` split into exactly `b/2` even and `b/2` odd values, so the
units digit alone contributes an unbiased ±1 to the digit-sum parity,
independent of the other digits. Therefore the integers in `[0, b^k)` are split
**exactly 50/50** by digit-sum parity. (Generating function: the parity count is
the coefficient symmetry of `(b/2 + (b/2) x)^k`.)

Both facts are exercised by `test_core.py` (odd-base triviality on all primes
< 50000; exact balance on complete ranges for `b = 2,4,6,8,10,16`).

## The hypothesis this repo tests

Fact 2 says the **integer** digit-sum parity excess is exactly 0 at value cutoffs
`x = b^k`, but at an arbitrary cutoff it inherits a deterministic "incomplete top
decade" boundary term. A **prime-count cutoff** `N = 5·10^7` lands at a value
(`p_N ≈ 9.82·10^8`) deep inside a decade — precisely where that boundary term is
largest and uncontrolled.

> **Hypothesis H.** The "surprisingly large / unexplained" base-10 prime
> digit-sum parity bias of Gupta & Sridhar is dominated by the choice of a
> prime-count cutoff. Measured at value cutoffs `x = 10^k`, the normalized bias
> `z(x) = (odd − even)/√π(x)` is small/oscillating (consistent with
> Mauduit–Rivat), while at prime-count cutoffs it is inflated by the boundary
> term.

Honest either/or:
- If `z(10^k)` is small while count-cutoff `z` is large → the mystery is a
  **sampling artifact**, and the fix is to report value cutoffs.
- If `z(10^k)` is *also* large and one-signed → a genuine arithmetic bias in
  prime digit-sum parity survives at round cutoffs, which would be a real and
  more interesting positive result.

## Method

`fast.py` streams primes to `10^9` with a segmented sieve and computes, for bases
`{2,4,6,8,10,16}`, the odd-minus-even excess at value cutoffs `10^k` and at
prime-count cutoffs every `10^6` primes. Base 2 is a **control**: its excess is
known (positive, Shevelev/Chaffin), so matching the published sign validates the
pipeline. `job.py` reports `z = excess/√n`.

**Status: results pending.** This repository ships the proved facts, the tested
code, and a queued computation. The verified numbers and the verdict on
Hypothesis H will be added once the compute job returns. No numeric claim about
primes is asserted here before it is computed.

## References

- C. Mauduit, J. Rivat, *Sur un problème de Gelfond: la somme des chiffres des
  nombres premiers*, Ann. of Math. 171 (2010) 1591–1646.
- V. Shevelev, arXiv:0706.0786; arXiv:0707.1761.
- A. Gupta, A. Sridhar, *A Peculiarity in the Parity of Primes*, arXiv:1812.11841.
- E. Navarrete-Pedraza, OEIS A380192 (2025).
- OEIS A027697, A027699, A130911, A000069, A001969.
