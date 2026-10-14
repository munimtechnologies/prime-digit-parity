# The base-10 prime digit-sum parity "peculiarity": what theory settles, and what it does not

*Munim, Inc. — research note. Produced with AI. Apache-2.0 (code) / CC BY 4.0 (text).*

**Status: partial. No numeric claim about primes is made here until the sandbox
returns it** (job `prime-digit-parity-v3`, primes to 10^9 and 10^10).

**Summary.** Gupta & Sridhar (arXiv:1812.11841, "A Peculiarity in the Parity of
Primes") report that the base-10 digit sum of primes is odd more often than even,
"about an order of magnitude farther from the mean than expected" over the first
50,000,000 primes, and know of no explanation. This note does three things:

1. It shows the *limit* statement is not mysterious: the parity classes have
   density 1/2 by the Mauduit–Rivat theorem, so any bias must be a lower-order
   (finite-sample) effect.
2. It proves elementary bounds (Facts 1–6) that rule out the ambient integers as a
   source of the bias and localize it to the high-order digits of the primes.
3. It does **not** claim the finite-sample bias is explained. Mauduit–Rivat gives
   only an upper bound on the error term, far too weak to predict a bias of a few
   standard deviations at 5·10^7 primes. Whether the bias is a one-signed
   secondary term (a "Chebyshev-type" lean) or fair-coin-like noise is exactly what
   the queued computation measures.

## 1. The claim

Gupta & Sridhar (MIT), *A Peculiarity in the Parity of Primes*, arXiv:1812.11841,
December 2018, posted under **General Mathematics (math.GM)** with no journal
reference. Over the **first 50,000,000 primes** (a prime-count cutoff, not a value
cutoff) the base-10 digit sum is significantly more often odd than even, measured
against a fair-coin standard deviation. They report that the effect persists across
bases (sometimes flipping which parity leads) and that the explanations they tried
(primes modulo random numbers, a Chebyshev-bias adjustment) did not satisfy them.

## 2. What the theory settles: density 1/2

Mauduit & Rivat, *Sur un problème de Gelfond: la somme des chiffres des nombres
premiers*, **Ann. of Math. 171 (2010) 1591–1646**, answering a 1968 question of
Gelfond, proved that for a base `q` and modulus `m` with `gcd(m, q−1) = 1` the
digit-sum `s_q(p)` over primes `p` is equidistributed in residue classes mod `m`.
Digit-sum **parity** is `m = 2`, so the hypothesis is `q−1` odd, i.e. **`q` even**.
This matches our Fact 5: by casting out `(q−1)`'s, `s_q(n) ≡ n (mod q−1)`, so for odd
`q` parity is just `n`'s parity and is degenerate (Fact 1). **Base 10 has `q−1 = 9`,
so the theorem applies:**

```
#{ p ≤ x : s_10(p) even }  =  π(x)/2  +  o(π(x)).
```

So the two classes each have density exactly 1/2: there is no *limiting* bias, and
the finite-sample excess `E(x) = #odd − #even` is a lower-order term. (The
Mauduit–Rivat proof gives a power saving, `E(x) ≪ x^{1−σ}` with a small `σ > 0`; we
use only the statement `E(x) = o(π(x))`.)

## 3. What the theory does NOT settle

- An upper bound `x^{1−σ}` is compatible with `z = E/√π(x)` being anywhere from
  O(1) to a power of `x`. It does not say `z` is large, small, one-signed, or
  growing. **It does not explain a bias of a few sigma at 5·10^7 primes.** An
  earlier draft of this note over-read the theorem as doing so; that was wrong and
  is corrected here.
- The comparison to a fair-coin `√n` is a modelling choice, not a theorem. The
  digit sum of a prime is not an independent fair coin: the last digit is forced
  odd (Fact 4), and the high digits follow Benford-like and prime-counting
  structure. A lean of order `√n · c` with small `c` would be perfectly consistent
  with Mauduit–Rivat (it is `o(n)`), and would be a genuine secondary term of the
  same flavour as Chebyshev's bias in prime races (e.g. base 2, where odious primes
  lead; OEIS A130911, verified to 10^19). The Drmota–Mauduit–Rivat central limit
  theorem for `s_q(p)` is too coarse to resolve a single residue bit.
- We did **not** find a published note connecting the Gupta–Sridhar "peculiarity"
  to Mauduit–Rivat or to a secondary-term analysis. The individual ingredients are
  standard; the controls and elementary bounds below are the contribution.

## 4. What is proved (stdlib-tested in `core.py` / `test_core.py`)

- **Facts 1 / 5** — odd bases are degenerate (`s_b(n) ≡ n mod b−1`); only even bases
  carry a signal; base 10 is non-degenerate.
- **Facts 2 / 3** — the *integer* excess over `[0, X)` is exactly 0 at `X = b^k` and
  lies in {−1, 0, +1} at every cutoff, so the count-vs-value cutoff convention
  cannot inflate a bias through the ambient integers (this corrected v0.1.0).
- **Fact 4** — every prime `p > 5` ends in an odd digit, so the base-10 prime bias
  equals (minus) the bias of the "beheaded" primes `⌊p/10⌋`.
- **Fact 6 (new in v0.1.2)** — the integers coprime to 10 (last digit 1, 3, 7, 9),
  a superset of the primes `> 5`, have digit-sum parity excess `|E| ≤ 8` at **every**
  cutoff, with an explicit closed form. So the control set is balanced to within a
  constant, and any prime excess that grows like `√n` is a property of the primes
  *inside* that balanced set. This is the clean baseline the fair-coin comparison
  lacked: the right question is "do primes deviate from the coprime-to-10 integers",
  and that deviation is the whole excess.

## 5. The measurement (queued)

`fast.py` / `job.py` stream primes with a segmented sieve and record, for bases
`{2, 4, 6, 8, 10, 16}`, the odd−even excess and `z = excess/√n` at value cutoffs
`x = 10^k` (k up to 9, then up to 10), at prime-count cutoffs every `10^6` primes
(including `N = 5·10^7`), and at the end of every sieve segment (a trajectory).
Pre-registered readings (written before the data exist):

- **Fair-coin-like:** `|z|` stays O(1) and changes sign along the trajectory →
  the Gupta–Sridhar value is a single-cutoff fluctuation.
- **Secondary term:** `z` keeps one sign and grows roughly like `√(log)`-to-`√n`
  scale → a genuine prime-specific lean, to be compared against base 2 (control).
- The cross-base table shows which parity leads and how `z` scales with the base.

## References

- C. Mauduit, J. Rivat, *Sur un problème de Gelfond : la somme des chiffres des
  nombres premiers*, Ann. of Math. (2) 171 (2010) 1591–1646.
- M. Drmota, C. Mauduit, J. Rivat, central limit theorem for the sum-of-digits
  function of primes.
- A. D. Gupta, M. Sridhar, *A Peculiarity in the Parity of Primes*,
  arXiv:1812.11841 (math.GM, 2018).
- A. O. Gelfond, *Sur les nombres qui ont des propriétés additives et
  multiplicatives données*, Acta Arith. 13 (1968) 259–265.
- V. Shevelev, arXiv:0706.0786, arXiv:0707.1761 (odious-prime race); OEIS
  A130911.
