# Prior art and novelty: the serial gap-parity statistic

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. This note was
written to satisfy an obligation this project had not yet discharged in one
place: a careful search of the literature for the exact object we study, with a
precise statement of what is and is not new. Nothing here is a theorem; the
empirical results it refers to live in `RESULTS.md`, `GAP-CORR.md`,
`LAGK-RESULTS.md` and the companion `squarefree-gap-parity` repository. No money
has been earned and no external referee has checked this; treat it as a
reputation-and-curiosity contribution.*

## The statistic, stated precisely

Let `a_1 < a_2 < a_3 < ...` be a strictly increasing integer sequence produced
by some sieve (primes; squarefree numbers; sums of two squares; lucky numbers).
Form the consecutive gaps

    g_n = a_{n+1} - a_n.

Fix an integer base `b >= 2`. Let `s_b(g)` be the sum of the base-`b` digits of
`g`, and define the +/-1 character

    chi_b(g) = (-1)^(s_b(g)).

The object of study is the **lag-k serial correlation of `chi_b` along the
sequence of gaps**:

    r_k = ( 1 / (N-k) ) * sum_n chi_b(g_n) * chi_b(g_{n+k}),

standardised to a z-score under the null that the signs `chi_b(g_n)` are i.i.d.
fair coins. The empirical finding of the arc is that `r_1` is strongly and
robustly **negative** (anti-correlation), with a magnitude that depends on the
base but not on which of the four sieves produced the gaps.

Three exact facts about `chi_b` make the construction concrete. All three are
verified by direct computation in `test_prior_art.py`:

1. **Base 2 is Thue-Morse.** `s_2(g) mod 2` is the parity of the number of 1
   bits of `g`, so `chi_2(g) = (-1)^(t(g))` where `t` is the Thue-Morse sequence
   (OEIS [A010060](https://oeis.org/A010060)). The base-2 statistic is therefore
   literally the serial correlation of the Thue-Morse transform evaluated along
   consecutive gaps.
2. **Odd bases collapse.** For any odd `b`, `s_b(g) = g (mod 2)`, so `chi_b`
   reduces to the parity of the gap length and carries no base-specific
   information. This is why only even bases are studied (see `THEORY.md`).
3. **Powers of two are strided Thue-Morse.** For `b = 2^k`,
   `s_{2^k}(g) mod 2` equals the parity of the binary digits of `g` at
   positions `0, k, 2k, 3k, ...`. So comparing bases 2, 4, 8, 16 is comparing
   serial correlations of *strided subsamples* of the binary expansion -- a
   precise version of the informal "base 2 uses all bits, base 4 every other
   bit" framing used in the dev.to write-ups.

## What the literature contains, and how it differs

A search in October 2026 (Google Scholar, arXiv, OEIS, PyPI/GitHub) turned up
rich work on neighbouring objects but nothing on this one. The relevant strands,
and the precise gap in each:

- **Mauduit & Rivat, *Sur un probleme de Gelfond: la somme des chiffres des
  nombres premiers*, Annals of Math. 171 (2010), 1591-1646.** Equidistribution
  of `s_b(p) mod m` over primes `p`, with power-saving error. This is a *single*
  digit-sum, taken over the primes themselves, not over their gaps, and it is a
  marginal (first-order) statement, not a serial correlation.

- **Gupta & Sridhar, *A Peculiarity in the Parity of Primes*,
  [arXiv:1812.11841](https://arxiv.org/abs/1812.11841).** The base-10 excess of
  `s(p) mod 2`. Again a single value over the primes, not a serial statistic and
  not about gaps. (This is the anomaly `THEORY.md` starts from.)

- **Correlation of the Rudin-Shapiro sequence along primes,
  [arXiv:2306.10629](https://arxiv.org/abs/2306.10629); The Thue-Morse and
  Rudin-Shapiro sequences at primes in principal number fields,
  [arXiv:2001.07017](https://arxiv.org/abs/2001.07017).** These study
  correlation sums of a digit-parity sequence indexed by the primes themselves
  (the automatic sequence evaluated at `p`, or products over consecutive
  primes). The index set is the primes, not the *gaps between* them, and the
  transform is applied to `p`, not to `g_n`. Proving Thue-Morse is normal along
  primes is noted there as still open.

- **Baake & Coons, and *On the absolute value of the autocorrelations of the
  Thue-Morse sequence*, [arXiv:2511.06386](https://arxiv.org/abs/2511.06386).**
  Autocorrelations `eta(m)` of the *pure* Thue-Morse sequence `t(n)` over all
  `n`. This is the autocorrelation of `t` on its natural domain, never composed
  with a sieve's gap sequence.

- **Lemke Oliver & Soundararajan (2016) and the sum-of-two-squares analogue,
  [arXiv:2105.05048](https://arxiv.org/abs/2105.05048).** Correlations between
  the *last digits* (residue classes) of consecutive elements. This is a
  congruence statistic about the elements' residues, not a digit-sum-parity
  statistic, and not about the digit sums of the gaps.

- **Gap-distribution literature: Freiberg, Kurlberg & Rosenzweig, Commun.
  Number Theory Phys. 11 (2017) 837-877 (gaps between sums of two squares);
  Hooley, Acta Math. 127 (1971) 279-297 (intervals between sums of two
  squares).** These describe the *distribution* and *pair correlation* of gap
  *lengths* -- spacing/marginal statistics -- not the serial correlation of a
  digit transform of consecutive gaps.

## The novel gap

The object "apply a fixed digit-sum-parity transform `chi_b` to the consecutive
**gaps** of a sieve, then measure its **serial** correlation across consecutive
gaps" does not appear in any of the above. The two moves that distinguish it are
(i) transforming the *gaps* rather than the elements, and (ii) correlating
*serially* across consecutive terms rather than studying a marginal or a
pair-correlation of lengths. To the best of a thorough search, the construction
and its empirical signal are new. The one place it currently appears online is
this project's own dev.to posts.

## Honest limitations

- **Bounded gaps weaken base comparisons for squarefree numbers.** Below 1e8 the
  squarefree gaps are small (see `squarefree-gap-parity/RESULTS.md`); once the
  base exceeds the largest gap, `s_b(g) = g` and `chi_b` collapses to gap-length
  parity. Base comparisons there are only meaningful among small bases. Sums of
  two squares (gaps unbounded, comparable density) were added as the control for
  exactly this.
- **An i.i.d.-sign null is not a model of the gaps.** A strong `r_k` for many
  `k` does *not* by itself prove higher-order memory: a first-order Markov sign
  process already has `r_k = r_1^k` at every lag. Distinguishing genuine
  higher-order structure from that trivial decay is the job of the
  pre-registered Markov-order test (`notes/markov-order-prereg.md`), whose
  compute result is still pending.
- **No proof, no external validation.** These are empirical serial
  correlations. Nothing here claims a mechanism or a theorem.

## Reproducing the identities

    python -m unittest -v test_prior_art.py

verifies, by direct computation, the Thue-Morse identity (and its agreement with
OEIS A010060), the odd-base collapse, and the power-of-two strided-subsample
identity used above.

-- MV2, Munim, Inc.
