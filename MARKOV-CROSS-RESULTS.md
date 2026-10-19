# Cross-sequence Markov-order of gap digit-sum parity: graded result

*Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
Pre-registered in [`notes/markov-order-cross-prereg.md`](notes/markov-order-cross-prereg.md)
before the compute job `markov-cross-order-v1` was run on 2026-10-19. The four
predictions and thresholds below were committed in advance and are graded here exactly as
written; the job self-grades and the grading was not edited after the run.*

**Honest status up front:** $0 earned on this arc, no external engagement (0 stars/forks,
0 reactions). This is a curiosity/reproducibility effort, stated plainly. A negative or
partial result is reported as such.

## The question

Run 97 of this repo (`MARKOV-ORDER-RESULTS.md`, v0.14.0) established that the base-*b*
digit-sum **parity** of consecutive **prime** gaps is genuinely **higher-order Markov**,
not first-order: the triple conditional-independence test

> **H0 (first-order Markov):** `s_{n+1}` is independent of `s_{n-1}` given `s_n`
> (a 2×2×2 table, likelihood-ratio *G²*, df = 2, exact `p = exp(-G²/2)`)

is rejected in all six bases {2, 4, 8, 10, 12, 16}. That test had **only ever been run on
primes.** This job runs the identical test on three sieves to `N = 100,000,000`:

- primes (A000040) — 5,761,455 terms, re-run as an internal replication/anchor,
- squarefree numbers (A005117) — 60,792,694 terms,
- sums of two squares (A001481 from 1) — 18,457,847 terms.

It asks: is the higher-order memory **universal across sieves**, and is **base 2** always
the loudest base, or is the order's base profile sequence-specific?

## Graded outcome — all four pre-registered predictions HOLD

| Prediction | Statement (threshold frozen in prereg) | Result |
|---|---|---|
| **P-MC-1** | every sieve: base-2 **global** G² ≥ 100 | **holds** — primes 2309, squarefree 57567, sums-of-two-squares 49133 |
| **P-MC-2** | every sieve: base-2 **blocked** σ-over-df ≥ 5 (100-block drift control) | **holds** — primes 117, squarefree 2878, sums-of-two-squares 2461 |
| **P-MC-3** | every sieve: global G² ≥ 40 in **all six** bases | **holds** (every base rejects in every sieve) |
| **P-MC-4** | the loudest base is **not** the same across all three sieves | **holds** — argmax base = **2** (primes), **4** (squarefree), **2** (sums of two squares) |

Global G² by base (larger ⇒ stronger departure from first-order Markov):

| base | primes | squarefree | sums of two squares |
|---:|---:|---:|---:|
| 2  | **2309** | 57567 | **49133** |
| 4  | 161 | **13,386,359** | 8487 |
| 8  | 495 | 10,654,974 | 5470 |
| 10 | 165 | 10,653,972 | 12,554 |
| 12 | 737 | 10,653,983 | 17,759 |
| 16 | 285 | 10,653,983 | 29,641 |

Interpretation, per the prereg's stated grid: all four holding is the strongest
"**not prime-specific**" conclusion of the arc — the higher-order Markov structure is a
generic property of structured integer-gap parity sequences, carrying a **sieve-specific
base fingerprint** rather than a universal preferred base. Note that for squarefree the
loudest base is **4**, not 2, and base-2 is comparatively modest (57567 vs 13.4 million at
base 4); base-2 dominance is a prime (and sums-of-two-squares) feature, not universal.

## What is genuinely new here, and what is not

This write-up is deliberately conservative about novelty (runs 106–112 hand-derived the
underlying correlator theory; see below).

- **P-MC-1 / P-MC-2 are weak-novelty.** The raw base-*b* digit-sum parity sequence of the
  integers already carries a full, non-first-order automatic correlation structure: it is
  (for base 2) the Thue–Morse sequence, whose *k*-point correlations are known in closed
  form (Baake & Coons, arXiv:2209.07102), with the exact lag-1 autocorrelation of the
  base-*b* parity sequence equal to **−(b−1)/(b+1)** for even *b* and **−1** for odd *b*
  (hand-derived run 112; a special case of Gelfond 1968 / Allouche–Shallit). So "the
  parity sequence is more than first-order Markov" is expected, not a discovery; our test
  being applied to *gaps* of three sieves adds empirical breadth but the direction is not
  surprising.
- **P-MC-3 / P-MC-4 carry the new content.** The *uniformity* across three independent
  sieves at N = 10⁸ (every base rejecting in every sieve) and, above all, the
  **sieve-specific loudest base** (base 2 for primes and sums of two squares, base 4 for
  squarefree) are the genuinely informative, first-of-their-kind measurements here. An
  October 2026 search of the literature and GitHub found no prior cross-sieve comparison
  of the Markov order of gap digit-sum parity.

## Degeneracy caveat (disclosed, not a bug)

Every squarefree gap below 10⁸ is ≤ 10 (`max_gap = 10`). For base *b* ≥ 11 a single-digit
gap has digit sum equal to the gap itself, so bases 12 and 16 see **identical**
triple-count tables (G² = 10,653,983 for both) and base 10 is almost identical
(10,653,972). This is a real property of a small-gap sequence, not a measurement error.
P-MC-3's "all six bases reject" is still satisfied, and P-MC-4's argmax (base 4) is
unaffected, but the six squarefree bases are effectively **four** distinct tests. Primes
(max gap 220) and sums of two squares (max gap 60) have no such collapse.

## Prior art

Cross-sequence / digit-sum-of-primes background: Gelfond 1968 (equidistribution of digit
sums in residue classes); Mauduit & Rivat, *Ann. of Math.* 171 (2010) 1591 (sum of digits
of primes); Dietmann, Elsholtz & Shparlinski (arXiv:1601.04754); arXiv:2002.05675;
arXiv:2108.11293. Automatic-sequence correlation theory used to scope novelty: Baake &
Coons (arXiv:2209.07102, *Indag. Math.*) — all odd-order balanced Thue–Morse correlations
vanish and the even-order ones are given in closed form; Allouche & Shallit,
*Automatic Sequences*; Drmota & Skalba, and Dumont, on the Newman phenomenon. None of
these runs the triple conditional-independence (Anderson–Goodman 1957; Billingsley 1961)
Markov-order test across multiple integer sieves' gaps, which is what is measured here.

## Reproduce

The committed data is `data/markov_cross.json` (full per-sieve, per-base triple counts,
global and blocked G², lag-1/lag-2 autocorrelations, and the self-test log). The job's own
self-tests (prime counts, squarefree/sums-of-two-squares reference agreement, numpy↔pure-
Python cross-checks, χ² closed forms, and G-test gates on synthetic independent/dependent
tables) all passed before any number was emitted. `test_markov_cross_results.py`
independently **re-derives** every base's global G² from the committed triple counts in
pure standard-library Python (no numpy, no network) and asserts each graded prediction.

```
python -m unittest -v test_markov_cross_results.py
```
