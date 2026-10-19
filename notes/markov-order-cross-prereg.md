# Pre-registration — Cross-sequence Markov-order test of gap digit-sum parity

*Pre-registered by MV2, Munim, Inc. Committed before the compute job
`markov-cross-order-v1` was run. Produced with AI. Code Apache-2.0; text CC BY 4.0.*

## Background and the open question

This repository's run-97 result (v0.14.0, `MARKOV-ORDER-RESULTS.md`) established that the
base-*b* digit-sum **parity** of consecutive **prime** gaps is genuinely **higher-order
Markov**, not first-order: the triple conditional-independence test

> **H0 (first-order Markov):** `s_{n+1}` is independent of `s_{n-1}` given `s_n`

is rejected globally in all six bases {2, 4, 8, 10, 12, 16} (base-2 global *G²* ≈ 2309 on
π(10⁸)=5,761,455 primes), and the rejection survives a 100-block drift control. That test
has **only ever been run on primes.**

Earlier in the arc two things were already shown for *other* integer sequences:
1. the lag-1 serial memory itself is **not** prime-specific — it is present, and often
   stronger, in squarefree numbers and sums of two squares (and in lucky numbers); and
2. there is **no universal dominant base** — primes and lucky numbers (base-16-dominated)
   and sums of two squares (base-4-loud, base-2-quiet) each have a different base profile.

The **Markov order** has not been compared across sequences. This study asks: is the
higher-order (beyond-first-order) memory **universal across sieves**, and is **base 2**
always the loudest base for it, or is the order's base profile sequence-specific too?

## Data and method (locked)

- **N = 100,000,000.** Three sequences, all to N:
  - primes (A000040) — re-run as an internal replication / anchor,
  - squarefree numbers (A005117),
  - sums of two squares (A001481, taken from 1).
- Observable: for each consecutive gap `g_n`, the parity of the base-*b* digit sum of
  `g_n`, for *b* ∈ {2, 4, 8, 10, 12, 16}, giving one 0/1 sequence per (sequence, base).
- Test statistic: the likelihood-ratio **G²** of the triple conditional-independence test
  (Anderson–Goodman 1957; Billingsley 1961), df = 2, exact `p = exp(-G²/2)`.
- **Global** reading: one 2×2×2 table over the whole sequence.
- **Blocked** reading: the sequence is cut into **100** contiguous blocks, the per-block G²
  summed (df = 2×blocks), controlling for the slow drift of the gap distribution so a
  rejection cannot be a non-stationarity artifact. We report how many σ the blocked G²
  sits above its null mean (`sigma_over_df = (G² − df)/√(2·df)`).
- The numpy fast paths (sieves, triple counts, autocorrelation) are cross-checked inside
  the job's self-tests against the proven pure-Python `markov_order.py` reference and, for
  the two sieves, against proven stdlib references (the arc's `build_bad` and a brute
  squarefree check). The job exits non-zero if any self-test fails, so only self-validated
  numbers are emitted.

## Pre-registered predictions (graded pass/fail)

- **P-MC-1 — higher-order memory is universal in base 2.** For *every* sequence, the base-2
  **global** G² ≥ 100 (i.e. first-order Markov is rejected overwhelmingly, p ≪ 10⁻²⁰).
- **P-MC-2 — it survives the drift control.** For *every* sequence, the base-2 **blocked**
  reading has `sigma_over_df` ≥ 5 (≥ 5σ above the null mean of a true first-order chain).
- **P-MC-3 — every base rejects.** For *every* sequence, the global G² ≥ 40 in **all six**
  bases (as primes showed).
- **P-MC-4 — the base profile of the order is sequence-specific.** The base with the
  largest global G² is **not the same** across all three sequences (at least two distinct
  argmax bases). This mirrors the arc's "no universal dominant base" finding, now for the
  higher-order component rather than the lag-1 component.

Thresholds (locked): reject = global G² ≥ 40; base-2 "overwhelming" = global G² ≥ 100;
blocking survival = `sigma_over_df` ≥ 5. These are fixed here and must not be edited after
the job runs.

## Interpretation grid (stated in advance)

- All four hold → the higher-order Markov structure is a **generic property of structured
  integer-gap parity sequences**, not a prime phenomenon, with a sequence-specific base
  fingerprint. (Strongest "not prime-specific" conclusion of the arc.)
- P-MC-1/2 hold but P-MC-3 fails for some sequence → higher-order memory is universal but
  not present in every base; which bases carry it is itself informative.
- P-MC-1 fails for squarefree or sums of two squares → the higher-order memory **is**
  (partly) prime-specific after all — a genuine surprise worth its own write-up.
- P-MC-4 fails (base 2 loudest everywhere) → there **is** a universal preferred base for
  the order, contradicting the lag-1 finding — also a reportable result.

A negative or partial result is a real result and will be published as such.
