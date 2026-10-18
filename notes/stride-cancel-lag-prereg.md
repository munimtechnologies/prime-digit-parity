# Pre-registration: is the base-2 stride-cancellation a lag-1 artefact, or does it persist at higher lags?

*Pre-registered by MV2 of Munim, Inc. Produced with AI. Written and committed before the
measurement is run. It builds directly on the graded result in `STRIDE-CANCEL-RESULTS.md`
(v0.13.0), which confirmed at lag-1 that the base-2 quiet of sums of two squares is
destructive XOR interference between two loud base-4 strided streams, not an absence of
signal.*

## What is already known (lag-1, graded)

At `n <= 1e8`, lag-1, exact permutation null, for the three streams
`a_n = P_4(g_n)`, `b_n = P_4(floor(g_n/2))`, `x_n = a_n XOR b_n = P_2(g_n)`:

| sieve | z(a) | z(b) | z(x) |
|---|---|---|---|
| sums of two squares | -278.6 | -280.0 | +23.4 |
| primes | +17.9 | -55.8 | -112.8 |
| squarefree | -118.6 | -2286.1 | -502.4 |

Cancellation (streams loud, XOR quiet) is confirmed for sums of two squares at lag-1.

## The new question

The confirmed cancellation is a lag-1 statement only. In +/-1 form,
`chi_x,n chi_x,n+1 = (chi_a,n chi_a,n+1)(chi_b,n chi_b,n+1)`, so the base-2 lag-1
autocorrelation is the mean over n of the *product* of the two streams' lag-1
agreement indicators. That this product is near-null at lag-1 says nothing about higher
lags. **Is the cancellation a general structural feature of the sums-of-two-squares gap
sequence, recurring at lag-2 and lag-3, or is it specific to adjacent gaps (lag-1)?**

## Pre-registered predictions (graded exactly as written, thresholds fixed now)

Let `z_k(.)` be the magnitude of the lag-k z-score at `n = 1e8`. Thresholds: a stream is
"loud" at `|z| >= 30`; the XOR is "quiet" at `|z| <= 30`; a base-2 signal "persists" at
`|z| >= 10`.

- **PRIMARY — P-L2 (cancellation persists for sums of two squares at lag-2).**
  `max(z_2(a), z_2(b)) >= 30` AND `z_2(x) <= 30`.
  - **Falsified (cancellation is lag-1-specific)** if `z_2(x) > 30`: the base-2 observable
    becomes loud at lag-2 even though it is quiet at lag-1.
  - **Inconclusive (no lag-2 stream signal)** if `max(z_2(a), z_2(b)) < 30` AND
    `z_2(x) <= 30`: there is no lag-2 stream signal to cancel, so cancellation is neither
    confirmed nor falsified at lag-2.

- **SECONDARY — P-L3 (same shape at lag-3).** `max(z_3(a), z_3(b)) >= 30` AND
  `z_3(x) <= 30`, with the same falsified / inconclusive branches.

- **SECONDARY — P-PERSIST (loud sieves stay loud in base-2 at lag-2).** For primes and
  for squarefree, `z_2(x) >= 10`. If either drops below 10, the base-2 serial signal does
  not survive to lag-2 for that sieve.

All z-scores for all three sieves, all three streams, at lags 1, 2, 3, are reported
regardless of outcome, including negatives and inconclusives.

## Method (fixed in advance)

Identical to `stride_cancel_job.py` but with `LAGS = [1, 2, 3]`. Sieve primes
(Eratosthenes), squarefree (segmented mu^2 sieve) and sums of two squares (p=3 mod 4 even
exponent) up to 1e8; form consecutive gaps; map each gap to the three stream bits via the
exact formulas; compute `S^(k) = sum_n chi_n chi_{n+k}` and convert to a z-score with the
repository's exact marginal-preserving permutation-null variance at lag k. Grade P-L2,
P-L3, P-PERSIST against the fixed thresholds.

## Self-test gate

The job re-runs the arc's self-tests before emitting any number and exits non-zero on any
failure: digit-sum parity values; the stride identity `P_2(g) = P_4(g) XOR P_4(floor(g/2))`
and `P_2(g) = popcount(g) mod 2` for all small g; `S^(k) = M - 2*D_k` checked directly at
k = 1, 2, 3 on a known string; the permutation-null mean and variance against brute-force
averaging over all permutations of small +/-1 multisets at k = 1, 2, 3; and prime /
squarefree / sums-of-two-squares prefix counts against known values. Expectations match the
finite-sample estimator exactly (no idealised limits), so a correct implementation passes.

## Honesty

Thresholds (30, 30, 10) are set before running and will not be moved. A null, mixed or
inconclusive result is a valid outcome and will be published as such. This is a
reputation/curiosity project: $0 earned, no stars. Nothing here has been measured yet.
