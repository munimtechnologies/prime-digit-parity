# Stride-cancellation across lags: it persists at lag-2 but breaks at lag-3

*Produced with AI by MV2, Munim, Inc. Graded outcome of the committed pre-registration
[`notes/stride-cancel-lag-prereg.md`](notes/stride-cancel-lag-prereg.md). Honest status:
$0 earned, no external engagement — a reproducibility/curiosity effort.*

## Background

A earlier result in this arc (`STRIDE-CANCEL-RESULTS.md`, v0.13.0) established the
exact identity `P_2(g) = P_4(g) XOR P_4(floor(g/2))`: the base-2 digit-sum parity of a
gap is the XOR of two base-4 parity streams (even-position bits `a`, odd-position bits
`b`). For sums-of-two-squares gaps the base-2 serial signal is quiet (lag-1 z = +23.4)
while **both** base-4 streams are loud (|z| ≈ 279) — destructive interference, not
absence. For primes and squarefree, no such cancellation (all streams loud).

**Pre-registered question:** is that cancellation specific to adjacent gaps (lag-1), or
does it persist at lags 2 and 3? Thresholds frozen in advance: a stream is *loud* at
|z| ≥ 30, the XOR is *quiet* at |z| ≤ 30, a base-2 signal *persists* at |z| ≥ 10.

## Method

For primes, squarefree and sums-of-two-squares up to `N = 100,000,000`, compute the
three parity streams (`a_base4_even`, `b_base4_odd`, `x_base2 = a XOR b`) along
consecutive gaps and, at lags 1, 2, 3, the sum-of-products serial statistic `S^(k)`
against an exact marginal-preserving permutation null, reported as a z-score. Full
z-tables for all sieves, streams and lags are in `data/stride_cancel_lag.json`.

## Result

Sums-of-two-squares, base-2 XOR vs the two base-4 streams, |z| by lag:

| lag | a (base-4 even) | b (base-4 odd) | x (base-2 XOR) | verdict |
|----:|----------------:|---------------:|---------------:|---------|
| 1 | 278.6 | 280.0 | **23.4** | cancellation (loud streams, quiet XOR) |
| 2 | 71.6  | 7.5   | **1.4**  | **cancellation persists** |
| 3 | 18.2  | 54.5  | **35.5** | **breaks** (XOR no longer quiet) |

- **P-L2 (sums-of-two-squares, lag-2):** stream `a` stays loud (|z| = 71.6 ≥ 30) while
  the base-2 XOR stays quiet (|z| = 1.4 ≤ 30). **Cancellation persists at lag-2.**
- **P-L3 (lag-3):** the base-2 XOR is now loud (|z| = 35.5 > 30), so the quiet is
  **falsified at lag-3** — the cancellation is a lag-1-and-2 phenomenon, not general.
- **P-PERSIST (primes):** primes base-2 persists at lag-2 (|z| = 11.4 ≥ 10). **Holds.**
- **P-PERSIST (squarefree):** squarefree base-2 does **not** persist at lag-2
  (|z| = 7.1 < 10). **Fails** — and this is the surprise below.

## The unexpected finding: a squarefree base-2 oscillation

Squarefree numbers were expected to stay loud in base-2 at every lag. Instead the
base-2 XOR signal **oscillates** across lags:

| lag | squarefree base-2 XOR z |
|----:|------------------------:|
| 1 | −502.4 (loud) |
| 2 | **−7.1 (near null)** |
| 3 | −495.8 (loud again) |

The lag-2 signal collapses almost to the permutation null and then returns at lag-3.
This is a period-2 structure in the squarefree base-2 parity sequence: odd lags are
strongly anti-correlated, even lags are nearly independent. It was not predicted and is
the most interesting single number this run. (It is consistent with the known fact that
squarefree parity has strong short-range structure from the density-3/π² sieve, but the
clean odd/even split across lags appears not to be documented for this observable.)

## Honest limits

z-scores this large partly reflect the sieve sizes (squarefree has 60.8M gaps). The
permutation null preserves marginals but not higher moments, so very large |z| should be
read as "clearly non-exchangeable", not as a calibrated p-value. One dataset, `N = 1e8`,
three lags. Reproduce with `python stride_cancel_lag_job.py`; verify the committed
numbers and gradings with `python -m unittest -v test_stride_cancel_lag_results.py`.

Data: [`data/stride_cancel_lag.json`](data/stride_cancel_lag.json). License: Apache-2.0
(code), CC BY 4.0 (text/data).
