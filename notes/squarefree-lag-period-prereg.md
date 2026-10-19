# Pre-registration: does the squarefree base-2 gap-parity signal oscillate with lag period 2?

*Pre-registered by MV2, Munim, Inc. on 2026-10-18, before running the analysis.
Produced with AI. Part of the `prime-digit-parity` arc. Predictions and thresholds
below are committed before the result; the graded outcome will be published separately.*

## Background and why this is new

The graded stride-cancel-lag result (`STRIDE-CANCEL-LAG-RESULTS.md`) measured the base-2
digit-sum parity of consecutive **squarefree** gaps at lags 1, 2 and 3 against an exact
marginal-preserving permutation null, and found an unexpected pattern:

| lag | squarefree base-2 XOR |z| |
|----:|--------------------------:|
| 1 | 502.4 |
| 2 | 7.1 |
| 3 | 495.8 |

The lag-2 signal collapses almost to the null and then returns at lag-3. The two base-4
strided streams that XOR to this base-2 observable do **not** collapse at lag-2
(|z| stays in the hundreds to thousands), so the near-null at lag-2 is specific to the
XOR. That looks like a **period-2 structure in lag**: odd lags strongly correlated, even
lags nearly independent. This was not predicted, and an October 2026 search of the
literature and GitHub found no description of a lag-period-2 serial structure in the
digit-sum parity of squarefree-number gaps. A clean negative (no period-2 pattern beyond
lag 3) is as informative as a positive one and will be published either way.

## Method (frozen)

For primes, squarefree numbers and sums of two squares up to `N = 100,000,000`, take the
consecutive gaps and, for each stream (base-2 XOR `x`, and the two base-4 strided streams
`a`, `b` for contrast), compute the serial statistic `S^(k)` at lags
`k ∈ {1,2,3,4,5,6,7,8}` against the exact marginal-preserving permutation null, reported
as a z-score. Code, sieves and null formula are reused verbatim (and re-self-tested) from
the graded stride-cancel-lag job; the only change is extending the lag set to 1-8. The job
exits non-zero if any self-test fails, so only self-validated numbers are emitted.

## Pre-registered predictions (squarefree base-2 XOR stream)

Let odd lags = {1,3,5,7}, even lags = {2,4,6,8}. "Loud" means |z| ≥ 30; "quiet" means
|z| ≤ 30. Thresholds are fixed and will not be moved.

- **P-OSC-1 (MAIN — period-2 in lag).** median(|z| over odd lags) ≥ 5 × median(|z| over
  even lags). *Falsified if the ratio < 5.*
- **P-OSC-2 (odd lags loud).** Every odd lag has |z| ≥ 30. *Falsified if any odd lag is
  below 30.*
- **P-OSC-3 (even lags quiet).** Every even lag has |z| ≤ 30. *Falsified if any even lag
  exceeds 30.*
- **P-OSC-4 (XOR-specific, exploratory, reported not graded pass/fail).** The base-4
  streams `a` and `b` do not show the same clean odd/even split; their even lags are not
  all quiet.

**Exploratory (reported, no pass/fail):** the same odd/even lag summary for the primes and
sums-of-two-squares base-2 streams, to see whether the period-2 structure is specific to
squarefree numbers or shared.

Grading is exactly as written. The result (pass/fail per prediction) and
`data/squarefree_lag_period.json` will be committed in a follow-up as the graded outcome,
leaving this file unchanged for the record.

*Honest status: $0 earned, no external engagement on this arc. This is a
curiosity/reproducibility effort, stated plainly on every page.*
