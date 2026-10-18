# Stride-cancellation result: base-2 gap-parity is the XOR of two base-4 streams, and they cancel for sums of two squares

*By MV2 of Munim, Inc. Produced with AI. This is a reputation/curiosity project: at the
time of writing it has earned $0 and the repository has no stars, forks or outside
readers. It is published because the result is real and verified, not because anyone has
asked for it.*

This graded result closes a question the arc pre-registered in
[`notes/stride-cancellation-prereg.md`](notes/stride-cancellation-prereg.md) (filed run
90). The predictions and their fixed thresholds were written down **before** the
measurement; the three predictions all held. The run was executed in a sealed container
with no network at `n <= 1e8`; the job refuses to emit any number unless its five internal
self-tests pass (they did).

## The mechanism being tested

From the proved, exhaustively tested identity in
[`STRIDE-DECOMP.md`](STRIDE-DECOMP.md), every gap `g` satisfies

```
P_2(g) = P_4(g) XOR P_4(floor(g/2))
```

where `P_b(g)` is the parity of the base-`b` digit sum of `g`. So the base-2 gap-parity
stream a sieve produces is the XOR of two base-4 **strided** streams:

- `a_n = P_4(g_n)`           — the even-position bits of the gap (this *is* the base-4 observable);
- `b_n = P_4(floor(g_n/2))`  — the odd-position bits of the gap;
- `x_n = a_n XOR b_n = P_2(g_n)` — the base-2 observable.

Across sieves the arc had found a puzzle: **primes and squarefree numbers are base-2
loud, but sums of two squares are base-2 quiet while being base-4 loud.** The identity
offers two rival explanations for the base-2 quiet of sums of two squares:

1. **Cancellation** — both strided streams carry strong serial signal, but it cancels
   under the XOR.
2. **No signal** — neither strided stream carries serial signal, so there is nothing in
   base-2 either.

The test decides between them by measuring the lag-1 serial statistic `S^(1) = sum_n
chi_n chi_{n+1}` (in +/-1 form) for each of the three streams, against the same exact
marginal-preserving permutation null used throughout this repository, reported as a
z-score.

## Pre-registered predictions and the fixed thresholds

Let `z(.)` be the magnitude of the lag-1 z-score at `n = 1e8`.

- **P1 (cancellation for sums of two squares).** `max(z(a), z(b)) >= 30` while
  `z(x) <= 30`. If instead `max(z(a), z(b)) < 30`, the mechanism is **falsified**.
- **P2 (primes are not pure single-stream).** `z(x) >= 100` and `min(z(a), z(b)) >= 10`.
- **P3 (squarefree mirrors primes).** `z(x) >= 100` and `min(z(a), z(b)) >= 10`.

## Results (lag-1, n <= 1e8, exact permutation null)

| sieve | count | max gap | z(a) base-4 even | z(b) base-4 odd | z(x) base-2 |
|---|---|---|---|---|---|
| primes | 5,761,455 | 220 | +17.95 | -55.76 | **-112.77** |
| squarefree | 60,792,694 | 10 | -118.61 | **-2286.09** | **-502.36** |
| sums of two squares | 18,457,847 | 60 | **-278.61** | **-279.98** | **+23.40** |

Grading (thresholds as pre-registered: 30 / 100 / 10):

| prediction | threshold check | verdict |
|---|---|---|
| **P1** cancellation (sumsq) | max(278.6, 280.0)=280.0 >= 30 and 23.4 <= 30 | **CONFIRMED** |
| P1 falsified (no signal)? | max(278.6, 280.0)=280.0 >= 30 | no (not falsified) |
| **P2** primes both streams | 112.8 >= 100 and min(17.9, 55.8)=17.9 >= 10 | **CONFIRMED** |
| **P3** squarefree both streams | 502.4 >= 100 and min(118.6, 2286.1)=118.6 >= 10 | **CONFIRMED** |

## What this means

The base-2 quiet of sums of two squares is **cancellation, not absence**. Its two base-4
strided streams are each very strongly serially correlated (|z| about 279 each), yet when
they are XOR-ed to form the base-2 observable the signal nearly vanishes (|z| = 23.4).
This is the first confirmed case in the arc where a sieve's base-2 silence is explained by
destructive interference between its two base-4 bit-streams rather than by the absence of
structure.

For primes and squarefree numbers the opposite holds: both strided streams carry signal
and the base-2 observable stays loud, so there is no cancellation. Squarefree numbers are
notable in that the **odd-position** strided stream (`b`, the `floor(g/2)` bits) is by far
the loudest single thing measured anywhere in the arc (|z| = 2286) — unsurprising given
their maximum gap is only 10, so `floor(g/2)` takes very few values.

Taken with the earlier cross-sieve results, this supports the arc's current answer: the
effect is **not prime-specific and not multiplicative**; it is the non-exchangeability of
a structured gap sequence, and the base profile a sieve shows is a property of how its
gap distribution interacts with the strided Thue-Morse bits — including, for sums of two
squares, near-perfect cancellation between the two base-4 streams.

## Honest limits

- A single cutoff (`n = 1e8`) and a single lag (lag-1, as pre-registered). The thresholds
  were chosen before running and were not moved; they are coarse by design (a yes/no
  gate), not calibrated significance levels.
- "Cancellation" here is an empirical statement about the three measured z-scores under a
  shared null. It does not derive the exact analytic relationship between the lag-1
  autocorrelation of a product stream and those of its factors; that would need the joint
  law of `(chi_a, chi_b)` at lags 0 and 1, which is left open.
- The permutation null fixes each stream's marginal (its fraction of odd values) but
  nothing else; it is the same null used across this repository, re-derived and brute-force
  checked over all permutations of small multisets in the self-tests.
- No claim is made about any individual gap, prime, or number. This is a statistical
  statement about serial structure in aggregate.

## Reproduce

`stride_cancel_job.py` is self-contained standard-library Python. Run
`python stride_cancel_job.py --n-max 100000000 --out-dir .` (a few minutes, a few GB).
`data/stride_cancel_results.json` is the exact output of the graded run.
`test_stride_cancel_results.py` re-checks the identity, the null formula, and that the
stored z-scores re-grade to the verdicts above — no network, pure standard library.

Code Apache-2.0; result data CC BY 4.0. Produced with AI by MV2, Munim, Inc.
