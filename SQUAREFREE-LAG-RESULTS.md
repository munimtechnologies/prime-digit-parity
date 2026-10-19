# Does the squarefree base-2 gap-parity signal oscillate with lag period 2? — graded result

*Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
Pre-registered in [`notes/squarefree-lag-period-prereg.md`](notes/squarefree-lag-period-prereg.md)
before the compute job `squarefree-lag-period-v1` was run on 2026-10-19. The predictions
and thresholds below were committed in advance and are graded here exactly as written.*

**Honest status up front:** $0 earned on this arc, no external engagement. This is a
curiosity/reproducibility effort. **The main pre-registered hypothesis is FALSIFIED** — a
genuine negative result, reported in full because a clean negative is as informative as a
positive one.

## The question

The graded stride-cancel-lag result (`STRIDE-CANCEL-LAG-RESULTS.md`) found that the base-2
digit-sum parity of consecutive **squarefree** gaps, measured against an exact
marginal-preserving permutation null, is loud at lag 1, nearly silent at lag 2, and loud
again at lag 3:

| lag | 1 | 2 | 3 |
|---|---:|---:|---:|
| squarefree base-2 XOR \|z\| | 502.4 | 7.1 | 495.8 |

That looked like it might be a **period-2 structure in lag** (odd lags loud, even lags
quiet). This job extended the lag set to 1–8 to test that directly, for three sieves to
`N = 100,000,000`, re-using the proven sieves and null from the stride-cancel-lag job.

## Graded outcome — the period-2 hypothesis is FALSIFIED

Squarefree base-2 XOR stream, |z| by lag:

| lag | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| \|z\| | 502.4 | **7.1** | 495.8 | 425.9 | 441.3 | 438.7 | 427.5 | 333.9 |

| Prediction (threshold frozen in prereg) | Result |
|---|---|
| **P-OSC-1 (MAIN)** median(odd \|z\|) ≥ 5 × median(even \|z\|) | **FALSIFIED** — ratio = 468.5 / 379.9 = **1.23**, far below 5 |
| **P-OSC-2** every odd lag \|z\| ≥ 30 | holds (502, 496, 441, 428) |
| **P-OSC-3** every even lag \|z\| ≤ 30 | **FALSIFIED** — lags 4, 6, 8 are loud (426, 439, 334); only lag 2 is quiet |
| **P-OSC-4** (exploratory) base-4 streams `a`,`b` not all-quiet on even lags | holds |

**The real structure is not a period-2 oscillation — it is that lag 2 is *uniquely*
quiet.** Every other lag 1–8 is loud (|z| ≥ 300); only lag 2 collapses to the null (|z| =
7.1). So the clean odd-loud/even-quiet picture suggested by lags 1–3 does not continue: it
was an artifact of only looking as far as lag 3. This is exactly the kind of overreach
that pre-registration is meant to catch, and it did.

## Why lag 2 is special (mechanism, consistent with the data)

The base-2 observable is `x = a XOR b`, where `a = P₄(g)` (even-position base-4 bits of the
gap) and `b = P₄(⌊g/2⌋)` (odd-position bits) — the stride decomposition from
`STRIDE-DECOMP.md`. In the committed data the two base-4 streams are loud at essentially
every lag (|z| in the hundreds to thousands, e.g. squarefree `b` at lag 1 = 2286). At lag
2 the two streams' contributions to the XOR nearly cancel (squarefree `x` lag-2 |z| = 7.1
while `a` lag-2 |z| = 1891 and `b` lag-2 |z| = 914), and this cancellation is specific to
lag 2 — it does not recur at lags 4, 6, 8. The lag-2 near-null is therefore a genuine but
**isolated** XOR cancellation, not a periodic one.

## Exploratory cross-sieve comparison (reported, not graded)

The same odd/even summary for the other two sieves' base-2 streams, to see whether even a
weak period-2 tendency is squarefree-specific:

- **primes:** odd median |z| = 1.39, even median |z| = 1.79 — both near null; base-2 prime
  gaps carry almost no lag structure beyond lag 1 (lag-1 |z| = 112.8, all other lags < 12).
- **sums of two squares:** odd median |z| = 31.7, even median |z| = 27.4 — comparable, no
  odd/even split (lag 1 = 23.4 is actually the *weakest* odd lag here).

No sieve shows the predicted period-2 pattern. The hypothesis fails across the board.

## Prior art

Lag/serial-correlation structure of digit-sum sequences and the Thue–Morse family:
Mahler (1927) on the Thue–Morse spectrum; Baake & Coons (arXiv:2209.07102) on Thue–Morse
*k*-point correlations (all odd-order balanced correlations vanish); arXiv:2511.06386 on
automatic-sequence autocorrelation. A run-105 hand-check established that **raw** Thue–Morse
prefers the **opposite** (even) lag parity to the one hypothesized here, so neither the
hypothesized pattern nor its failure is a trivial Thue–Morse artifact — this is a property
of the *squarefree-gap* parity sequence specifically. An October 2026 search of the
literature and GitHub found no prior description of a lag-period-2 (or isolated-lag-2)
serial structure in the digit-sum parity of squarefree-number gaps.

## Reproduce

The committed data is `data/squarefree_lag_period.json` (per-sieve, per-stream serial
statistic S, null mean, null variance and z at lags 1–8). The job's self-tests passed
before any number was emitted. `test_squarefree_lag_results.py` independently re-derives
every z from the committed `(S, E_null, var_null)` in pure standard-library Python and
asserts each graded prediction, including the falsifications.

```
python -m unittest -v test_squarefree_lag_results.py
```
