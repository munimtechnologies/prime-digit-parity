# Markov-order result: prime-gap digit-sum parity is higher-order, not first-order

*Produced with AI by MV2, Munim, Inc. Graded outcome of the pre-registration in
[`notes/markov-order-prereg.md`](notes/markov-order-prereg.md). Honest status: $0
earned, no external engagement on this arc — a reproducibility/curiosity effort.*

## The question

Earlier files in this repo showed the base-*b* digit-sum parity of consecutive
**prime gaps** carries serial memory: its lag-1 autocorrelation is strongly non-zero
(base-2 lag-1 z ≈ −120.6 at 1e8) and several higher lags are individually significant.
But a **first-order Markov** process already has non-zero correlation at every lag
(`r_k = r1^k`). So "many lags are significant" does not prove memory beyond one step.

**Open question (pre-registered):** is the serial memory fully explained by one-step
dependence, or is it genuinely higher-order? We apply the classical triple
conditional-independence / Markov-order test (Anderson–Goodman 1957; Billingsley 1961)
to this sequence. A literature and GitHub search in October 2026 found no prior
Markov-order test of prime-gap digit-sum parity; a negative result counts.

## Method (frozen before running)

Primes up to `N = 100,000,000` (π(N) = 5,761,455; 5,761,454 gaps; max gap 220).
For each base `b ∈ {2,4,8,10,12,16}`, parity bit `s_n = digitsum_b(g_n) mod 2`.
Null **H0:** `s_{n+1}` ⊥ `s_{n-1}` given `s_n`. Over the eight triples
`(s_{n-1}, s_n, s_{n+1})`, sum the likelihood-ratio `G²` across the middle state
(df = 2, survival `p = e^{−G²/2}`). A **blocked** reading (100 contiguous blocks,
df = 200) controls for the slow drift of the gap distribution. Companion: `r1`, `r2`
and the first-order prediction `r1²` (first-order requires `r2 = r1²`).

Correctness gates (built into the job, non-zero exit on failure): π(x) checks;
numpy paths cross-checked against a pure-Python reference; chi-square survival checked
against closed forms; a synthetic first-order chain must NOT be flagged (`G²` = 0.115,
gate < 40) and a synthetic second-order chain MUST be (`G²` = 42,570, gate > 100).
All gates passed.

## Result: all five predictions confirmed — memory is higher-order

| Base | global G² (df 2) | blocked G² (df 200) | r1 | r2 | r1² | (r2−r1²)/r1² |
|------|-----------------:|--------------------:|------:|------:|------:|-------------:|
| 2  | 2309.24 | 2531.10 | −0.04698 | +0.004763 | 2.207e−3 | **+1.158** |
| 4  |  160.90 |  358.73 | +0.00748 | −0.001503 | 5.59e−5  | −27.9 |
| 8  |  495.19 |  722.29 | −0.01375 | −0.004333 | 1.89e−4  | −23.9 |
| 10 |  165.18 |  376.23 | −0.00170 | −0.004635 | 2.88e−6  | −1609 |
| 12 |  737.10 |  975.52 | −0.00501 | −0.003299 | 2.51e−5  | −132.5 |
| 16 |  285.29 |  520.67 | −0.02327 | −0.005043 | 5.42e−4  | −10.3 |

- **P-M1 (main):** base-2 global `G²` = 2309.24 ≥ 30, `p ≈ 0`. First-order Markov is
  rejected. **Confirmed.**
- **P-M2 (robust):** base-2 blocked `G²` = 2531.10 (df 200), `p ≈ 0 ≤ 1e−3`. The
  rejection survives blocking, so it is not an artifact of the drifting gap
  distribution. **Confirmed.**
- **P-M3 (base-2 dominance):** base 2 has the largest global `G²` of the six bases.
  **Confirmed.**
- **P-M4 (first-order relation violated):** base-2 `|r2 − r1²| / |r1²|` = 1.158 > 0.10.
  **Confirmed.**
- **P-M5 (sanity):** gates pass; `G² ≥ 0`; `0 < p ≤ 1`; triple counts sum to
  `n_gaps − 2`. **Confirmed.**

**Beyond the pre-registration:** every one of the six bases rejects first-order Markov
(all global `G² ≥ 160`, all blocked `p` far below 1e−3), so the higher-order memory is
not a base-2 peculiarity — base 2 is merely the loudest.

## What this means, and honest limits

The serial dependence in prime-gap digit-sum parity is **genuinely second-order or
higher**: knowing `s_{n-1}` tells you something about `s_{n+1}` beyond what `s_n`
already carries. The effect is statistically overwhelming because `n` is large, but
the **effect size is small**: the one-step correlation is only `r1 ≈ −0.047` in base 2
and smaller in other bases. The result is a statement about structure, not a strong
predictor. This is a frequentist conditional-independence test on one dataset (primes
to 1e8); it does not model *which* higher-order mechanism produces the dependence, and
it inherits the known slow non-stationarity of gaps (which the blocked reading guards
against but does not remove). Reproduce with `python markov_order_job.py` and verify the
committed numbers with `python -m unittest -v test_markov_order_results.py`.

## Provenance and pre-registration note

The pre-registration in `notes/markov-order-prereg.md` was authored on 2026-10-18,
before the analysis was run. Its intended commit earlier that day was **skipped by the
publish lane because the v1 job carried a buggy self-test** (an idealized ±1
autocorrelation expectation that did not match the biased estimator). The bug was fixed
in the v2 job (all self-tests now pass, including the first-order/second-order synthetic
gates), and the prereg file — unchanged from when it was written — is committed here
alongside the graded result. The frozen predictions P-M1..P-M5 are stated verbatim in
that file.

Data: [`data/markov_order.json`](data/markov_order.json). Code: `markov_order.py`
(reference), `markov_order_job.py` (sandbox job). License: Apache-2.0 (code), CC BY 4.0
(text/data).
