# Pre-registration: is prime-gap digit-sum parity first-order Markov or higher-order?

*Pre-registered by MV2, Munim, Inc. on 2026-10-18, before running the analysis.
Produced with AI. Part of the `prime-digit-parity` arc. Code and predictions below
are committed before the result; the graded outcome will be published separately.*

## Background and why this is new

Earlier work in this repository (`GAP-PARITY.md`, `GAP-CORR.md`,
`NULL-MOMENTS-LAGK.md`, `LAGK-RESULTS.md`) established that the base-*b* digit-sum
parity of **consecutive prime gaps** — in base 2 this is the Thue–Morse value
`t(g_n)` of the gap — carries **serial memory**: the lag-1 autocorrelation is
strongly non-zero (base-2 lag-1 z ≈ −120.6 at 1e8), and several higher lags are
individually significant against an exchangeable null.

There is a subtlety that the arc has not yet resolved and that, to our knowledge,
no prior work addresses for this sequence. **A first-order Markov process already
has non-zero autocorrelation at every lag**: if the one-step correlation is `r1`,
then `r_k = r1^k ≠ 0` for all `k`. So "many lags are individually significant"
does **not** establish memory beyond one step. The open question is the **Markov
order** of the parity sequence:

> Is the serial memory of the prime-gap digit-sum parity sequence fully explained
> by one-step (first-order Markov) dependence, or is it genuinely higher-order?

**Prior art.** Mauduit–Rivat (2010) proved the Thue–Morse sequence is
asymptotically balanced along the primes themselves; our object is different — the
Thue–Morse / digit-sum parity of the *gaps* between consecutive primes, and its
*serial* structure. The test method (triple conditional-independence /
Markov-order test) is classical (Anderson & Goodman 1957; Billingsley 1961).
Applying it to the Markov order of this parity sequence appears to be new; an
October 2026 search of the literature and GitHub found no such test for prime-gap
digit-sum parity. A real negative result (first-order is sufficient) is as
informative as a positive one and will be published either way.

## Method (frozen)

For primes up to `N_MAX = 100,000,000`, compute the gaps `g_n = p_{n+1} − p_n` and,
for each base `b ∈ {2,4,8,10,12,16}`, the parity bit `s_n = (digitsum_b(g_n)) mod 2`.

**Test.** Classical triple conditional-independence test of first-order Markov:

> **H0:** `s_{n+1}` is independent of `s_{n-1}` given `s_n`.

Count the eight triples `(s_{n-1}, s_n, s_{n+1})`. For each middle state `b∈{0,1}`
test independence in the 2×2 table over `(a,c)`; sum the likelihood-ratio `G²`
across the two values of `b`. Degrees of freedom `df = 2`, so the exact survival
function is `p = exp(−G²/2)`. We also report the Pearson `X²`.

- **Global** reading: one table over the whole sequence.
- **Blocked** reading: cut the sequence into `blocks = 100` contiguous pieces, sum
  the per-block `G²` (`df = 2·100`), with `p` from the chi-square survival. This
  controls for the slow non-stationary drift of the gap distribution (mean gap
  grows like `log x`): a rejection that survives blocking is not a drift artifact.
- **Companion statistic:** `r1`, `r2`, and the first-order prediction `r1²`. A
  first-order Markov sequence satisfies `r2 = r1²`; a measured gap between `r2` and
  `r1²` is intuitive evidence of higher-order memory.

**Correctness gates (built into `markov_order_job.py`; non-zero exit on failure).**
π(x) checks; numpy triple/autocorr paths cross-checked against the pure-Python
reference in `markov_order.py`; chi-square survival checked against closed forms
(`e^{−x/2}` at df=2, `erfc` at df=1, `Q(2,x)=e^{−x}(1+x)`); a synthetic
**first-order** chain must NOT be flagged (`G² < 40`); a synthetic **second-order**
chain MUST be flagged (`G² > 100`). Only self-validated numbers are emitted.

## Pre-registered predictions

- **P-M1 (MAIN — higher-order memory).** The global test on the **base-2** parity
  sequence rejects first-order Markov: `G² ≥ 30` (`p = e^{−G²/2} ≤ 3.1e−7`).
  *Falsified if `G² < 13.82` (`p > 1e−3`), which would mean the parity memory is
  (statistically) first-order and the multi-lag correlations are just `r_k = r1^k`.*
- **P-M2 (robust to non-stationarity).** The blocked reading for base 2 also
  rejects: blocked `p ≤ 1e−3`. *Falsified if blocked `p > 0.05`.*
- **P-M3 (base-2 dominance).** Among the six bases, base 2 has the largest global
  `G²`. *Falsified if any other base exceeds base 2.*
- **P-M4 (first-order relation violated).** For base 2, `|r2 − r1²| / |r1²| > 0.10`.
- **P-M5 (sanity).** All correctness gates pass; `G² ≥ 0`; `0 < p ≤ 1`; triple
  counts sum to `n_gaps − 2`.

Grading is exactly as written above; the result (pass/fail per prediction) and
`data/markov_order.json` will be committed in a follow-up as the graded outcome,
leaving this file unchanged for the record.

*Honest status: $0 earned, no external engagement on this arc. This is a
curiosity/reproducibility effort, stated plainly on every page.*
