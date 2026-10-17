# Lag-k serial correlation of prime-gap digit-sum parity

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.*

## The question

v0.7.0 (`GAP-CORR.md`) showed that consecutive prime-gap digit-sum parities
χ_b(g_n), χ_b(g_{n+1}) carry strong, base-structured **lag-1** serial memory
beyond their fixed parity histogram. It left one question open — the one its own
"honest limits" named: **does that memory persist, decay, or oscillate at longer
range?** For lag k ≥ 1, base b ∈ {2,4,6,8,10,16}, cutoff 10^6/10^7/10^8 define

    S_b^(k) = Σ_{n=1}^{N−k} χ_b(g_n)·χ_b(g_{n+k}),   χ_b(m) = (−1)^{s_b(m)},  N = n_gaps,

the lag-k signed adjacent-pair sum (same-parity minus different-parity count at
separation k). The signal is the connected deviation z_b^(k) of S_b^(k) from a
marginal-preserving (uniform-permutation) null that fixes the exact multiset of
χ values. We measure k = 1,…,5. k = 1 must reproduce v0.7.0 exactly (a regression
gate), so every number here is on the same footing as the shipped lag-1 result.

## Why this is new

- The nearest art — Rudin–Shapiro / Thue–Morse correlation **along primes**
  (arXiv:2306.10629; sums Σ_{p≤N} r_p r_{p+k}) and the classical Thue–Morse
  autocorrelation η(m) = lim (1/N)Σ t_k t_{k+m} (Baake–Grimm, math/0610791) —
  correlates the TM/RS *word itself* (or its evaluation at primes), a different
  object from χ_b(g_n), the digit-sum parity of the prime **gaps** correlated
  across gap index.
- Lemke Oliver–Soundararajan (arXiv:1604.02443) and the sums-of-two-squares
  analogue (arXiv:2105.05048) bias the **residue classes** of consecutive primes,
  not the digit-sum parity of the gaps at any lag.
- This repo's v0.7.0 established lag-1 only. The lag-k profile — whether the
  base-structured memory is range-1, decaying, or oscillating — is the new result.

Eight fresh prior-art searches (runs 47, 48, 50, 52) plus one more at this
shipping run found no source computing the lag-k serial correlation of χ_b(g_n). → novel.

## Exact permutation-null mean (hand-derived — gates the empirical null)

Fix the multiset of the N values χ_b(g_n); Σχ = −E_b^gap, Σχ² = N. For a uniform
random permutation, any lag-k pair (positions n, n+k, k ≥ 1) is a uniform random
*ordered pair of distinct* positions, so each contributes (E_b²−N)/(N(N−1)) and
there are N−k of them:

    **E_perm[S_b^(k)] = (N−k)·(E_b² − N)/(N(N−1)).**

At k = 1 this reduces to (E_b²−N)/N, the v0.7.0 lag-1 form. `test_lagk_results.py`
brute-forces this general closed form over *every* permutation of small ±1
multisets for k = 1,2,3. In the sealed run the empirical K-shuffle mean matched it
to within 5·std/√K for every k (a hard gate, passed).

## Anchors (gate — all must pass or the number is a bug, not a result)

The sealed CPU-sandbox run printed, **before any new number**:

> `method gate PASSED (pi, telescoping, single odd gap, half-gap identity,`
> `two-routine guard, v0.6.2 excess reproduced).`
> `null-mean gate PASSED (empirical shuffle mean == (N-k)(E^2-N)/(N(N-1)) for every lag).`
> `v0.7.0 lag-1 regression PASSED (S_b^(1) reproduces GAP-CORR exactly).`

π(x) = 78498 / 664579 / 5761455; telescoping Σg_n = p_M − 2; one odd gap; the
half-gap Thue–Morse identity χ_2(g) = χ_2(g/2); a two-routine digit guard on the
first 200 000 gaps every base; exact reproduction of the v0.6.2 excess E_b^gap;
and — new here — the k=1 column reproducing v0.7.0's S_b^(1) per base and cutoff
(e.g. base 2 @10^8 = −166 721, base 16 @10^8 = 508 353).

## Result

Gate: **PASSED** (all three gates, sealed CPU-sandbox run log).

Connected signal z_b^(k) = (S_b^(k) − shuffle_mean)/shuffle_std at 10^8
(n_gaps = 5 761 454, K = 200 shuffles). Full per-cutoff S, E_perm, shuffle_mean,
shuffle_std and K for all three cutoffs are in `data/gap_parity_lagk.csv`.

| base |   k=1  |  k=2  |  k=3  |  k=4  |  k=5  |
|-----:|-------:|------:|------:|------:|------:|
| 2  | −120.58 | **+12.41** |  −0.63 |  −0.99 |  −0.40 |
| 4  |  **+18.10** |  −3.58 |  −5.11 |  −2.05 |  −0.97 |
| 6  |  −3.48 |  +4.07 |  +0.17 |  +0.03 |  −1.96 |
| 8  | −33.25 | −10.47 |  +1.88 |  −0.45 |  −1.32 |
| 10 |  −4.39 | **−11.77** |  −2.61 |  −1.17 |  −2.31 |
| 16 | −58.15 | −12.86 |  −4.20 |  −3.12 |  −1.10 |

The 10^7 column tells the same story (base 2: −47.25, +6.68, −0.23, …; base 16:
−20.50, −3.76, −1.78, −1.46, −3.32; base 8: −11.96, −7.26, …), so the lag profiles
below are stable across a decade of x, not a 10^8 fluctuation.

## Which hypothesis the data support (decided from the z-table; prereg run 49)

Registered **before** the numbers (`notes/gap-parity-lagk-prereg.md`):

- **D0 (null):** no memory beyond the histogram at any lag; |z| = O(1) for all b, k.
- **D1 (decaying memory):** lag-1 structure decays with k.
- **D2 (periodic / base-specific memory):** structure at a specific lag/base, not a
  monotone decay.

**Verdict: not D0 — a mix of D1 and D2. The serial memory extends beyond one step,
and its lag profile is strongly base-dependent; there is no single universal
decay.** Four distinct, reproducible shapes appear:

- **Base 2 — a clean sign-alternation of range exactly 2.** Adjacent gaps carry
  *opposite* parity (k=1: −120.6), gaps *two* apart carry the *same* parity
  (k=2: **+12.4**, a sign flip), then the signal collapses to noise (|z| < 1 for
  k = 3,4,5). The lag-2 reversal is itself a 12σ effect and reproduces at 10^7
  (+6.7). This is the sharpest new structure: binary-digit-sum parity of prime
  gaps behaves like a damped period-2 oscillation.
- **Base 8 — same-sign anti-correlation, range 2.** Both k=1 (−33.3) and k=2
  (−10.5) are anti-correlated (no sign flip), then noise. Memory of range ~2.
- **Base 16 — anti-correlation decaying monotonically, the longest range.** k=1
  (−58.2) → k=2 (−12.9) → k=3 (−4.2) → k=4 (−3.1), only dying by k=5. Clean D1
  decay reaching lag ~4.
- **Base 4 — positive then negative** (the opposite handedness to base 2): k=1
  (+18.1, clustering) flips to anti-correlation at k=2 (−3.6) and k=3 (−5.1).

A clean cross-base regularity: **the power-of-two bases (2, 8, 16) carry the
strongest multi-lag memory, and the memory *range* grows with the base** — base 2
effectively range 2 (with a sign flip), base 8 range 2 (same sign), base 16
range ~4. Non-power-of-two bases 6 and 10 stay weak at lag 1, and base 10's
signal actually *peaks at lag 2* (−11.8 vs −4.4) rather than lag 1 — though this
lag-2 dominance only appears at 10^8 (at 10^7 base 10 is −5.8, −3.8), so it is
reported as emerging, not settled.

So the order of prime-gap digit-sum parities is not a parity-independent reshuffle
at any lag ≤ 5: it carries real, finite-range, base-structured serial memory whose
*shape* (oscillating vs. decaying, how far it reaches) is set by the base.

## Honest limits

- Lags k ≤ 5 only; deeper memory not probed. The sign/shape interpretations are
  empirical descriptions, not an explanation of *why* a base oscillates (2, 4) vs.
  decays monotonically (16). Base 10's lag-2 peak is scale-dependent (clear only at
  10^8) — treat as suggestive.
- Measured to 10^8 (5.76 M gaps); K-shuffle null (2000/800/200) on a single
  machine, single sieve. The reproduce-v0.6.2, reproduce-v0.7.0 and null-mean
  identities are the internal cross-checks standing in for an external anchor; the
  exact permutation variance (v0.8.0 `null_moments_lagk.py`) gives a deterministic
  analytic z for any lag with zero Monte-Carlo error in the denominator.
- An empirical pattern in a finite range, not a theorem; the growth of |z| with the
  cutoff is suggestive but unproven asymptotically.

## Reproduce

```
python -m unittest -v test_lagk_results   # anchors + brute-force lag-k null mean, stdlib only
```

Full measurement (sealed CPU sandbox): `gap_lagk_job.py --out-dir out --limit 100000000`.
Data: `data/gap_parity_lagk.csv` (CC BY 4.0). Code Apache-2.0.

— MV2, Munim, Inc.
