# Pre-registration — Is the base-2 gap-parity "lag-2 unique-quiet" signature a k-free family property or specific to squarefree (k=2)?

Produced with AI by MV2, Munim, Inc. Locked before the compute job `kfree-lag-signature-v1`
runs (this commit precedes the job run). Code Apache-2.0; result data CC BY 4.0.

## Background (all graded and published in this repo)

For an integer sieve, form consecutive gaps `g_1, g_2, …` and the Thue–Morse stream
`x_n = P_2(g_n)` (parity of the base-2 digit sum). The exact identity
`P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` splits `x` into two base-4 stride streams. Serial memory is
graded by `S^(k) = Σ_n u_n u_{n+k}` against an exact marginal-preserving permutation null
(reported as `z`), lags `k = 1..8`.

The **lag-2 unique-quiet signature** on the `x_base2` stream (scale-free detector: present
iff `argmin_{k∈1..8}|z_x(k)| == 2` and `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was found
**present** for squarefree (A005117) and sums-of-two-squares (A001481), and **absent** for
primes (A000040) and, in the just-graded lucky-lag-signature result (v0.17.0), for the
positional lucky sieve (A000959). Conclusion so far: the lag-2 XOR cancellation is tied to
multiplicative sieves, not generic to any thinned sequence.

Squarefree is the **k = 2** member of the **k-free family**: `n` is *k-free* iff no prime
power `p^k` divides `n`. Its indicator `∏_p [v_p(n) < k]` is multiplicative for every `k ≥ 2`.

## Open question graded here

> Is the lag-2 unique-quiet signature a property of the whole **k-free family**
> (so **cubefree** k=3 and **4th-power-free** k=4 also show it), or is it **specific to
> squarefree (k=2)**?

## Method (fixed in advance)

- Four sieves to a common **N = 20,000,000**: squarefree (k=2, control), **cubefree (k=3,
  primary)**, 4th-power-free (k=4), and primes (negative control). All are dense
  (k-free density = 1/ζ(k): 0.608, 0.832, 0.924), giving >1.2·10⁷ gaps each — ample power.
- The observable, the exact permutation null `null_z`, the stride split, the stream analysis
  and the scale-free signature detector are reused **verbatim** and re-self-tested from the
  graded `lucky_lag_signature_job.py`. The only new code is a generalized k-free segmented
  sieve, re-self-tested against exact small counts (squarefree ≤100 = 61, cubefree ≤100 =
  85, 4th-power-free ≤100 = 93, squarefree ≤10⁶ = 607926) and the exact prime count
  (π(10⁶) = 78498). Self-tests run first; the job aborts (exit≠0) on any failure.

## Pre-registered predictions (graded by the job against the frozen detector)

- **P-K-2 (control).** Squarefree shows the signature (`present == True`), replicating
  v0.17.0 at this N. *Expected to HOLD.*
- **P-K-3 (PRIMARY).** Cubefree shows the signature (`present == True`): the signature
  generalizes to k=3, so it is a **k-free family property**, not mu²-specific. If
  **FALSIFIED** (cubefree lacks it), the lag-2 cancellation is **specific to squarefree
  (k=2)** — a sharper, equally publishable outcome.
- **P-K-4.** 4th-power-free shows the signature (`present == True`).
- **P-K-primes (negative control).** Primes do **not** show the signature (`present ==
  False`). *Expected to HOLD.*

## Honest limits

A positive result for P-K-3/P-K-4 shows the signature across k = 2,3,4, strong evidence it
is a k-free family property, but not a proof for all k or all multiplicative-indicator
sieves. The signature test is a ratio test (scale-free), so it is not invalidated by the
finite N, but a null means "no signature at this power". A negative or partial result is
reported as such. No money or external engagement is attached; this is a pre-registered
number-theory micro-result.
