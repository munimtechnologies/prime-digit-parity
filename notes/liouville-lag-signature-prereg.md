# Pre-registration — Does the base-2 gap-parity "lag-2 unique-quiet" signature need a multiplicative INDICATOR, or does it attach to multiplicative-function structure more broadly?

Produced with AI by MV2, Munim, Inc. Locked before the compute job
`liouville-lag-signature-v1` runs (this commit precedes the job run). Code Apache-2.0;
result data CC BY 4.0.

## Background (all graded and published in this repo)

For an integer sieve, form consecutive gaps `g_1, g_2, …` and the Thue–Morse stream
`x_n = P_2(g_n)` (parity of the base-2 digit sum). The exact identity
`P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` splits `x` into two base-4 stride streams. Serial memory is
graded by `S^(k) = Σ_n u_n u_{n+k}` against an exact marginal-preserving permutation null
(reported as `z`), lags `k = 1..8`.

The **lag-2 unique-quiet signature** on the `x_base2` stream (scale-free detector: present
iff `argmin_{k∈1..8}|z_x(k)| == 2` and `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was found
**present** for squarefree (A005117) and sums-of-two-squares (A001481) — whose set
indicators are multiplicative — and **absent** for primes (A000040) and the positional lucky
sieve (A000959, graded in v0.17.0). Working conclusion: the lag-2 XOR cancellation is tied to
**multiplicative** sieves, not generic to any thinned sequence.

A separate queued job (`kfree-lag-signature-v1`) asks whether the signature extends across
the **k-free family** (cubefree, 4th-power-free). This job probes an **orthogonal axis**.

## The object tested here

The **Liouville function** `λ(n) = (−1)^{Ω(n)}` is *completely* multiplicative, where `Ω(n)`
counts prime factors with multiplicity. Its **+1 level set**

> `A028260 = { n : λ(n) = +1 } = { n : Ω(n) even }`  (1, 4, 6, 9, 10, 14, 15, 21, 22, 25, …)

has natural density **1/2**. Crucially, the **indicator** `1_{Ω even}` is **NOT
multiplicative**: `1_{Ω even}(mn)` is not the product of the indicators of `m` and `n`
(e.g. `m=n=2`: `1(2)=0`, `1(4)=1 ≠ 0·0`). So this set has a clean multiplicative *origin*
(it is a level set of a completely multiplicative function) but a **non-multiplicative
indicator** — exactly the discriminator the earlier results could not separate.

## Open question graded here

> Does the lag-2 unique-quiet signature require a **multiplicative indicator** (as
> squarefree / SoTS / k-free sets have), or does it attach to **multiplicative-function
> structure** more broadly — appearing even for the Liouville +1 level set, whose indicator
> is not multiplicative?

This is informative **independently of** the k-free job's outcome: together the two jobs give
a 2-D picture {k-free family: yes/no} × {Liouville level set: yes/no}, and every cell is
interpretable (indicator-multiplicative vs. multiplicative-origin vs. squarefree-specific).

## Method (fixed in advance)

- Three sieves to a common **N = 20,000,000**: **Liouville +1 (A028260, primary)**,
  squarefree (A005117, multiplicative-indicator control, expect present), primes (A000040,
  negative control, expect absent). The Liouville +1 set has density 1/2 → ~10⁷ elements and
  ~10⁷ gaps — ample power; gaps are tiny (array("H") safe).
- The observable, the exact permutation null `null_z`, the stride split, the stream analysis,
  the scale-free signature detector, the squarefree sieve and the prime sieve are reused
  **verbatim** and re-self-tested from the graded `kfree_lag_signature_job.py` /
  `lucky_lag_signature_job.py`. The only new code is a **segmented exact Ω-parity (Liouville)
  sieve**, re-self-tested against **brute-force trial-division Ω** across many forced segment
  boundaries (full count *and* the complete gap stream over [1, 20000] must match), plus the
  exact small set (A028260 ≤ 10 = {1,4,6,9,10}, count 5, consistent with L(10)=0). Self-tests
  run first; the job aborts (exit≠0) on any failure, so only self-validated numbers are
  emitted.

## Pre-registered predictions (graded by the job against the frozen detector)

- **P-L-sf (control).** Squarefree shows the signature (`present == True`), replicating
  v0.17.0 at this N. *Expected to HOLD.*
- **P-L-primes (negative control).** Primes do **not** show the signature
  (`present == False`). *Expected to HOLD.*
- **P-L-liouville (PRIMARY).** The Liouville +1 level set does **not** show the signature
  (`present == False`), because its indicator is not multiplicative and the signature is
  hypothesized to require a multiplicative indicator. If **FALSIFIED** (the Liouville set
  **does** show it), the signature attaches to multiplicative-function structure more broadly
  than a multiplicative indicator — a surprising, sharper, equally publishable outcome.

## Honest limits

The signature test is a scale-free ratio test, so a finite N does not invalidate it, but a
null means "no signature at this power" rather than a proof of exact absence for all N. A
single level set (λ = +1) is one object; a positive or negative result constrains, but does
not by itself settle, the general multiplicative-vs-indicator question. A negative or partial
result is reported as such. No money or external engagement is attached; this is a
pre-registered number-theory micro-result by MV2 of Munim, Inc.
