# Pre-registration — Does the base-2 gap-parity "lag-2 unique-quiet" signature need the squarefree sieve to admit even numbers?

Produced with AI by MV2, Munim, Inc. Locked before the compute job
`odd-squarefree-lag-signature-v1` runs (this commit precedes the job run). Code Apache-2.0;
result data CC BY 4.0.

## Background (all graded and published in this repo)

For an integer sieve, form consecutive gaps `g_1, g_2, …` and the Thue–Morse stream
`x_n = P_2(g_n)` (parity of the base-2 digit sum). The exact identity
`P_2(g) = P_4(g) XOR P_4(⌊g/2⌋)` splits `x` into two base-4 stride streams. Serial memory is
graded by `S^(k) = Σ_n u_n u_{n+k}` against an exact marginal-preserving permutation null
(reported as `z`), lags `k = 1..8`.

The **lag-2 unique-quiet signature** on the `x_base2` stream (scale-free detector: present
iff `argmin_{k∈1..8}|z_x(k)| == 2` and `|z_x(2)| ≤ 0.1 × median{|z_x(k)| : k≠2}`) was found
**present** for squarefree (A005117) and sums-of-two-squares (A001481), and **absent** for
primes (A000040) and the positional lucky sieve (A000959, v0.17.0). The k-free follow-up
(v0.18.0) then showed it is **absent** for cubefree (k=3) and 4th-power-free (k=4): the
signature is **specific to squarefree (k=2 / μ²)**, not a k-free family property and not a
generic consequence of multiplicative-indicator density.

## The mechanism tested here

The full squarefree set admits **both even and odd** integers, so its gap sequence contains
**odd gaps** (e.g. 2→3 and 5→6 and 6→7 are gaps of 1) as well as even gaps. Two base-4
stride streams `a_n = P_4(g_n)` and `b_n = P_4(⌊g_n/2⌋)` combine by XOR to give `x`; the
lag-2 near-null was identified as an XOR cancellation between these two streams that depends
on their relative phase, which in turn depends on whether `g` is odd (⌊g/2⌋ = (g−1)/2) or
even (⌊g/2⌋ = g/2 exactly).

Restricting to **odd squarefree numbers** (A056911) removes every even element, so **every
gap becomes even**. For even `g`, `P_2(g) = P_2(g/2)`, and `b_n = P_4(g_n/2)` is exactly
aligned with `a_n = P_4(g_n)` rather than half-shifted. If the lag-2 cancellation depends on
the sieve admitting even numbers (hence producing odd gaps and the half-shifted stride
phase), this restriction should **destroy** the signature. The indicator `μ²(n)·1_{2∤n}`
remains a multiplicative indicator, so the sieve stays inside the class that *did* show the
signature — isolating the role of admitting even numbers from the role of multiplicativity.

## Open question graded here

> Is the squarefree lag-2 unique-quiet signature a property of the multiplicative squarefree
> indicator per se, or does it require the sieve to **admit even integers** (and therefore
> produce odd gaps)? Removing the even elements removes the odd gaps without leaving the
> multiplicative-indicator class — does the signature survive?

## Method (fixed in advance)

- Four sieves to a common **N = 20,000,000**: **squarefree (A005117, control, expect
  present)**, **odd squarefree (A056911, PRIMARY, expect absent)**, **squarefree coprime to
  6** (a second all-even-gap multiplicative variant that also drops multiples of 3, expect
  absent), and **primes (A000040, negative control, expect absent)**. Odd squarefree have
  density `3/π²` → ~6×10⁶ elements and gaps; ample power; gaps stay tiny (array("H") safe).
- The observable, the exact permutation null `null_z`, the stride split, the stream analysis,
  the scale-free signature detector and the prime sieve are reused **verbatim** and
  re-self-tested from the graded `kfree_lag_signature_job.py`. The only new code is a
  **segmented squarefree sieve with a coprimality filter**, re-self-tested against a
  brute-force trial-division reference: the full count, the **complete gap stream** and the
  max gap over `[1, 20000]` must match exactly, with a deliberately tiny segment forcing many
  segment-boundary crossings, for each of forbidden = (), (2,) and (2,3); plus the exact
  `Q(10⁶)=607926` squarefree count and `π(10⁶)=78498`. Self-tests run first; the job aborts
  (exit ≠ 0) on any failure, so only self-validated numbers are emitted.

## Pre-registered predictions (graded by the job against the frozen detector)

- **P-OS-sf (control).** Full squarefree shows the signature (`present == True`), replicating
  v0.18.0 at this N. *Expected to HOLD.*
- **P-OS-odd (PRIMARY).** Odd squarefree does **not** show the signature (`present == False`),
  because forcing all gaps even collapses the half-shifted stride phase hypothesized to drive
  the lag-2 XOR near-null. If **FALSIFIED** (the signature **survives** under the even-number
  restriction), it attaches to the multiplicative squarefree indicator itself rather than to
  the presence of odd gaps — a sharper, equally publishable outcome.
- **P-OS-coprime6 (variant).** Squarefree coprime to 6 does **not** show the signature
  (`present == False`), corroborating P-OS-odd with a second all-even-gap variant. Divergence
  between the two variants would itself be informative.
- **P-OS-primes (negative control).** Primes do **not** show the signature (`present ==
  False`). *Expected to HOLD.*

## Honest limits

The signature test is a scale-free ratio test, so a finite N does not invalidate it, but a
null means "no signature at this power" rather than a proof of exact absence for all N.
Removing even integers changes both the admitted set and the gap parity at once; the two are
inseparable (removing even integers necessarily makes every gap even), so this experiment
tests exactly that combined hypothesis and does not separate "the prime 2 in the sieve" from
"odd gaps in the stream". A negative or partial result is reported as such. No money or
external engagement is attached; this is a pre-registered number-theory micro-result by MV2
of Munim, Inc.
