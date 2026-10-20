# Pre-registration — Does the squarefree lag-2 signature track "admits even integers" (lift an odd prime, keep the even ones) rather than the full squarefree structure?

Produced with AI by MV2, Munim, Inc. Locked before the compute job
`admit-even-lag-signature-v1` runs (this commit precedes the job run). Code Apache-2.0;
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
(v0.18.0) showed it is **absent** for cubefree and 4th-power-free: the signature is **specific
to squarefree (k=2 / μ²)**, not a k-free family property and not a generic consequence of
multiplicative-indicator density.

A companion probe (`odd-squarefree-lag-signature-v1`, queued, pending compute) removes **the
prime 2**: odd squarefree (A056911) and squarefree coprime to 6 admit **no even integers**, so
**every gap is even** and the half-shifted stride phase (`⌊g/2⌋ = (g−1)/2` for odd `g`)
collapses to exact alignment (`⌊g/2⌋ = g/2` for even `g`). Its pre-registered prediction is
that removing even integers **destroys** the signature.

## The question isolated here (the complementary cell)

The odd-squarefree probe removes the even integers and, with them, both (a) the prime 2 from
the sieve's congruence structure and (b) the odd gaps / half-shifted stride phase — the two are
inseparable in that probe. This job fills the **complementary 2×2 cell**: lift an **odd**
prime from the squarefree sieve while **keeping all the even integers**. Squarefree coprime to
3 (and to 5, and to 15) still admit even numbers, so they still produce **odd gaps** and keep
the half-shifted stride phase intact, yet they are no longer the full squarefree set and have
had a prime removed from their congruence structure.

|                              | admits even integers? | congruence = full μ²? |
|------------------------------|-----------------------|-----------------------|
| squarefree (control)         | yes                   | yes                   |
| squarefree coprime to 3      | **yes**               | no (prime 3 lifted)   |
| squarefree coprime to 5      | **yes**               | no (prime 5 lifted)   |
| squarefree coprime to 15     | **yes**               | no (3 and 5 lifted)   |
| odd squarefree (coprime 2)   | no                    | no (prime 2 lifted)   |

> Is the squarefree lag-2 unique-quiet signature driven by the sieve **admitting even
> integers** (so odd gaps and the half-shifted base-4 stride phase survive), or does it require
> the **full** squarefree congruence structure with every prime present? Lifting an odd prime
> keeps the even integers but breaks the full-μ² structure — does the signature survive?

## Method (fixed in advance)

- Five sieves to a common **N = 20,000,000**, each analyzed with the identical observable,
  null and detector:
  - **squarefree (A005117, control, expect PRESENT)** — replicates v0.18.0 at this N.
  - **squarefree coprime to 3 (PRIMARY, expect PRESENT)** — keeps even integers, odd gaps
    survive; forbidden = (3,). Density `(6/π²)·(3/4) ≈ 0.456` → ~9.1×10⁶ elements, ample power.
  - **squarefree coprime to 5 (expect PRESENT)** — second odd-prime lift; forbidden = (5,).
  - **squarefree coprime to 15 (expect PRESENT)** — two odd primes lifted, still even-admitting;
    forbidden = (3, 5). A sharper break of full-μ² while keeping even integers.
  - **odd squarefree (A056911, anchor, expect ABSENT)** — removes even integers; forbidden =
    (2,). Replicates the companion probe's primary prediction inside this self-contained job.
- The observable, the exact permutation null `null_z`, the stride split, the stream analysis,
  the scale-free signature detector and the segmented squarefree-coprime sieve are reused
  **verbatim** and re-self-tested from the graded `odd_squarefree_lag_signature_job.py`. The
  self-test extends the exact brute-force reference match (full count, complete gap stream, max
  gap over `[1, 20000]`, tiny segment forcing many boundary crossings) to forbidden = (), (2,),
  (2,3), **(3,), (5,), (3,5)**, and additionally asserts the **mechanism contrast**: every
  odd-prime-lift sieve (forbidden (3,), (5,), (3,5)) **contains even members and ≥1 odd gap**,
  whereas the prime-2 lift (forbidden (2,)) has **all even gaps**. Plus `Q(10⁶)=607926`,
  `π(10⁶)=78498`. Self-tests run first; the job aborts (exit ≠ 0) on any failure.

## Pre-registered predictions (graded by the job against the frozen detector)

- **P-AE-sf (control).** Full squarefree shows the signature (`present == True`). *Expect HOLD.*
- **P-AE-c3 (PRIMARY).** Squarefree coprime to 3 **shows** the signature (`present == True`),
  because keeping the even integers keeps odd gaps and the half-shifted stride phase that drives
  the lag-2 XOR near-null. If **FALSIFIED** (signature ABSENT despite admitting even integers),
  the signature needs the full squarefree congruence structure, not merely even-admission — a
  sharper, equally publishable outcome.
- **P-AE-c5 (replication).** Squarefree coprime to 5 shows the signature (`present == True`).
- **P-AE-c15 (replication, double lift).** Squarefree coprime to 15 shows the signature
  (`present == True`). Divergence among c3/c5/c15 would itself be informative (prime-specific
  effect).
- **P-AE-odd (anchor).** Odd squarefree does **not** show the signature (`present == False`),
  replicating the companion probe's primary prediction. If this anchor and the c3/c5/c15
  predictions both hold, the driver is cleanly identified as **admitting even integers**, not
  the full μ² structure and not multiplicativity alone.

## Interpretation grid (locked before running)

- **c3/c5/c15 PRESENT and odd ABSENT (predicted):** the signature is driven by the sieve
  admitting even integers (odd gaps + half-shifted base-4 stride phase). "Multiplicative μ²"
  and "the prime 2 in the congruence structure" are both ruled out as the cause; parity of the
  admitted integers is the cause. Next: confirm with a non-squarefree even-admitting
  multiplicative set, and with a prime-2-lift that is *not* squarefree-derived.
- **c3/c5/c15 ABSENT (any of them):** admitting even integers is **not** sufficient; the full
  squarefree congruence structure (every prime present) is required. Next: test squarefree with
  a single large prime lifted (minimal structural change) to find the boundary.
- **Mixed (e.g. c3 present, c5 absent):** a prime-specific effect; next isolate which lifted
  prime matters and why (base-2 vs base-4 interaction with that modulus).

## Honest limits

The detector is a scale-free ratio test, so a finite N does not invalidate it, but a null means
"no signature at this power," not exact absence for all N. Lifting an odd prime changes the set
and its gap distribution together; the design keeps even-admission fixed across c3/c5/c15 so the
contrast with the odd-squarefree anchor isolates even-admission, but it does not separately vary
the base-4 residue of the lifted prime. A negative or partial result is reported as such. No
money or external engagement is attached; this is a pre-registered number-theory micro-result
by MV2 of Munim, Inc.
