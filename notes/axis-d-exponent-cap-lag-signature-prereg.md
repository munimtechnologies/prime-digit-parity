# Pre-registration — Axis D: does the lag-2 gap-digit-sum-parity signature require squarefreeness (μ², every exponent ≤ 1), or only a multiplicative indicator that admits even integers?

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the
Axis D experiment is run. The work on this repository is produced with AI agents working
in MV2's name. Code is licensed Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

Across this repository the base-2 Thue–Morse parity stream of consecutive **gaps** of a
sieved integer set carries a reproducible **"lag-2 unique-quiet" serial-memory signature**.
It is PRESENT for the full squarefree sieve (A005117) and for sums-of-two-squares, and
ABSENT for primes and for the positional "lucky" sieve. The k-free experiment (v0.18.0)
showed it is specific to **k = 2** (squarefree): cubefree and 4th-power-free lack it.

Two experiments now in the compute queue probe **what inside squarefreeness** drives it:

- **Axis B (`liouville-lag-signature`, pending):** the Liouville λ = +1 level set
  (A028260) — a completely-multiplicative-function level set whose *indicator* is itself
  non-multiplicative. Pre-registered ABSENT.
- **Axis C (`odd-squarefree`, pending; `admit-even`, designed):** a 2×2 that isolates
  "admits even integers / produces odd gaps" from "μ²". Odd-squarefree and
  squarefree-coprime-6 (both force all gaps even) pre-registered ABSENT; squarefree
  coprime to 3, 5, 15 (keep even integers) pre-registered PRESENT.

**Axis D asks the complementary question Axis C cannot:** if even-admission is the driver,
is the exact squarefree exponent cap (every prime exponent ≤ 1) still *necessary*, or does
**any multiplicative 0/1 indicator that admits even integers** carry the signature even
when it allows some square divisors?

## Detector (unchanged from the rest of the arc; grade base 2 only)

For a sieved set, form consecutive gaps `g_n`, the base-2 Thue–Morse parity
`x_n = P_2(g_n) = s_2(g_n) mod 2`, and the ±1 character `u_n = 1 - 2 x_n`. Measure serial
memory `S^(k) = Σ_n u_n u_{n+k}` for lags `k = 1..8` against an **exact marginal-preserving
permutation null** (fix the counts of each parity, permute order), reporting each as a
z-score `z(k)`.

**Lag-2 unique-quiet signature is PRESENT iff**
`argmin_{k∈1..8} |z(k)| == 2` **and** `|z(2)| ≤ 0.1 × median{ |z(k)| : k ≠ 2 }`.

Base 2 only: by the base-2 caveat recorded in `base2-even-admission-nontautology.md`, this
question is only non-tautological at base 2/4, so Axis D is graded in base 2 exclusively.

## Test sets (all multiplicative 0/1 indicators; natural density in brackets)

A multiplicative indicator is fixed by the allowed exponent set `A_p` at each prime; the
density that `v_p(n) ∈ A_p` is `(1 − 1/p) Σ_{a∈A_p} p^{-a}`, and the set density is the
product over primes.

- **S0 — full squarefree** (control). `A_p = {0,1}` at every p. Density `6/π² ≈ 0.6079`.
- **S-D1 — squarefree-except-2** (PRIMARY). `A_2 = {0,1,2}`, `A_p = {0,1}` for odd p.
  Admits the squares `4, 12, 20, 28, 36, …` (all with `v_2 = 2`), so it is **not**
  squarefree, yet it still admits even integers. Density `7/π² ≈ 0.7094`.
- **S-D2 — squarefree-except-3** (replication). `A_3 = {0,1,2}`, `A_p = {0,1}` otherwise.
  Admits `9, 18, 36, 45, …`; relaxes the cap at an **odd** prime, controlling that any
  effect of relaxing it is not specific to p = 2. Density `13/(2π²) ≈ 0.6586`.
- **S-anchor — odd squarefree** (anchor). Forces all gaps even → parity collapse.
- **S-neg — primes** (negative control).

Each set is generated to the same cutoff on the natural numbers with a linear sieve; the
same detector and the same exact permutation null are applied to each.

## Pre-registered predictions (PRESENT / ABSENT)

On the hypothesis that **even-admission (odd gaps), not the exact μ² exponent cap, is the
driver** — the branch Axis C is pre-registered to support:

| label | set | pre-reg |
|---|---|---|
| **P-D-sf** | full squarefree (control) | PRESENT |
| **P-D-ex2** (PRIMARY) | squarefree-except-2 | **PRESENT** |
| **P-D-ex3** (replication) | squarefree-except-3 | **PRESENT** |
| **P-D-odd** (anchor) | odd squarefree | ABSENT |
| **P-D-primes** (neg control) | primes | ABSENT |

**Decisive falsification.** If S-D1 (and/or S-D2) comes out **ABSENT** while S0 is PRESENT,
the "even-admission alone" hypothesis is falsified and a strict exponent-cap (μ²)
requirement is reinstated — i.e. allowing even a single squared prime destroys the
signature. Either direction is a real result and will be reported honestly, including a
clear statement if the primary prediction fails.

## Gating (do not queue out of order)

Axis D is **downstream of Axis C** and must be queued only **after** Axis C has landed and
been graded:

- If Axis C resolves **as pre-registered** (admit-even PRESENT at coprime-3/5/15,
  odd-squarefree ABSENT), queue Axis D as specified above.
- If Axis C **falsifies** even-admission (any coprime-k ABSENT), Axis D is **withdrawn** and
  replaced by the single-large-prime-lift probe recorded in the playbook, because the
  driver would then be the full squarefree congruence structure rather than even-admission.

One job per name; the detector, the null and the test-set generators are specified here in
full so the eventual job only executes this plan.

## Prior art (why this is new)

The serial correlation of digit-sum parity *of consecutive integers* is classical:
Gelfond (1968); Mauduit–Rivat (Ann. of Math. 171); Bésineau (Acta Arith. 20, 1972,
statistical independence of digit-sum sets); Drmota–Kauers–Spiegelhofer on the Cusick
conjecture; "Decomposing the sum-of-digits correlation measure" (arXiv 2411.07779). What
is new here is applying this parity detector to the **gap sequences of exponent-restricted
multiplicative sieves** and using the lag-2 unique-quiet detector to separate which single
structural feature (even-admission vs. the exact μ² exponent cap) the signature depends on.
No prior work tests a "squarefree-except-one-prime" multiplicative indicator this way.
