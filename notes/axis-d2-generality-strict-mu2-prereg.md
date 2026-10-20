# Pre-registration — Axis D generality: does the strict-squarefree (μ²) requirement for the lag-2 gap-parity signature hold when the exponent cap is relaxed at a *larger* prime?

Pre-registered by MV2 of Munim, Inc. This note is written and committed **before** the
experiment is run. The work on this repository is produced with AI agents working in MV2's name.
Code is licensed Apache-2.0; result data is licensed CC BY 4.0.

## Where this sits in the arc

The base-2 Thue–Morse parity stream of consecutive **gaps** of a sieved integer set carries a
reproducible **"lag-2 unique-quiet"** serial-memory signature, localized to the squarefree sieve
(A005117). Across four axes the repository has shown the signature needs **exactly** the strict
squarefree congruence structure:

- **Axis A (which k-free):** only k = 2 (v0.18.0).
- **Axis B (origin vs. indicator):** needs a multiplicative *indicator*, not merely a
  multiplicative origin — the Liouville +1 level set is flat (v0.20.0).
- **Axis C (thinning):** removing even integers (odd squarefree) *and* lifting one odd prime
  (coprime 3/5/15) both destroy it (v0.21.0 / v0.23.0).
- **Axis D (enlarging):** relaxing the squarefree exponent cap at a single small prime —
  squarefree-except-2 (allow 4|n) and squarefree-except-3 (allow 9|n) — destroys it too, even
  though both keep every even integer (v0.24.0).

Axis D used only the two smallest relevant primes, 2 and 3. **This replication asks whether that
destruction is a general property of relaxing the cap at *any* single prime, or an artefact of the
small primes 2/3.** If it is general, the strict-μ² account is strengthened; if a larger-prime
relaxation *keeps* the signature, the account does not generalize and the question reopens.

## Detector (unchanged from the rest of the arc; grade base 2 only)

For a sieved set, form consecutive gaps `g_n`, the base-2 Thue–Morse parity
`x_n = P_2(g_n) = s_2(g_n) mod 2`, and the ±1 character `u_n = 1 − 2 x_n`. Measure serial memory
`S^(k) = Σ_n u_n u_{n+k}` for lags `k = 1..8` against an **exact marginal-preserving permutation
null**, reporting each as a z-score `z(k)`.

**Lag-2 unique-quiet signature is PRESENT iff**
`argmin_{k∈1..8} |z(k)| == 2` **and** `|z(2)| ≤ 0.1 × median{ |z(k)| : k ≠ 2 }`.

The question is non-tautological only at base 2/4 (see `notes/base2-even-admission-nontautology.md`),
so it is graded in base 2 exclusively.

## Test sets (all to a common N = 20,000,000)

A multiplicative indicator is fixed by the allowed exponent cap `c_p` at each prime; `n` is a
member iff `v_p(n) ≤ c_p` for every prime `p`. Default `c_p = 1` is squarefree.

- **S0 — full squarefree** (control). `c_p = 1` everywhere. Density `6/π²`.
- **S-5 — squarefree-except-5** (PRIMARY). `c_5 = 2`, `c_p = 1` otherwise: admits `25|n` but not
  `125|n`; keeps every even integer.
- **S-7 — squarefree-except-7** (replication). `c_7 = 2`, `c_p = 1` otherwise: admits `49|n` but
  not `343|n`.
- **S-anchor — odd squarefree** (anchor). `c_2 = 0` (A056911): removes evens → all gaps even.
- **S-neg — primes** (negative control, A000040).

The exponent-cap sieve, observable, exact permutation null, stride split and detector are reused
verbatim from the graded Axis D job and re-self-tested (count / full gap stream / max gap matched
against full-factorization reference over [1, 20000] for every cap set used, plus the mechanism
contrast and the exact counts Q(10⁶)=607926, π(10⁶)=78498). Self-tests run first; the job aborts
non-zero on any failure, so only self-validated numbers are emitted.

## Pre-registered predictions (PRESENT / ABSENT)

Following the Axis D result (any single-prime exponent-cap relaxation destroys the signature):

| label | set | pre-reg |
|---|---|---|
| **P-D2-sf** | full squarefree (control) | PRESENT |
| **P-D2-ex5** (PRIMARY) | squarefree-except-5 | **ABSENT** |
| **P-D2-ex7** (replication) | squarefree-except-7 | **ABSENT** |
| **P-D2-odd** (anchor) | odd squarefree | ABSENT |
| **P-D2-primes** (neg control) | primes | ABSENT |

**Decisive read.** `strict_mu2_requirement_generalizes_to_larger_primes` is True iff the control
is PRESENT while both large-prime relaxations (and the anchor and negative control) are ABSENT —
i.e. relaxing the cap at 5 or 7 destroys the signature exactly as it did at 2 and 3. If instead
S-5 and/or S-7 comes out **PRESENT**, the Axis D destruction was specific to the small primes 2/3,
the strict-μ² account does **not** generalize, and that surprise will be reported honestly. Either
direction is a real result.

## Prior art (why this is new)

The serial correlation of digit-sum parity *of consecutive integers* is classical: Gelfond 1968;
Mauduit–Rivat (Ann. of Math. 171, 2010, 1591–1646); Baake–Coons (arXiv:2209.07102, "Correlations
of the Thue–Morse sequence"). What is new here is applying the base-2 gap-parity lag-2 unique-quiet
detector to **exponent-cap-relaxed multiplicative sieves at larger primes** to test whether the
strict-squarefree requirement found in Axis D is a general property of the μ² congruence structure
rather than an artefact of the smallest primes. No prior work tests a "squarefree-except-one-larger-
prime" multiplicative indicator this way.
