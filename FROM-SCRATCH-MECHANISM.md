# From-scratch reproduction of the two mechanism results (Liouville flat, odd-squarefree loud)

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. Code is Apache-2.0;
this text is CC BY 4.0.*

## What this adds

Two of the headline mechanism findings in this repository were each established by a
**single sealed compute job at N = 2×10⁷** and then audited by a test that re-derives
the reported numbers from that job's own committed JSON:

- **Axis B** — [`LIOUVILLE-LAG-SIGNATURE-RESULTS.md`](LIOUVILLE-LAG-SIGNATURE-RESULTS.md)
  (v0.20.0): the Liouville **+1 level set** A028260 (`{n : Ω(n) even}`, i.e. λ(n) = +1)
  does **not** carry the lag-2 unique-quiet signature. Its whole base-2 gap-parity
  z-profile is **flat**. Conclusion: the signature needs a *multiplicative indicator*,
  not merely a multiplicative origin — λ is multiplicative, but the indicator of its +1
  level set is not.
- **Axis C** — [`ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md`](ODD-SQUAREFREE-LAG-SIGNATURE-RESULTS.md)
  (v0.21.0): the **odd** squarefree numbers A056911 (the prime 2 lifted, so every gap is
  even) do **not** carry it either — they are **loud and uniform** at every lag.
  Conclusion: the driver is the sieve **admitting even integers**, not multiplicativity.

Re-deriving from committed JSON is good for auditing one specific run, but it trusts that
run's data. `test_from_scratch_mechanism.py` trusts **no committed data**. With the Python
standard library only it:

1. sieves four sequences in `[1, N]` from scratch — squarefree, **odd** squarefree,
   the **Liouville +1 level set**, and the primes;
2. forms each one's consecutive gaps and maps every gap `g` to its base-2 digit-sum
   parity `t(g)` (the Thue–Morse bit);
3. recomputes, for lags 1…12, the agreement z-score of each parity stream under an
   **exact marginal-preserving permutation null**; and
4. checks the four distinct fingerprints that underlie the published conclusions.

A reader can run `python -m unittest test_from_scratch_mechanism` and watch all four
fingerprints appear in a few seconds, with no downloads and no reliance on any earlier
committed result in this repo.

## The Ω sieve and the null are both proven, not assumed

- **Ω (big-omega) sieve.** The Liouville level set needs the parity of Ω(n) (prime
  factors with multiplicity). This file builds Ω by an **additive prime-power sieve** —
  a code path completely independent of the sealed job — and `TestOmegaSieve` checks it
  against independent **trial-division** factorisation for every `n ≤ 5000`, plus the
  textbook λ-signs (λ(1)=λ(4)=λ(6)=λ(60)=+1; λ(2)=λ(8)=λ(12)=λ(30)=−1).
- **Permutation null.** For a fixed bit multiset (n₁ ones among L bits) the lag-k
  agreement count has an exact mean and variance under a uniform random permutation;
  the closed form in falling factorials is checked against **exhaustive brute force**
  over all permutations for five small cases in `TestNullMomentFormula`.

## The four fingerprints (recomputed at N = 2,000,000)

| sequence | base-2 gap-parity fingerprint | verdict |
|---|---|---|
| squarefree (control) | lag-2 is the **uniquely quiet** lag; lag-1 louder than lag-2 | **signature fires** |
| odd squarefree (Axis C) | lag 2 is **loud** (louder than squarefree's lag 2); no unique quiet lag | absent |
| Liouville +1 set (Axis B) | **flat** — loudest lag far below the squarefree control's loudest lag | absent |
| primes (control) | lag 2 is **not** the quietest lag | absent |

The decisive contrast is one place (`test_contrast_is_sharp`): the **same** squarefree
sieve fires (quiet lag-2) when even integers are present, and does **not** fire — with a
**loud** lag-2 — once the prime 2 is lifted (odd squarefree), isolating "admits even
integers" as what the lag-2 cancellation needs. That is exactly the Axis C conclusion,
now reproduced independently.

## Why the test asserts fingerprints, not absolute |z| thresholds

The autocorrelation z-scores scale like √N, so their absolute magnitudes at this test's
N = 2×10⁶ are smaller than at the sealed job's N = 2×10⁷ (for example the squarefree
lag-1 peak is ≈ 225 at 2×10⁷ and ≈ 70 at 2×10⁶). A test calibrated with hard magnitude
constants at one N is therefore fragile at another. The published conclusions are
**qualitative fingerprints** — fires / loud-and-uniform / flat / not-quietest — and the
**relative ordering** between sieves, both of which are N-robust. So every mechanism
assertion in this file is one of:

- a **detector verdict**: does the published lag-2 detector fire, and which lag is the
  quietest (`lag2_signature` → `present`, `quiet`); or
- a **relative comparison** between two sieves at the same N: squarefree lag-1 louder
  than its own quiet lag-2; odd-squarefree lag-2 louder than squarefree lag-2; the
  Liouville profile's loudest lag below the squarefree control's loudest lag.

The exactness checks are unaffected and remain exact at any N: the Ω sieve versus
trial-division factorisation, and the null's closed-form moments versus exhaustive
brute force.

## Honest limits

- This **confirms**, it does not extend, the two published mechanism results. The
  absolute z magnitudes are smaller than at N = 2×10⁷; the *qualitative fingerprints*
  (fires / loud-uniform / flat / not-quietest) and their relative ordering are what is
  reproduced.
- "Flat" for the Liouville +1 level set is stated as a relative bound (its loudest lag
  is far below the squarefree control's loudest lag) together with the detector's
  non-firing, rather than a bare `max |z| < c`: for a genuinely null-like profile the
  detector's 10%-of-median ratio is itself noisy at small N, so the relative statement
  is the honest, robust one.
- The open questions (even-admission vs the squarefree exponent cap, Axis C keep-even
  and Axis D) remain with the pre-registered compute jobs in `notes/`.
