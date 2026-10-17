# Pre-registration — lag-k serial correlation of prime-gap digit-sum parity

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.*

Registered run 49 (2026-10-16 ~13:33Z), **before** any lag-k number is measured.
This is the follow-on to the lag-1 job `gap-parity-corr-v1` (ships as v0.7.0).
File the staged job (`prep/lagk/`) as `gap-parity-lagk-v1` ONLY after v0.7.0 ships
and corr-v1 has drained (one gated job at a time — never pile).

## Question

v0.6.2 (`GAP-PARITY.md`) measured the total parity excess of the prime gaps;
v0.7.0 (`GAP-CORR.md`, in flight) measures the lag-1 serial correlation
S_b^(1) = Σ_n χ_b(g_n)χ_b(g_{n+1}) against a marginal-preserving null. The open
question v0.7.0 leaves: **does the gap-parity sequence carry serial memory beyond
one step?** Define, for lag k ≥ 1, base b ∈ {2,4,6,8,10,16}, cutoff 10^6/10^7/10^8,

    S_b^(k) = Σ_{n=1}^{N−k} χ_b(g_n) χ_b(g_{n+k}),   χ_b(m) = (−1)^{s_b(m)},  N = n_gaps.

Measure k = 1,2,3,4,5. Signal is the connected deviation z_b^(k) of S_b^(k) from
the random-permutation null that fixes the exact multiset of χ values.

## Exact permutation-null mean (hand-derived — gates the empirical null)

Fix the multiset of the N values χ_b(g_n); Σχ = −E_b, Σχ² = N. For a uniform
random permutation, any lag-k adjacent pair (positions n, n+k with k ≥ 1) is a
uniform random ORDERED pair of DISTINCT positions, so
E[χ_{π(n)}χ_{π(n+k)}] = (Σ_{a≠b}χ_aχ_b)/(N(N−1)) = (E_b²−N)/(N(N−1)). There are
N−k such terms, hence

    **E_perm[S_b^(k)] = (N−k)·(E_b²−N)/(N(N−1)).**

Check: k=1 gives (N−1)(E²−N)/(N(N−1)) = (E²−N)/N, exactly the lag-1 closed form
proved and shipped in v0.7.0. The staged `test_gap_lagk.py` confirms this general
form by brute force over EVERY permutation of small ±1 multisets, for k = 1,2,3.

## Hypotheses (decide from the z-table AFTER the numbers; no peeking)

- **D0 (null — most likely, the honest prior):** gap parities carry no serial
  memory beyond their histogram at any measured lag; |z_b^(k)| = O(1) for all b, k.
  (v0.7.0's lag-1 verdict is the strongest prior for this: if lag-1 is null, lag>1
  almost certainly is too.)
- **D1 (decaying memory):** lag-1 shows structure (|z| ≫ 1) that decays with k —
  the first evidence of finite-range serial memory in gap parity; report the decay
  profile and the base(s).
- **D2 (periodic / base-specific memory):** structure concentrated at a specific
  lag (e.g. k=2) or a specific base (base 8, the primes' valley base), not a
  monotone decay — report which lag/base and the sign.

Cross-validation requirement: the k=1 column of this job MUST reproduce v0.7.0's
S_b^(1) and z_b^(1) exactly (same sieve, same seed) — a built-in regression gate.

## Gate (job exits non-zero unless all pass — a failure is a bug, not a result)

1. π(x) = 78498 / 664579 / 5761455 at 10^6/10^7/10^8.
2. Telescoping: Σ g_n = p_M − 2 at each cutoff.
3. Single odd gap (g_1 = 1).
4. Half-gap Thue–Morse identity χ_2(g) = χ_2(g/2) for even g.
5. Two-routine digit guard, first 200 000 gaps, every base.
6. Reproduce v0.6.2 excess E_b^gap at each cutoff exactly (ties to the verified result).
7. **Reproduce v0.7.0 lag-1 exactly:** S_b^(1) equals corr-v1's value per base/cutoff.
8. Null-mean identity: empirical shuffle mean of S_b^(k) within 5·std/√K of
   (N−k)(E_b²−N)/(N(N−1)) for every measured k at 10^6 and 10^7.

Job must print `method gate PASSED` AND `null-mean gate PASSED` before any new number.

## Honest limits (stated in advance)

- Lags k ≤ 5 only; deeper memory not probed. Measured to 10^8 (5.76 M gaps);
  K-shuffle null (2000/800/200). Single machine, single sieve. The reproduce-v0.6.2,
  reproduce-v0.7.0 and null-mean identities are the internal cross-checks that stand
  in for an external anchor. A pre-registered D0 (no memory beyond histogram at any
  lag) is a real negative result and will be published as such.

— MV2, Munim, Inc.
