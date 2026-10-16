# Lag-1 serial correlation of prime-gap digit-sum parity

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc.*

## The question

v0.6.2 (`GAP-PARITY.md`) measured the **total** digit-sum-parity excess
E_b^gap = −Σ_n χ_b(g_n) of the prime gaps g_n = p_{n+1} − p_n, and showed it is
large in every base. But that total is **fully fixed by the gap-size histogram**:
E_b^gap = −Σ_g count(g)·χ_b(g). It says nothing about the *order* in which the
gaps appear.

Here we ask the question v0.6.2 leaves open — the Lemke Oliver–Soundararajan
analogue for gap parity: are consecutive gap parities χ_b(g_n), χ_b(g_{n+1})
**serially correlated beyond what the fixed multiset of parities forces?** For
bases b ∈ {2,4,6,8,10,16} and cutoffs 10^6/10^7/10^8 define the signed
adjacent-pair sum (same-parity minus different-parity count),

    S_b = Σ_{n=1}^{N−1} χ_b(g_n)·χ_b(g_{n+1}),   χ_b(m) = (−1)^{s_b(m)},  N = n_gaps.

The signal is the *connected* deviation of S_b from a marginal-preserving
(random-permutation) null that fixes the exact multiset of χ values.

## Why this is new

- Lemke Oliver–Soundararajan (PNAS 2016, arXiv:1604.02443) and its analogues
  (sums of two squares: David–Devin–Nam–Schlitt, arXiv:2105.05048) study serial
  bias in the **residue classes** of consecutive primes / sums, not the
  digit-sum parity of the gaps.
- Prime-gap serial-dependence work uses windowed rank, spectral analysis, the
  Cramér i.i.d. model and Hardy–Littlewood correlations — none on χ_b(g_n).
- Correlation-of-Thue–Morse work (decomposing the sum-of-digits correlation
  measure, arXiv:2411.07779; Rudin–Shapiro along primes, arXiv:2306.10629;
  Spiegelhofer's gaps *within* the Thue–Morse word; the Baake–Grimm Thue–Morse
  autocorrelation *measure*, math/0610791 / S0019357723000186, a singular-
  continuous spectrum) treats the TM/RS sequence itself or its evaluation
  *along* a subsequence of ℕ — not the serial correlation of a derived
  prime-gap sequence χ_b(g_n).
- This repo's own v0.6.2 measured only the *total* excess (the histogram
  readout); the lag-1 serial structure is the one thing that total cannot see.

Six fresh prior-art searches (two at run 47, two at run 48, two at run 50,
2026-10-16) and one more at this shipping run found no source computing the
serial correlation of χ_b(g_n). → novel.

## Exact permutation-null mean (hand-derived — gates the empirical null)

Fix the multiset of the N values χ_b(g_n); then Σχ = −E_b^gap and Σχ² = N.
For a uniform random permutation, summing over ordered distinct position pairs
gives Σ_{a≠b} χ_aχ_b = (Σχ)² − Σχ² = E_b² − N over N(N−1) ordered pairs, and each
of the N−1 adjacent pairs is a uniform random distinct ordered pair, so

    **E_perm[S_b] = (N−1)·(E_b² − N)/(N(N−1)) = (E_b² − N)/N.**

`null_moments.py` (shipped in v0.6.3) and `test_gap_corr.py` confirm this closed
form by brute force over *every* permutation of small ±1 multisets. In the sealed
run the empirical K-shuffle mean matched it to within 5·std/√K at 10^6 and 10^7 —
a hard gate (passed).

## Anchors (gate — all must pass or the number is a bug, not a result)

1. **π(x)** = 78498 / 664579 / 5761455 at 10^6/10^7/10^8.
2. **Telescoping:** Σ g_n = p_M − 2 at each cutoff.
3. **Single odd gap:** exactly one odd gap (g_1 = 1).
4. **Half-gap Thue–Morse identity:** χ_2(g) = χ_2(g/2) for even g.
5. **Two-routine digit guard** for the first 200 000 gaps, every base.
6. **Reproduce v0.6.2 exactly:** E_b^gap at each cutoff equals the published
   values (e.g. E_8@10^8 = −437046) — ties this job to the verified result.
7. **Null-mean identity:** empirical shuffle mean of S_b within 5·std/√K of
   (E_b²−N)/N at 10^6 and 10^7.

The sealed CPU-sandbox run printed, before any new number:
> `method gate PASSED (pi, telescoping, single odd gap, half-gap identity,`
> `two-routine guard, v0.6.2 excess reproduced).`
> `null-mean gate PASSED (empirical shuffle mean == (E_b^2-N)/N).`

## Result

Gate: **PASSED** (both gates, from the sealed CPU-sandbox run log).

S_b, the closed-form permutation mean E_perm = (E_b²−N)/N, and the connected
signal z_b = (S_b − shuffle_mean)/shuffle_std, by base and cutoff:

| base | S_b @10^8 | E_perm @10^8 | z @10^6 | z @10^7 | z @10^8 |
|-----:|----------:|-------------:|--------:|--------:|--------:|
| 2  | −166 721 |   99 300.8 | −19.37 | −47.25 | −120.58 |
| 4  |  233 701 |  192 058.9 |  −1.23 |  +3.73 |  +18.10 |
| 6  |   12 329 |   20 245.1 |  −2.42 |  −5.11 |   −3.48 |
| 8  |  −45 597 |   33 152.0 |  −4.34 | −11.96 |  −33.25 |
| 10 |   67 285 |   76 936.8 |  −3.09 |  −5.83 |   −4.39 |
| 16 |  508 353 |  627 829.4 |  −6.88 | −20.50 |  −58.15 |

(n_gaps = 78 497 / 664 578 / 5 761 454; K shuffles = 2000 / 800 / 200. Full
per-cutoff S_b, E_perm, shuffle_mean, shuffle_std and K are in
`data/gap_parity_corr.csv`.)

## Which hypothesis the data support (decided from the z-column; prereg run 47)

Registered **before** the numbers (`notes/gap-parity-correlation-prereg.md`):

- **C0 (null):** gaps serially independent given their parity marginal;
  |z_b| = O(1) for all b.
- **C1 (real lag-1 structure):** some base shows |z_b| ≫ 1 — the first such
  finding; sign tells whether equal parities cluster (+) or alternate (−).
- **C2 (base-8 special again):** only base 8 (the primes' valley base) shows
  serial structure.

**Verdict: C1, decisively — and it is *not* base-8-specific (so not C2).**
Consecutive prime-gap digit-sum parities carry strong lag-1 serial structure,
and the *sign is base-dependent*:

- **Bases 2, 8, 16 are strongly *anti-correlated* (z < 0): adjacent gaps tend to
  carry *opposite* digit-sum parity.** The effect strengthens monotonically with
  the cutoff — base 2 runs −19.4 → −47.2 → −120.6, base 16 −6.9 → −20.5 → −58.1,
  base 8 −4.3 → −12.0 → −33.3 — so it is a persistent, growing structure, not a
  fluctuation. Base 8 is in the middle of the pack, not special: bases 2 and 16
  are far stronger.
- **Base 4 is the lone *positively* correlated base (z > 0: equal parities
  cluster), and the effect *emerges at scale*** — not significant at 10^6
  (−1.2), then +3.7 at 10^7 and +18.1 at 10^8.
- **Bases 6 and 10 show only weak, non-monotone anti-correlation** (|z| ≲ 6),
  the borderline of the effect.

So prime gaps are *not* a parity-independent reshuffle of their own histogram:
the order in which the digit-sum parities of consecutive gaps appear carries a
real, large, base-structured lag-1 memory — anti-correlation in most bases,
clustering in base 4.

## Analytic null variance (exact cross-check, zero Monte-Carlo error)

Because the null is a uniform permutation of the fixed ±1 multiset, the empirical
shuffle std is an unbiased estimate of the *exact* permutation std, whose closed
form was derived and exhaustively verified in v0.6.3 (`NULL-MOMENTS.md`,
`null_moments.perm_var`). With s = Σχ = −E_b and every χ² = 1, the
cancellation-free form is

    **Var[S_b] = (N−1) − 2·E_b²/N + E_b⁴/(N²(N−1)) − 1/(N−1).**

(`test_gap_corr.py` brute-forces this against every permutation of small
multisets with exact `fractions.Fraction`.) The exact std and the empirical
K-shuffle std agree within Monte-Carlo error; e.g. base 2 at 10^6 gives analytic
std 277.6 vs empirical 273.9 (0.8σ, K=2000), and base 2 at 10^7 gives 804.1 vs
783.7 (~1σ, K=800). Reporting the analytic z = (S_b − E_perm)/√Var removes the
Monte-Carlo error in the denominator (material only for the most imbalanced
bases, e.g. base 16, where the K=800 shuffle std carries a few-percent error);
it leaves every verdict above unchanged. The full analytic-z column is a
deterministic post-process of the shipped CSV via `null_moments.perm_var`.

## Honest limits

- Lag-1 only; k > 1 (does the anti-correlation persist or oscillate at longer
  range?) is deferred to the lag-k write-up. The sign interpretation is for
  adjacent gaps; it does not by itself explain *why* base 4 clusters while base 2
  alternates.
- Measured to 10^8 (5.76 M gaps); K-shuffle null (2000/800/200 shuffles) on a
  single machine, single sieve. The reproduce-v0.6.2 and null-mean identities are
  the internal cross-checks that stand in for an external anchor; the exact
  permutation variance (above) is the analytic cross-check on the significance.
- This is an empirical pattern in a finite range, not a theorem; the growth of
  |z| with the cutoff is suggestive but unproven asymptotically.

## Reproduce

```
python -m unittest -v test_gap_corr   # anchors + brute-force null mean & variance, stdlib only
```

Full measurement (sealed CPU sandbox, needs numpy):
`gap_corr_job.py --out-dir out --limit 100000000`. Data:
`data/gap_parity_corr.csv` (CC BY 4.0). Code Apache-2.0.

— MV2, Munim, Inc.
