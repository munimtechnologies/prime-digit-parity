# Preregistration: serial (lag-1) correlation of prime-gap digit-sum parity

**Author:** MV2 of Munim, Inc. Work by Munim, Inc., produced with AI.
**Status:** preregistered design + exact null-mean derivation. NOT yet measured.
**Written:** run 47 (2026-10-16), discovery field. Follow-up to v0.6.2 (GAP-PARITY.md).

## The gap left by v0.6.2
v0.6.2 established that the **total** digit-sum-parity excess E_b^gap of prime gaps
is large in every base. But that total is **fully determined by the gap-size
histogram**: E_b^gap = −Σ_g (count of gap size g)·χ_b(g). It says nothing about the
*order* of the gaps. The genuinely open question — the Lemke Oliver–Soundararajan
analogue for gap parity — is whether χ_b(g_n) and χ_b(g_{n+1}) are **serially
correlated beyond what the fixed multiset of gap parities forces**.

## Prior art (searched 2026-10-16)
Prime-gap serial-dependence work uses windowed rank (XuanJi-ISA issue #174),
spectral analysis, the Cramér i.i.d. model, and Hardy–Littlewood correlations —
none on the **digit-sum parity** of gaps. Lemke Oliver–Soundararajan (1604.02443)
studies serial bias in the *last digits* of consecutive primes, not gap digit-sum
parity. → the observable below is novel.

## Observable
χ_b(m) = (−1)^{s_b(m)} ∈ {+1,−1}. Over the N = n_gaps gaps below a cutoff define
the signed adjacent-pair sum (same-parity minus different-parity count):

    S_b = Σ_{n=1}^{N−1} χ_b(g_n)·χ_b(g_{n+1})  = (#same-parity pairs) − (#different).

## Exact permutation-null mean (hand-derived — gates the empirical null)
Fix the multiset of the N values χ_b(g_n); Σχ = −E_b^gap, Σχ² = N. For a uniformly
random permutation, any two distinct positions have
E[χ_aχ_b] = ((Σχ)² − Σχ²)/(N(N−1)) = (E_b² − N)/(N(N−1)), and there are N−1
adjacent pairs, so

    **E_perm[S_b] = (E_b² − N)/N.**     (closed form; a hard gate on the shuffle mean)

The empirical shuffle (K marginal-preserving shuffles, fixed seed) must reproduce
this mean within a few standard errors. The **signal** is the connected deviation
Δ_b = S_b(actual) − E_perm[S_b], standardized z_b = (S_b − mean_shuffle)/std_shuffle.

## Hypotheses (decided BEFORE seeing numbers)
- **C0 (null — most likely):** gaps are serially independent given their parity
  marginal; |z_b| = O(1) for all b. Prime gaps are widely modeled as pseudorandom,
  so no lag-1 parity correlation beyond the multiset.
- **C1 (repulsion/attraction):** some base shows |z_b| ≫ 1 — a real sequential
  structure in gap parities (would be the first such finding), with sign telling
  whether equal parities cluster (+) or alternate (−).
- **C2 (base-8 special again):** only base 8 (the primes' valley base) shows serial
  structure — would tie the gap ordering back to the prime base-8 effect.

A clean C0 null is a real, publishable negative (gap parities carry no lag-1 memory
beyond their histogram).

## Gate (job exits non-zero unless all pass)
1. Original v0.6.2 anchors: π(x), telescoping Σg=p_M−2, single odd gap, half-gap
   Thue–Morse identity, two-routine digit guard.
2. **Reproduce v0.6.2 exactly:** E_b at each cutoff equals the published values
   (e.g. E_8@10^8 = −437046) — ties this job to the verified result.
3. **Null-mean identity:** empirical shuffle mean of S_b within 5·std_shuffle/√K of
   the closed form (E_b²−N)/N, at 10^6 and 10^7.

## Output
`gap_parity_corr.json`: per base/cutoff E_b, S_b, E_perm[S_b], mean/std shuffle,
z_b. Ships as v0.7.0 with a CC BY 4.0 CSV + self-contained test iff the gate passes
and a search confirms the correlation observable is still unclaimed.
