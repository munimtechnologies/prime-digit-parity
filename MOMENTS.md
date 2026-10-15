# MOMENTS.md — The base-2 prime digit-sum bias as a four-point function of base-16 parities: the complete 15-partition decomposition, and a derivation of the purity boundary

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. Code is
Apache-2.0; this text and the numbers in it are CC BY 4.0. This document
continues the study in `HIERARCHY.md` and `CORRELATION-RESULTS.md`.*

## Summary

`HIERARCHY.md` showed that the single base-2 (Thue–Morse) prime digit-sum bias
is, for every binary block width `b`, a `b`-point joint correlation of shifted
base-`2^b` parities, and *stated* a sharp boundary: that correlation is a pure
top-order connected cumulant only for `b ∈ {2, 3}`, and for `b ≥ 4` lower
pairings enter at leading order too. It did not say **why** the boundary sits
exactly there.

This document closes that gap. For the first width past the boundary, `b = 4`
(base 16), it writes the **complete, exact moment→cumulant expansion** — all 15
set partitions of the four parities — and then shows that the purity boundary
is nothing more than an elementary fact about those partitions: **a block of
size 1 carries a factor equal to a single-parity marginal mean, which the prime
data makes small, so only the partitions whose blocks all have size ≥ 2
survive at leading order.** The number of such "singleton-free" partitions is
`1` for `b = 2`, `1` for `b = 3`, and `4` for `b = 4` — which is precisely
"pure for `b ∈ {2, 3}`, not pure for `b ≥ 4`."

Everything here is an exact identity, checked from scratch by
`test_moments.py` over the real distribution on primes below `10^6` (no compute
queue, no floating point, `fractions.Fraction` throughout).

## 1. The character factorisation (recap, exact)

For integer `n ≥ 0`, with `χ_q(n) = (-1)^{s_q(n)}` the base-`q` digit-sum parity
character (values in `{+1, -1}`), the `b = 4` rung of the tower is the pointwise
identity

```
    χ_2(n) = χ_16(n) · χ_16(n>>1) · χ_16(n>>2) · χ_16(n>>3).
```

Each binary position of `n` lands in exactly one of the four residue classes
mod 4, and `χ_16(n>>r)` reads the bits at positions `≡ r (mod 4)`; the product
reads every bit, which is `χ_2(n)`. `test_moments.py` checks this for every
`n` in a dense low range and for every prime below `10^6`.

Write the four parities, as functions on the primes, as

```
    A = χ_16(n),  B = χ_16(n>>1),  C = χ_16(n>>2),  D = χ_16(n>>3).
```

Following the repo convention (`CORRELATION-RESULTS.md`), the base-2 excess over
the first `N` primes is `E_2 = #odd − #even = −Σ_p χ_2(p) = −N·E[ABCD]`, so the
four-point correlation `E[ABCD] = Σ_p A·B·C·D / N = −E_2/N` carries the same
anomaly up to this fixed sign (`E_2(10^6) = 4764`, so `N·E[ABCD] = −4764`).

## 2. The complete moment→cumulant expansion for four parities

Let `m_S = E[∏_{X∈S} X]` be the joint moments and `κ_S` the joint cumulants,
with singletons `κ_X = m_X = E[X] =: μ_X`. For four variables the
moment–cumulant relation runs over the 15 set partitions of `{A,B,C,D}`,
grouped by type:

```
  E[ABCD] =  κ_ABCD                                                  (1)   type 4
           + κ_ABC μ_D + κ_ABD μ_C + κ_ACD μ_B + κ_BCD μ_A           (4)   type 3+1
           + κ_AB κ_CD + κ_AC κ_BD + κ_AD κ_BC                       (3)   type 2+2
           + κ_AB μ_C μ_D + κ_AC μ_B μ_D + κ_AD μ_B μ_C
           + κ_BC μ_A μ_D + κ_BD μ_A μ_C + κ_CD μ_A μ_B              (6)   type 2+1+1
           + μ_A μ_B μ_C μ_D                                         (1)   type 1+1+1+1
```

`1 + 4 + 3 + 6 + 1 = 15 = B_4` (the fourth Bell number). The cumulants are
obtained from the moments by Möbius inversion on the partition lattice,
`κ_S = Σ_π (-1)^{|π|-1} (|π|-1)! ∏_{B∈π} m_B`. `test_moments.py` computes all
cumulants from the 16 measured sign-pattern counts on primes below `10^6`,
reconstructs `E[ABCD]` through the 15-partition sum above, and asserts it
equals the directly measured moment **exactly** — and that
`N · E[ABCD] = −4764`, i.e. the published base-2 excess `E_2 = 4764` at `10^6`.
The `b = 2`
(base-4 pair) and `b = 3` (base-8 triple) sub-tuples are reconstructed the same
way as a cross-check.

## 3. Singleton suppression, and why the boundary is at b ∈ {2, 3}

A partition block of size 1, say `{X}`, contributes the factor `κ_X = μ_X`, the
mean of one shifted base-16 parity. The prime data makes these marginals small:
the published base-16 anchor gives `|μ_A| = E_16(10^6)/N = 244/78498 ≈ 0.0031`
(the test pins `−N·μ_A = 244`), and the three shifted marginals `μ_B, μ_C, μ_D`
— excesses of the digit-sum parity of `⌊p/2^r⌋` over the primes, quantities not
previously tabulated, measured here by the test — are of the same small order,
consistent with Mauduit–Rivat balance. A partition with `s` singleton blocks
therefore carries a factor of order `μ^s`, suppressed relative to the leading
terms.

So at leading order only the **singleton-free** partitions survive — those whose
every block has size ≥ 2. Counting them is elementary:

| `b` | singleton-free partitions of a `b`-set | leading structure |
|----|----|----|
| 2 | `{AB}` — **1** | `E[AB] ≈ κ_AB` → pure pairwise cumulant |
| 3 | `{ABC}` — **1** | `E[ABC] ≈ κ_ABC` → pure 3rd cumulant (the CORRELATION result) |
| 4 | `{ABCD}` and the three 2+2 pairings — **4** | `κ_ABCD + κ_AB κ_CD + κ_AC κ_BD + κ_AD κ_BC` → **not** pure |
| 5 | `{ABCDE}` and the ten 3+2 splits — **11** | `κ_5` + `κ_3·κ_2` terms → not pure |

A set of size `b` can be split into blocks all of size ≥ 2 by some partition
other than the whole set **iff** `b ≥ 4` (you need room for a second block of
size ≥ 2; `b = 2, 3` leave none). That single combinatorial statement **is** the
purity boundary of `HIERARCHY.md`: the top cumulant stands alone at leading
order exactly for `b ∈ {2, 3}`, and the first thing to break it at `b = 4` is
the three `κ_2·κ_2` pairings — recovering, and now explaining, the earlier
pre-registration. `test_moments.py` verifies the counts `1, 1, 4` directly and
that the three surviving non-trivial `b = 4` partitions are exactly the 2+2
pairings.

This also re-reads the confirmed `b = 3` result (`CORRELATION-RESULTS.md`: the
third base-8 cumulant carried 99.84% of the base-2 excess at `10^8`) as a
special case: at `b = 3` there is *no* singleton-free competitor to the top
cumulant, so near-total dominance is the generic expectation, not a surprise —
whereas at `b = 4` the pairings are structurally guaranteed a leading-order
share.

## 4. What is verified, and the honest limits

Checked exactly by `test_moments.py` (stdlib, `Fraction`, from scratch):

* the pointwise factorisation `χ_2 = A·B·C·D` (all `n` in a dense range; all
  primes `< 10^6`);
* the full 15-partition moment→cumulant round trip over the real prime
  distribution at `10^6`, reproducing `N·E[ABCD] = −4764` (the published base-2
  excess `E_2 = 4764`);
* the independent-distribution control: on four independent parities every
  cumulant spanning ≥ 2 variables vanishes and `E[ABCD] = ∏ μ`, isolating the
  cumulant code from the round trip;
* the singleton-free partition counts `1, 1, 4` for `b = 2, 3, 4`;
* the anchors `π(10^6)=78498`, `E_2=4764`, `E_16=244`.

Honest limits. (i) The *suppression* of singleton partitions is an exact
algebraic factor (`= μ_X`) times a small, empirically measured marginal; the
claim that `μ_X → 0` fast enough to make the leading-order picture asymptotic
rests on Mauduit–Rivat balance, which bounds but does not by itself order these
finite-`N` pieces — the per-partition split at `10^8` for `b = 4` needs the base-16
shifted-parity counts, a sandbox job still in the queue (`prime-parity-pairwise-v1`).
(ii) This document is the exact *bookkeeping and derivation*; the measured
`b = 4` partition magnitudes at `10^6` are computed by the test but not frozen
into prose here, and are a finite-`N` snapshot, not an asymptotic statement.
(iii) The decomposition is a mathematical identity for any distribution; its
content for primes is which partitions the prime data loads — settled only in
the small-`b` limit and for the marginals' order, as stated.

## 5. Prior art

Gupta–Sridhar (arXiv:1812.11841) report the base-10 excess at one cutoff and its
cross-base persistence as unexplained; Mauduit–Rivat (Annals 2010) and Drmota
and collaborators give limiting densities with upper-bound error terms. Within
this project, `REDUCE.md` gave the identity without cumulants, `CORRELATION.md`
the single `b = 3` rung, and `HIERARCHY.md` the tower and the *statement* of the
purity boundary. The complete 15-partition expansion for `b = 4` and the
reduction of the purity boundary to the count of singleton-free set partitions
are, to a prior-art search across the literature, arXiv, OEIS and GitHub, new
here, and are machine-checked as exact `Fraction` identities.

— MV2, Munim, Inc.
