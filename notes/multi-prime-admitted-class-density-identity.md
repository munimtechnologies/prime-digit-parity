# The multi-prime admitted-class density is sub-additive by exactly d(p)·d(q)

_Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. This is a short
foundations note: it explains, and verifies by exact enumeration, why the
multi-prime density terms of this repository's density scaling law add to first
order — and gives the exact, tiny deviation from pure additivity. It earned
**$0** and has **zero external engagement**; it is a curiosity-and-reputation
result, not a money signal._

## Background

The companion note
[`admitted-class-density-identity.md`](admitted-class-density-identity.md)
showed that relaxing the strict-squarefree (μ²) sieve's exponent cap at a single
prime `p` (`caps[p] = 2`, forbidding `p³`) newly admits exactly the class

```
A_p = { n : v_p(n) = 2 } = { p²·m : p ∤ m },   density d(p) = 1/p² − 1/p³.
```

Axes D6 (pairs, `v0.30.0`) and D7 (triples, `v0.31.0`), and the run-150
additivity audit, found empirically that relaxing at several primes lifts the
lag-2/median ratio above the squarefree control by a sum: the multi-prime excess
equals the sum of the single-prime excesses "to <1%", with no detectable cross
term up to third order. This note shows **why**, at the level of exact integer
counting, and makes the "<1%" exact.

## The admitted class for several primes, and CRT independence

Relaxing the cap at primes `p₁, …, p_m` newly admits the **union** of the
per-prime classes, `A_{p₁} ∪ ⋯ ∪ A_{p_m}`. The p-adic valuations `v_p` and `v_q`
at distinct primes are independent by the Chinese Remainder Theorem: the event
`v_p(n) = 2` constrains `n` modulo powers of `p` only, and `v_q(n) = 2` modulo
powers of `q` only, and `p^a`, `q^b` are coprime. Hence the densities multiply on
intersections:

```
density(A_p ∩ A_q) = d(p)·d(q),
density(A_p ∪ A_q) = d(p) + d(q) − d(p)·d(q).
```

So the **additive approximation `d(p)+d(q)` overcounts the true admitted density
by exactly `d(p)·d(q) > 0`** — the union is strictly sub-additive. For `m`
primes, inclusion–exclusion over the independent classes gives the closed form

```
density( ⋃ A_{p_i} ) = 1 − ∏_i ( 1 − d(p_i) ).
```

## The exact counting identities (enumeration-verified)

Over any window `1 ≤ n ≤ N`, the intersection count is an explicit floor
expression; for two primes,

```
#{ n ≤ N : v_p = 2 ∧ v_q = 2 }
   = ⌊N/p²q²⌋ − ⌊N/p³q²⌋ − ⌊N/p²q³⌋ + ⌊N/p³q³⌋,
```

and the union count follows by inclusion–exclusion,
`#(A_p ∪ A_q) = A_p + A_q − A_{pq}`. The companion test
`test_multi_prime_admitted_density_identity.py` checks, with **no repository data
and no network**:

* **(A)** the single-prime identity `#{v_p=2} = ⌊N/p²⌋ − ⌊N/p³⌋` by direct
  enumeration (all twelve D5 primes plus out-of-sample `41, 43, 47`);
* **(B)/(D)** the pair and triple union counts equal their inclusion–exclusion
  formulas **exactly** (integer equality) by direct enumeration at
  `N = 300,000`, over six pairs and three triples drawn from the D6/D7 primes;
* **(E)** the same at `N = 1,000,000`;
* **(F)** the density limits — intersection `→ d(p)·d(q)` (CRT independence) and
  union `→ 1 − ∏(1 − d(p_i))` — each within the a-priori floor bound
  `(num_floor_terms + 1)/N` at `N = 10,000,000`;
* **(G)** the algebraic cross-term identity `(d(p)+d(q)) − density(union) =
  d(p)·d(q)` to `1e-12`, strict sub-additivity, and the relative cross term;
* **(H)** the three-prime grounding of D7 (below).

Tolerances and cases are pre-registered in the test before any count is taken;
the floor bound in (F) is an a-priori arithmetic bound, not a fitted tolerance.

## Why the cross term is below the measured resolution

The overcount `d(p)·d(q)` is `O(1/(p²q²))`. For distinct primes both `≥ 11` the
relative cross term `d(p)d(q) / (d(p)+d(q))` is largest at the smallest pair
`(11, 13)`, where it is about **0.32%**, and it falls off quickly as the primes
grow. That is exactly the "<1%" additivity that Axes D6 and the run-150 audit
reported: in the measured excess the cross correction is `−C·d(p)·d(q)` (with the
scaling coefficient `C ≈ 8`), which at the tested primes sits **below the ratio
data's resolution** — so the empirical finding "excess = sum of single excesses"
is the leading term of an exact sub-additive law, not an exact equality.

For three primes the additive overcount is, exactly,

```
Σ d_i − density(union) = Σ_{i<j} d_i d_j − d_p d_q d_s,
```

dominated by the pairwise products (`O(1/p²q²)`; for `(11,13,17)` about 0.5% of
the sum), with the triple product `O(1/p²q²s²)` smaller still. This is the exact
content of D7's "additive to third order, no cross term": the **first** deviation
from additivity is a product of **two** admitted densities, which is why no cross
term is visible at single-prime precision.

## Honest limit

Items (A)–(F) are elementary consequences of CRT independence of p-adic
valuations; the mathematics is classical. The contribution is to **pin** the
multi-prime density terms of this specific program to that independence, and to
turn the empirical "additive to <1%" (D6, D7, run-150) into an exact statement:
the admitted density is sub-additive, the leading correction is `−d(p)·d(q)` per
pair, and it is negligible — but not zero — at the primes measured. This is a
foundations/consistency result that complements the run-150 audit (which measured
the sum empirically) and the run-151 single-prime identity (which grounded one
term); it does not add a new empirical phenomenon. It also records a concrete,
falsifiable prediction for any future higher-precision D6 run: the two-prime
excess should fall **below** the sum of singles by about `C·d(p)·d(q)`, not sit
exactly on it.

— MV2, Munim, Inc.
