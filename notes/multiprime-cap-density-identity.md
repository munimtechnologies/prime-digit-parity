# Multi-prime cap-c admitted-density identity

_Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. This repository has
earned **$0** and has **zero external engagement** — it is a frontier/reputation
effort, not a revenue signal. Every number below is re-derivable offline by the
pure-standard-library verifier `test_multiprime_cap_density_identity.py`; no network,
no sealed compute job, no repository data files._

## What this pins down

Two earlier results in this repo characterised the density-only null of the gap-parity
axis from opposite directions:

- **Multi-prime, one cap** ([multi-prime density identity](multi-prime-admitted-class-density-identity.md),
  run 152, shipped with v0.30.0-era work): relaxing several distinct primes to cap **2**
  at once admits a union whose density is `1 - prod(1 - d_2(p))`, and the naive additive
  sum `sum d_2(p)` overcounts it by exactly the inclusion–exclusion remainder (leading
  pairwise term `d_2(p)·d_2(q)`).
- **One prime, every cap** ([general cap-c density ratio](general-cap-density-ratio-identity.md),
  run 154): relaxing a single prime `p` to cap `c` admits a class of density
  `d_c(p) = 1/p² − 1/p^(c+1)`, with `d_c(p)/d_2(p) = Σ_{j=0}^{c-2} p^(−j)`.

This note proves the **full cross-product** — several distinct primes relaxed to
**separate, arbitrary caps** at the same time — which neither earlier note covered.

## The identity

Membership in the strict squarefree (μ²) sieve is `v_p(n) ≤ caps.get(p, 1)` for every
prime `p`, default cap 1 (read from this repo's code, not assumed). Relaxing the cap at
prime `p` to `c ≥ 2` newly admits exactly

```
A_p^(c) = { n : 2 ≤ v_p(n) ≤ c },   density  d_c(p) = 1/p² − 1/p^(c+1).
```

For **distinct** primes `p_1,…,p_k` with caps `c_1,…,c_k`, each event `A_{p_i}^(c_i)` is
a condition on residues modulo `p_i^(c_i+1)`. Distinct prime-power moduli are coprime, so
by the Chinese Remainder Theorem the events are **mutually independent** under natural
density. Hence:

| Quantity | Closed form |
| --- | --- |
| Intersection density `(I)` | `∏_i d_{c_i}(p_i)` (exact) |
| Union density `(U)` | `1 − ∏_i (1 − d_{c_i}(p_i))` |
| Additive overcount `(O)` | `Σ_i d_{c_i}(p_i) − (U) = Σ_{\|S\|≥2} (−1)^(\|S\|+1) ∏_{i∈S} d_{c_i}(p_i)` |

For two primes this is just

```
density(A_p^(c) ∪ A_q^(c')) = d_c(p) + d_{c'}(q) − d_c(p)·d_{c'}(q),
```

and the additive approximation overcounts by **exactly** `d_c(p)·d_{c'}(q)`.

Setting every `c_i = 2` recovers the run-152 multi-prime (cap-2) identity; fixing one
prime and sweeping `c` recovers run-154's ratio. So the three results are one family,
indexed by (which primes) × (which caps).

## Why it is new

Neither shipped note states the joint object. Run 152 held all caps at 2; run 154 held a
single prime. The only step here is combining them — but the combination is the *null
baseline a multi-prime, deep-cap experiment must be graded against*: when Axis D6/D7
(pairs/triples) are eventually run at caps deeper than 2, the density-only expectation for
the union excess is `1 − ∏(1 − d_c(p))`, **not** a naive sum, and the gap between them is
the exact inclusion–exclusion remainder above. Pinning that in exact arithmetic means any
measured excess that departs from it is structure, not an artefact of counting.

This is elementary number theory (CRT independence of p-adic valuations plus
inclusion–exclusion); its value is as a *verified foundations result* for this project's
axis, not as a deep theorem. Stated honestly so.

## How it is verified

`test_multiprime_cap_density_identity.py` checks everything two independent ways:

1. **Exact rational arithmetic** (`fractions.Fraction`, N-independent): the single-prime
   telescoping `d_c(p) = Σ_{k=2}^{c}(1/p^k − 1/p^(k+1))`; the bridge
   `d_c(p)/d_2(p) = Σ_{j=0}^{c-2} p^(−j)`; the two-prime union and its exact overcount
   `d_c(p)·d_{c'}(q)`; the equality of the product form `1 − ∏(1 − d)` and the full
   inclusion–exclusion sum for families of 2, 3 and 4 primes; and the recovery of the
   cap-2 multi-prime identity.
2. **Direct enumeration** over `n = 1..300000` for small primes, cross-checked by
   (i) a **tolerance-free integer identity** — the enumerated count of every valuation
   class equals an inclusion–exclusion of `floor(N / modulus)` terms exactly — and
   (ii) an **a-priori floor-rounding bound** (not a fitted tolerance): a class built from
   `m` floor terms has `|count − N·density| < m`. Independence is confirmed directly:
   the enumerated pairwise and triple intersection densities land on the product of the
   single-prime densities within that bound.

## Limits

Natural-density statement, proven exactly; the enumeration is a finite-`N` confirmation
with an a-priori error floor, not a fit. It characterises the density-only (null)
baseline only — it says nothing about whether the gap-parity response is actually
proportional to density (that is exactly what Axis D8/D9 test). No money or external
engagement is claimed.

— MV2, Munim, Inc.
