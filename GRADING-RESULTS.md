# Measured singleton grading of the base-2 prime digit-parity excess at 10^8

*Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Part of
munimtechnologies/prime-digit-parity. Code: Apache-2.0. Text and data:
CC BY 4.0.*

[GRADING.md](GRADING.md) (v0.4.11) predicted that the base-2^b prime
digit-sum parity correlation splits into suppression orders indexed by the
number of **singleton** blocks *k*, with exactly `P(b,k) = C(b,k)*A000296(b-k)`
terms at order *k*. [CORRELATION-RESULTS.md](CORRELATION-RESULTS.md) (v0.4.8)
measured, for base 8 (b = 3) over the 5,761,455 primes below 10^8, that the
connected 3-point cumulant alone carries "99.84 %" of the base-2 excess, with a
single lumped remainder of "114". This note closes the loop between the two:
it reads the **full graded table** off the same eight octal-triple counts, so
every suppression order *k = 0, 1, 2, 3* gets a measured magnitude, and the
order-*k* term count matches `P(3,k)` exactly.

## Setup (recap)

For a prime *p*, let the three base-8 digit-sum parity columns be

```
A = chi_8(p),   B = chi_8(p>>1),   C = chi_8(p>>2),      each in {+1,-1}.
```

Because a base-8 digit contributes only its lowest bit to the digit-sum parity,
column *r* is the parity of the bits of *p* at positions ≡ *r* (mod 3), and the
tower reduction gives `chi_2(p) = A*B*C`. Hence the base-2 (Thue-Morse) excess
is `E_2 = -sum_p A*B*C`. The eight counts of the sign pattern (A,B,C) over the
primes below 10^8 (index `i = 4a+2b+c`, bit 0 = even/`+1`):

```
693845  736247  734665  716635  736921  719439  715554  708149   (sum = 5761455)
```

These gate on the published anchors: `-sum(A*B*C) = 70509 = E_2(10^8)` and
`-sum(A) = -1329 = E_8(10^8)`.

## The moment-cumulant (linked-cluster) expansion

The whole-correlation mean is the sum over set partitions of {A,B,C} of the
product of the connected cumulant of each block:

```
E[ABC] = kappa(A,B,C)                                   # 1 partition, 0 singletons
       + <A>kappa(B,C) + <B>kappa(A,C) + <C>kappa(A,B)  # 3 partitions, 1 singleton
       + <A><B><C>                                      # 1 partition, 3 singletons
```

Grading by the number of singleton blocks *k* is exactly GRADING.md's row:
`P(3,0)=1, P(3,1)=3, P(3,2)=0, P(3,3)=1` (sum = Bell(3) = 5). Each singleton is
a bare single-column mean `<.>`, and those means are tiny
(`<A>=2.31e-4, <B>=1.99e-3, <C>=8.94e-5`), so each extra singleton suppresses a
term by one factor of a small marginal — the "small-marginal expansion" of
GRADING.md, now with numbers.

## Measured table (all exact rationals; see `data/base8_graded.csv`)

| k | # terms = P(3,k) | blocks | contribution to E_2 = 70509 | share |
|---|---|---|---|---|
| 0 | 1 | {A,B,C} (connected 3-point) | **+70394.68** | **99.8379 %** |
| 1 | 3 | one singleton × a connected pair | +114.32 | 0.1621 % |
| 2 | 0 | — (P(3,2) = 0, structurally) | 0 | 0 % |
| 3 | 1 | {A}{B}{C} (fully disconnected) | −0.000236 | ≈ 0 % |
| | 5 = Bell(3) | total | 70509 | 100 % |

The three order-1 terms, individually (contribution to E_2):

| singleton | paired columns | contribution |
|---|---|---|
| B (largest marginal) | A·C | +97.89 |
| A | B·C | +11.62 |
| C (smallest marginal) | A·B | +4.82 |

So the famous base-10/base-2 bias at 10^8 is, to 99.84 %, a single **irreducible
3-point correlation** among the three Thue-Morse sub-columns of the primes; the
entire reducible remainder is 0.16 %, and within it the order where "two columns
are marginal and one pairs" **cannot occur** (k = 2 is empty because
`A000296(1) = 0`, i.e. you cannot leave exactly one element un-singletoned). The
fully-disconnected k = 3 piece is ~2×10^-4, numerically negligible.

## Why this is new

GRADING.md gives the *prediction* (the `P(b,k)` structure); CORRELATION-RESULTS.md
gives one *lumped* measured number (κ_3 = 99.84 %, remainder 114). Neither splits
the measured remainder by singleton order, reports the three individual order-1
terms, or notes the structural vanishing of the k = 2 order in the data. A prior-
art search (Thue-Morse / digit-sum-parity / set-partition literature, and this
project's own releases) turns up no measured singleton-graded decomposition of a
prime digit-sum parity correlation. The result is a measured validation of
GRADING.md's leading-order law (the k = 0 singleton-free term dominates) with the
exact rational magnitudes of every lower-order correction.

## Honest limits

- **One base, one cutoff.** This is b = 3 (base 8) at 10^8 only. The analogous
  tables for b = 2 (base 4) and b = 4, 5 need the per-column counts from the
  queued sandbox jobs (`prime-parity-pairwise-v1`, `evenbase-factors-v1`), which
  have not returned; when they do, the same code grades them.
- **No mechanism.** We measure that the connected 3-point cumulant dominates; we
  do not prove *why* it is non-zero along the primes (that is the Mauduit-Rivat /
  Drmota circle of questions, where only asymptotic balance is known).
- **Drift with x is open.** The 99.84 % share is reported at a single x; how it
  moves with the cutoff is not settled here.

Everything above is recomputed from the eight counts as exact `Fraction`s in
`test_grading_results.py`, which also re-sieves the primes below 10^6 from
scratch, rebuilds the three columns from the bit positions, and checks they
reproduce the independent anchors `E_2(10^6) = 4764` and `E_8(10^6) = -510`
before any 10^8 number is asserted.
