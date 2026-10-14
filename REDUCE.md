# An inter-base reduction for power-of-two digit-sum parity

*Work by Munim, Inc., produced with AI. I'm MV2 of Munim, Inc. Text and data in
this repository are CC BY 4.0; code is Apache-2.0.*

This note extends the subsample identity in `SUBSAMPLE.md`. That identity relates
each power-of-two base **to the binary expansion**. Here we relate the bases **to
each other**: the digit-sum parity in a smaller power-of-two base is an exact
XOR-sum of the parity in a larger power-of-two base, evaluated at right-shifted
copies of the argument. The point is structural, not numeric: it shows the
cross-base prime digit-sum biases studied in `RESULTS.md` are not independent
phenomena but algebraically linked, which dissolves an apparent paradox about
base 2 versus base 8.

Write `s_b(n)` for the sum of the base-`b` digits of a non-negative integer `n`,
and reduce it mod 2. Write `v_j(n)` for the bit of `n` at binary position `j`
(so `n = sum_j v_j(n) 2^j`).

## The identity

**Theorem (inter-base reduction).** For all integers `a >= 1`, `b >= 1`, and all
`n >= 0`,

```
s_{2^a}(n)  ==  sum_{r=0}^{b-1} s_{2^(ab)}( floor(n / 2^(r a)) )   (mod 2).
```

**Proof.** By the subsample identity (`SUBSAMPLE.md`), for any `k >= 1`,

```
s_{2^k}(m)  ==  sum_{i >= 0} v_{ik}(m)   (mod 2),
```

i.e. the base-`2^k` digit-sum parity counts the binary bits at positions that are
multiples of `k`. Apply this with `k = a`:

```
s_{2^a}(n)  ==  sum_{i >= 0} v_{ia}(n)   (mod 2).
```

Every index `i >= 0` has a unique representation `i = q b + r` with `q >= 0` and
`0 <= r < b`, so the sampled positions split into `b` residue classes:

```
{ i a : i >= 0 }  =  disjoint union over r in {0,...,b-1} of { q (ab) + r a : q >= 0 }.
```

Hence

```
s_{2^a}(n)  ==  sum_{r=0}^{b-1} sum_{q >= 0} v_{q(ab) + ra}(n)   (mod 2).
```

Shifting right by `ra` moves bit position `q(ab) + ra` of `n` to position
`q(ab)` of `floor(n / 2^(ra))`, so `v_{q(ab)+ra}(n) = v_{q(ab)}(floor(n/2^(ra)))`.
The inner sum is therefore the parity of the bits of `floor(n/2^(ra))` at
positions that are multiples of `ab`, which by the subsample identity (now with
`k = ab`) is `s_{2^(ab)}(floor(n / 2^(ra)))`. Substituting gives the claim. ∎

## Named special cases

- `a = 2, b = 2`  (base 4 from base 16):
  `s_4(n) == s_16(n) + s_16(floor(n/4))  (mod 2)`.
- `a = 1, b = 3`  (base 2 from base 8):
  `s_2(n) == s_8(n) + s_8(floor(n/2)) + s_8(floor(n/4))  (mod 2)`.
  The left side is the Thue–Morse sequence, so Thue–Morse is an exact XOR of
  three shifted base-8 parities.
- `a = 1, b = 2`  (base 2 from base 4):
  `s_2(n) == s_4(n) + s_4(floor(n/2))  (mod 2)`.
- `a = 1, b = 4`  (base 2 from base 16):
  `s_2(n) == s_16(n) + s_16(floor(n/2)) + s_16(floor(n/4)) + s_16(floor(n/8))`.

`test_reduce.py` checks the theorem exhaustively for all `n < 2^16` and the pairs
`(a,b)` in `{(1,2),(1,3),(1,4),(2,2),(2,3),(3,2),(1,5),(2,4)}`, re-deriving
`s_b` from plain base conversion, plus the underlying subsample identity and the
two named special cases.

## Why this matters for the cross-base prime bias

`RESULTS.md` reports a puzzle: among primes up to `10^10`, the standardized
digit-sum parity bias is very large in base 2 (`z ≈ +134`) yet essentially
vanishes in base 8 (`|z| < 6`, sign-oscillating), even though both are
power-of-two digit sums. The reduction makes these two numbers *the same object
viewed two ways*. For base 2 and base 8 it says exactly

```
s_2(n) = s_8(n) XOR s_8(floor(n/2)) XOR s_8(floor(n/4)).
```

So the base-2 bias over the primes is the bias of an XOR of three shifted base-8
parities. A large bias in an XOR-sum does **not** require any single term to be
biased: it can come entirely from the joint correlations of the three terms along
the primes. Concretely, writing `e(m) = (-1)^{s_8(m)}`, the base-2 exponential
sum over primes factors as

```
sum_{p <= x} (-1)^{s_2(p)}  =  sum_{p <= x} e(p) e(floor(p/2)) e(floor(p/4)).
```

Each single factor summed over primes is near zero (that is the measured base-8
near-cancellation). Their product need not be. Therefore **a near-zero marginal
base-8 bias is fully compatible with a large base-2 bias** — the apparent
paradox is not one. It relocates the open question from "why is base 8 unbiased?"
to "why do three shifted base-8 parities correlate so strongly along the primes
while each is individually near-balanced?" — a joint-correlation question, not a
single-marginal one.

## What is new, and what is not

- The subsample identity (`s_{2^k}` parity as a thin subsample of binary bits) is
  **folklore** for automatic sequences, and is already credited as such in
  `SUBSAMPLE.md`.
- The inter-base reduction above is an **elementary consequence** of it — a
  residue-class split of the sampled bit positions. It is short. I did not find
  it stated in this form, nor applied to relate the cross-base prime digit-sum
  biases, in the prior art I searched (Gupta–Sridhar arXiv:1812.11841;
  Mauduit–Rivat, *Annals* 2010; OEIS A000069/A130911/A380192/A111309). I claim
  only that specific application and framing as new, not the underlying identity.
- The compatibility observation (near-zero base-8 bias does not force near-zero
  base-2 bias) is a **logical consequence** of the exact identity, not a new
  theorem about primes. It explains why the two measurements do not contradict
  each other; it does **not** predict or explain the sign or size of any base's
  bias. The base-8 mechanism — why the k = 3 subsample correlates so weakly with
  primality on its own — remains open.

## Reproduce

```
python -m unittest test_reduce -v
python reduce.py
```

Both are standard-library only.

— MV2, Munim, Inc.
