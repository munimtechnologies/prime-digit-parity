# Why the base-2 admit-even discriminator is not tautological

*Methodology note, MV2 of Munim, Inc. Produced with AI. Written run 121
(2026-10-20), before the `admit-even-lag-signature-v1` compute job has run. It
does **not** alter the locked predictions in
`notes/admit-even-lag-signature-prereg.md`; it records an interpretive guard so
the eventual grade is not over-read in either direction.*

## The objection this note answers

The admit-even experiment contrasts B-free sets that **admit even integers**
(squarefree coprime to 3 / 5 / 15 — these keep even members, so odd gaps occur)
against sets that **forbid even integers** (odd squarefree, squarefree coprime
to 6 — every member is odd, so every consecutive gap is even). The
pre-registered reading is: if the lag-2 "unique-quiet" signature is PRESENT for
the even-admitting sets and ABSENT for the odd-only sets, the driver is
*admitting even integers* (the parity of `n`), not the full `mu^2` structure.

A fair referee will object: **isn't that circular?** If you remove every even
integer, all gaps are even, so the gap-parity character `chi` is constant, and
of course a constant signal carries no serial correlation. If so, "signature
absent for odd-only sets" would be an artifact of `chi` degenerating, not
evidence about any mechanism.

This note shows the objection is **correct for large bases but false for the
base actually registered (base 2)**, and pins down exactly why.

## The character and the gap-size fact

For a sieved set with consecutive members `..., n_j, n_{j+1}, ...`, the gap is
`g_j = n_{j+1} - n_j`. The character studied across this arc is

    chi_b(g) = (digit sum of g in base b)  mod 2.

In base 2 this is the Thue-Morse value `t(g) = popcount(g) mod 2` (A010060).

Squarefree gaps are **small in the computed range**: every gap between
consecutive squarefree numbers below 1e8 is at most about 10 (this repo's own
committed squarefree gap data; summarized in the dev.to write-up). Rigorously,
gaps between squarefree numbers are known to be sublinear (Filaseta-Trifonov,
*J. Number Theory* 1992, gap `<< x^{1/5} log x`); the concrete bound `<= 10` is
empirical **within** `[1, 1e8]`, not a claim for all `x`.

## Large bases (b >= 11): the objection holds — tautology

If every gap `g <= 10 < b`, then `g` is a single base-`b` digit, so its digit
sum is `g` itself and

    chi_b(g) = g mod 2   (the parity of the gap length).

For an odd-only set every gap is even, hence `chi_b(g) = 0` for **all** gaps:
the character is constant and necessarily shows no serial structure. At any base
`b >= 11` the "signature absent" outcome for odd-only sets is therefore
**tautological** and carries no mechanistic content. Any admit-even experiment
run at a large base would be uninformative on the odd-only side.

## Base 2: the objection fails — a genuine structural claim

The registered experiment runs on the **base-2** stream only
(`x_n = P_2(g_n) = t(g_n)`). Here the collapse does **not** happen, because the
Thue-Morse parity of an *even* number is not constant. Write an even gap as
`g = 2h`. Appending a zero bit does not change the number of 1 bits, so

    popcount(2h) = popcount(h)   =>   t(2h) = t(h).            (identity A)

Thus for an odd-only set (all gaps even, `g_j = 2 h_j`),

    chi_2(g_j) = t(2 h_j) = t(h_j),                            (reduction B)

the Thue-Morse parity of the **halved** gaps `h_j = g_j / 2`. This sequence is
*not* constant: among the even gaps that actually occur, `g = 2,4,6,8,10,12`
give `t(g) = 1,1,0,1,0,0`, i.e. both values appear. So the base-2 character
remains a live, two-valued signal even when every gap is even.

Consequently the pre-registered **P-AE-odd = ABSENT** prediction is a real,
falsifiable structural claim: *the Thue-Morse parity of the halved gaps of an
odd squarefree set has no lag-2 unique-quiet signature.* It is not an artifact
of `chi` going constant. Symmetrically, **P-AE-c3 = PRESENT** is a claim that
admitting even integers — which reintroduces odd gaps and the half-shifted
base-4 stride phase `P_2(g) = P_4(g) XOR P_4(floor(g/2))` — restores the
signature. The experiment can genuinely come out either way.

## What this changes

- It bounds where the phenomenon can even be asked about: for squarefree gaps in
  the computed range the question is only non-trivial at **small bases (2 and
  4)**. Any future "admit-even" style probe must stay at base 2/4; a large-base
  version would be vacuous.
- It tells us how to read the forthcoming grade honestly:
  - If P-AE-odd is ABSENT in base 2, that is informative (halved-gap Thue-Morse
    parity lacks the signature), **not** a degeneracy artifact.
  - If, against the prereg, the odd-only set *did* show the signature, the
    "even-admission" story would be wrong and the driver would sit in the
    halved-gap structure itself — equally publishable.

## Verification

`test_base2_nontautology.py` (pure standard library, runs in milliseconds)
proves the exact facts this note relies on:

1. `t(2h) = t(h)` for all `h` in `[0, 200000]` (identity A) — hence
   `chi_2(2h) = t(h)` (reduction B).
2. For `b >= 11` and `1 <= g < b`, `digit_sum(g, b) = g` and
   `chi_b(g) = g mod 2`; over all even `g < b` the character is constant `{0}`
   (the large-base tautology).
3. In base 2, over the even gaps `{2,4,6,8,10,12}` that actually occur for
   squarefree sets, `chi_2` takes **both** values `{0, 1}` (base 2 does not
   collapse).

---

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
