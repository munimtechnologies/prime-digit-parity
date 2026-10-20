# Why the base-2 admit-even discriminator is not tautological

*Methodology note, MV2 of Munim, Inc. Produced with AI. Written run 121
(2026-10-20), before the `admit-even-lag-signature-v1` compute job has run. It
does **not** alter the locked predictions in
`notes/admit-even-lag-signature-prereg.md`; it records an interpretive guard so
the eventual grade is not over-read in either direction.*

*Refined run 127 (2026-10-20): the load-bearing gap bound below is now exact and
independently cited (OEIS A045882), replacing the earlier empirical "about 10".*

*Corrected run 131 (2026-10-20): the Filaseta-Trifonov sublinear-gap citation
below was mis-attributed to "J. Number Theory 1992". The `x^{1/5} log x`
squarefree-gap bound is from their 1990 Allerton Park paper (Progr. Math. 85,
Birkhauser, pp. 235-253), checked against the Birkhauser record; their 1992
"II" paper is in the J. London Math. Soc., and the J. Number Theory papers are
Filaseta's solo short-interval works (1988, 1990). Corrected here.*

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

Squarefree gaps are **small in the computed range**, and the bound is now exact.
The maximal gap between consecutive squarefree integers below 1e8 is **exactly
10**: a gap of 10 — equivalently, a run of 9 consecutive non-squarefree integers —
first occurs at 8,870,024, whereas the next larger gap, 11 (a run of 10
consecutive non-squarefree integers), first occurs only at 221,167,422 > 1e8.
This is independently confirmed by OEIS **A045882** ("smallest term of first run
of at least n consecutive integers which are not squarefree"): a(9) = 8,870,024
< 1e8 and a(10) = 221,167,422 > 1e8, so no squarefree gap below 1e8 exceeds 10.
Rigorously, gaps between consecutive squarefree numbers are known to be
sublinear: Filaseta & Trifonov proved that for some `c > 0` the interval
`(x, x + c*x^{1/5} log x]` contains a squarefree integer for all large `x`
(M. Filaseta and O. Trifonov, "On gaps between squarefree numbers," in
*Analytic Number Theory* (Allerton Park, IL, 1989), Progr. Math. 85, Birkhauser,
Boston, 1990, pp. 235-253). The exact value `= 10` is specific to the window
`[1, 1e8]`, not a claim for all `x`.

## Large bases (b >= 11): the objection holds — tautology

If every gap `g <= 10 < b`, then `g` is a single base-`b` digit, so its digit
sum is `g` itself and

    chi_b(g) = g mod 2   (the parity of the gap length).

For an odd-only set every gap is even, hence `chi_b(g) = 0` for **all** gaps:
the character is constant and necessarily shows no serial structure. At any base
`b >= 11` the "signature absent" outcome for odd-only sets is therefore
**tautological** and carries no mechanistic content. Any admit-even experiment
run at a large base would be uninformative on the odd-only side. (The threshold
is exactly `b >= 11` precisely because the below-1e8 gap maximum is exactly 10:
`b = 11` is the smallest base in which every occurring squarefree gap is a
single digit.)

## Base 2: the objection fails — a genuine structural claim

The registered experiment runs on the **base-2** stream only
(`x_n = P_2(g_n) = t(g_n)`). Here the collapse does **not** happen, because the
Thue-Morse parity of an *even* number is not constant. Write an even gap as
`g = 2h`. Appending a zero bit does not change the number of 1 bits, so

    popcount(2h) = popcount(h)   =>   t(2h) = t(h).            (identity A)

Thus for an odd-only set (all gaps even, `g_j = 2 h_j`),

    chi_2(g_j) = t(2 h_j) = t(h_j),                            (reduction B)

the Thue-Morse parity of the **halved** gaps `h_j = g_j / 2`. This sequence is
*not* constant: among the even gaps that actually occur below 1e8
(`g in {2,4,6,8,10}`; see the exact bound above), the Thue-Morse values are
`t(g) = 1,1,0,1,0`, i.e. both values appear. (The next even gap, `g = 12`, first
occurs only at 221,167,422 > 1e8; including it changes nothing — `t(12) = 0` — so
the live two-valued signal is already present within the computed range.) So the
base-2 character remains a live, two-valued signal even when every gap is even.

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
3. In base 2, over the small even gaps `{2,4,6,8,10,12}` (the even gaps up to and
   just past the exact below-1e8 maximum of 10), `chi_2` takes **both** values
   `{0, 1}` (base 2 does not collapse).

The exact below-1e8 gap maximum of 10 is a lookup against OEIS A045882
(a(9) = 8,870,024; a(10) = 221,167,422), not a claim produced by this repo's own
data alone.

---

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
