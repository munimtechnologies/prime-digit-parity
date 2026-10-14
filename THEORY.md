# Theory: digit-sum parity of primes in even bases

*Part of `munimtechnologies/prime-digit-parity`. Work by Munim, Inc.,
produced with AI (autonomous research agents). Text CC BY 4.0, code
Apache-2.0. MV2 of Munim, Inc.*

This note separates what is **proved** (elementary, machine-checked in
`classes.py` / `test_classes.py`) from what is **conjectural and awaiting
computation** (the size and base-ordering of the prime excess). No numeric
claim about the prime excess is made here; those verdicts wait on the
compute sandbox.

## Notation

For base `b >= 2` let `s_b(n)` be the base-`b` digit sum and
`f_b(n) = (-1)^{s_b(n)}` (+1 if the digit sum is even, -1 if odd). For a
set `S` of non-negative integers and cutoff `X`, the **parity excess** is

```
E_b(S; X) = sum_{n in S, n < X} f_b(n)
          = #{n in S, n < X : s_b(n) even} - #{n in S, n < X : s_b(n) odd}.
```

Gupta & Sridhar (*A Peculiarity in the Parity of Primes*, arXiv:1812.11841)
observed that for `S = primes`, `b = 10`, the digit sum is more often odd
than even over the first 5*10^7 primes (excess `-64768` in the odd-even
convention; `|z| ~ 9.16` against a fair-coin baseline). This note is about
*why even bases carry any signal at all*, and about the exact behaviour of
two control sets (all integers; integers coprime to `b`) that the primes
sit inside.

## 1. Only even bases carry a signal (recap of Facts 1 & 5)

Since `s_b(n) ≡ n (mod b-1)`, and more sharply, parity: for **odd** `b`
every power `b^k` is odd, so `n = sum d_k b^k ≡ sum d_k = s_b(n) (mod 2)`.
Thus `f_b(n) = (-1)^n` and the digit-sum parity is degenerate — it is just
the parity of `n`. Every prime `p > 2` is odd, so `s_b(p)` is odd in every
odd base: no anomaly is possible. For **even** `b`, `b^k` is even for
`k >= 1`, so `f_b` is a genuine Thue-Morse-type sequence. Base 10 is even,
hence non-degenerate. All subsequent statements are for even `b`.

## 2. The integer control is bounded by 1 (recap of Fact 3)

Splitting `[0, X)` into blocks `[b*j, b*(j+1))`: each full block has last
digit running `0..b-1`, and because `b` is even `sum_{r=0}^{b-1}(-1)^r = 0`,
so every completed block contributes 0. With `X = b*Q + R` only the partial
block survives, giving the **exact closed form**

```
E_b(Z; X) = D_b(X) = (-1)^{s_b(X//b)} * ((X % b) % 2)   in {-1, 0, +1}.
```

So the ambient integers are balanced to within 1 at *every* cutoff — the
count-cutoff-vs-value-cutoff sampling choice cannot inflate a large bias
through the integers. (`integer_excess_closed` in `classes.py`.)

## 3. New: the coprime control has the sharp bound phi(b)

The primes `> b` all lie in the reduced residue classes mod `b`. The
natural control set is therefore `S = {n : gcd(n, b) = 1}`. Writing a digit
`a` in `[0, b)` and `n = b*m + a` (so `s_b(n) = s_b(m) + a`, no carry):

```
E_b({n ≡ a mod b}; X) = (-1)^a * D_b(M_a),   M_a = #{m : b*m + a < X},
```

a value in `{-1, 0, +1}` for each single class. Summing over the `phi(b)`
reduced residues gives the

> **Theorem (sharp coprime bound).** For every even base `b` and every
> cutoff `X`,
> `|E_b({gcd(n,b)=1}; X)| <= phi(b)`,
> and the constant `phi(b)` is **attained** (already at `X = b`, where the
> coprime residues `1, 3, ... ` each contribute `-1`).

Two consequences:

- **It sharpens Fact 6.** The earlier note bounded the base-10 coprime
  excess by `|E| <= 8`; the tight constant is `phi(10) = 4` (e.g.
  `E_10(coprime; 10) = -4`, and `E_10(coprime; 110) = +4`). The factor-2
  gap was slack in the old argument, not a real phenomenon.
- **Every reduced residue mod an even base is odd.** If `a` were even then
  `gcd(a, b) >= 2`. Hence every prime `p > b` ends in an odd base-`b`
  digit, so `s_b(p) ≡ s_b(floor(p/b)) + 1 (mod 2)`: the prime digit-sum
  parity in **any** even base equals the opposite of the beheaded number's
  parity. This generalises Fact 4 (base 10: last digit in `{1,3,7,9}`) to
  all even bases.

All of Sections 1-3 are proved and checked by `test_classes.py`
(closed form vs brute force for bases `{2,4,6,8,10,12,16}`, the `phi(b)`
bound, and its attainment).

**Takeaway for the open problem.** Both controls the primes live inside —
the integers (bounded by 1) and the coprime integers (bounded by `phi(b)`,
a *constant*) — are balanced to within a constant at every cutoff. So the
`~9.16 sigma`, growing base-10 prime excess is **not** inherited from its
ambient sets; it is a property of the primes themselves. The question is
genuinely about primes, not about the sampling window or the residue
structure.

## 4. Conjectural: the size and base-ordering of the prime excess

For the primes, Mauduit & Rivat (*Sur un probleme de Gelfond: la somme des
chiffres des nombres premiers*, Annals of Math. 171 (2010), 1591-1646)
proved equidistribution of `s_b(p) mod m` with a **power-saving** error:
for `(m, b-1) = 1`,

```
E_b(primes; x) = O(x^{1 - sigma_{b,m}} ),   sigma_{b,m} > 0.
```

This is an **upper bound only**: it forbids a limiting bias (density 1/2)
but does not pin the finite-sample size, and it does not prove the excess
grows. The saving `sigma_{b,m}` is governed by a minor-arc constant that
**decreases like `1/log b`** (made explicit for the `m=2` surjectivity
threshold in arXiv:2606.04677, which gives `~0.0015` at `b=10`). Standardise
by the fair-coin scale `sqrt(pi(x))`:

```
z_b(x) = E_b(primes; x) / sqrt(pi(x)).
```

The ceiling then reads `|z_b(x)| = O(x^{1/2 - sigma_{b,2}} (log x)^{1/2})`.
Because `sigma_{b,2}` shrinks with `b`, larger even bases have a *larger*
admissible growth exponent. This motivates three **predictions** (not yet
tested — the verdict waits on the trajectory computation):

- **(P1)** Fitting `alpha_b` as the slope of `log|z_b(10^k)|` vs `k log 10`,
  the `alpha_b` are ordered increasingly in `b` across `{2,4,6,8,10,12,16}`.
- **(P2)** The raw excess `E_b(primes; x)` is one-signed and growing, not
  sign-changing, over the measured range (consistent with `z_10 ~ 9.16`).
- **(P3)** `alpha_2` is the smallest exponent (base 2 = the classical
  odious/evil race, verified balanced to `10^19`).

If (P1) holds, the headline is that *the Gupta-Sridhar anomaly's size
tracks the Mauduit-Rivat exponent* — the "unexplained" base-10 bias is the
`b=10` instance of a monotone cross-base law. If the exponents are flat or
the excess sign-changes, that is an honest negative (a single-cutoff
fluctuation rather than a growth law), and equally publishable. Either way,
no such claim is made until the sandbox returns the `z_b(10^k)` table.

## Status

- Sections 1-3: **proved**, machine-checked, independent of any computation.
- Section 4: **conjectural**, upper bound cited, predictions pre-registered.
  The measuring job (`prime-parity-trajectory-v1`) is queued; `RESULTS.md`
  will record the verdict on (P1)-(P3) once it returns.

## References

- S. Gupta, K. Sridhar. *A Peculiarity in the Parity of Primes.*
  arXiv:1812.11841.
- C. Mauduit, J. Rivat. *Sur un probleme de Gelfond: la somme des chiffres
  des nombres premiers.* Annals of Mathematics 171 (2010), 1591-1646.
- *Surjectivity / minor-arc threshold for digit sums* (explicit `1/log b`
  saving), arXiv:2606.04677.
