# prime-digit-parity

Reproducible study of the **digit-sum parity bias of primes** across even bases,
starting from a published "unexplained" result.

Gupta & Sridhar (arXiv:1812.11841, *A Peculiarity in the Parity of Primes*)
reported that over the *first 50,000,000 primes* the base-10 digit sum is more often
odd than even, by about an order of magnitude more standard deviations than chance
(fair coin), with no known explanation.

**Where this stands (v0.1.2; see [`RESOLUTION.md`](RESOLUTION.md)).**

- *Settled by prior art:* by the Mauduit–Rivat theorem (*Annals of Mathematics* 171,
  2010) the digit-sum parity of primes has density 1/2 in base 10, so there is no
  limiting bias; any bias is a finite-sample, lower-order effect.
- *Not settled by prior art:* the theorem only upper-bounds the error term. It does
  **not** explain a few-sigma bias at 5·10^7 primes, and an earlier draft of this
  repo that said it did was wrong (corrected in v0.1.2). Whether the base-10 `z`
  stays O(1) or keeps one sign and grows is the open question.
- *Being measured:* a sandbox job computes the excess for bases `{2,4,6,8,10,16}` at
  value cutoffs `10^k`, at prime-count cutoffs (including `N = 5·10^7`), and along a
  per-segment trajectory, to `10^9` and then `10^10`. Readings are pre-registered in
  `RESOLUTION.md`. **No numeric claim about primes is made before it returns.**

Six things are proved outright (see [`NOTES.md`](NOTES.md), checked in
[`test_core.py`](test_core.py)):

1. **Odd bases are trivial** — every prime `p > 2` has an odd digit sum in every
   odd base, so only even bases carry information.
2. **Even-base ranges are exactly balanced** — `[0, b^k)` splits exactly 50/50 by
   digit-sum parity for even `b`.
3. **The integer excess is bounded by 1 at *every* cutoff** — for even `b`, the
   odd−even digit-sum parity excess over `[0, X)` is always in `{-1, 0, +1}`
   (closed form in `core.integer_parity_excess_closed`). So the count-vs-value
   choice **cannot** inflate a large bias through the ambient integers.
4. **The base-10 bias lives in the high digits** — every prime `p > 5` ends in an
   odd digit `{1,3,7,9}`, so `s_10(p) mod 2 = 1 XOR s_10(⌊p/10⌋) mod 2`.
5. **Casting out (b−1)'s pins the degeneracy** — `s_b(n) ≡ n (mod b−1)`, so digit-sum
   parity is degenerate exactly when `b` is odd. Base 10 is non-degenerate.
6. **A bounded control (new)** — the integers coprime to 10 (a superset of the
   primes `> 5`) have base-10 digit-sum parity excess `|E| ≤ 8` at *every* cutoff
   (`core.coprime10_parity_excess_closed`). Any `√n`-scale excess of the primes is
   therefore a property of the primes, not of the balanced ambient set.

## Layout

| file | what |
|------|------|
| `core.py` | stdlib-only reference: digit sums, sieve, brute-force + closed-form excess, controls |
| `fast.py` | numpy segmented sieve + vectorised parity + checkpoints |
| `job.py` | runs the study, writes `results.json` / `summary.json` / `summary.txt` |
| `test_core.py` | correctness + the six proved facts (no numpy) |
| `NOTES.md` | the research note, prior art, and the six facts |
| `RESOLUTION.md` | what Mauduit–Rivat does and does not settle; pre-registered readings |

## Run

```bash
python -m unittest -v                          # correctness (no dependencies)
pip install -r requirements.txt
python job.py --out-dir out --limit 1000000000 # the computation
```

Base 2 (the Shevelev/Chaffin odious-prime race, verified to `10^19`) is a control
that validates the pipeline's sign.

## Provenance and license

Work by **Munim, Inc.** This project was **produced with AI** (autonomous
research agents) and is published openly under that disclosure. Nothing here is
claimed to be written or verified by a human unless a human says so.

- Code: **Apache-2.0** (`LICENSE`).
- Text and data (`NOTES.md`, `RESOLUTION.md`, result tables): **CC BY 4.0**.

Corrections and prior-art pointers are welcome via issues.
