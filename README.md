# prime-digit-parity

Reproducible study of the **digit-sum parity bias of primes** across even bases —
and a specific, testable question about a published "unexplained" result.

Gupta & Sridhar (arXiv:1812.11841) reported that over the *first 50,000,000
primes* the base-10 digit sum is far more often odd than even, by about an order
of magnitude more than chance, with no known explanation. That measurement uses a
**prime-count cutoff**. This repo asks whether the effect is an artifact of that
choice versus a **value cutoff** (`primes ≤ x`), and proves the integer-level
facts that make the question sharp.

Four things are proved outright (see [`NOTES.md`](NOTES.md), checked in
[`test_core.py`](test_core.py)):

1. **Odd bases are trivial** — every prime `p > 2` has an odd digit sum in every
   odd base, so only even bases carry information.
2. **Even-base ranges are exactly balanced** — `[0, b^k)` splits exactly 50/50 by
   digit-sum parity for even `b`.
3. **The integer excess is bounded by 1 at *every* cutoff** — for even `b`, the
   odd−even digit-sum parity excess over `[0, X)` is always in `{-1, 0, +1}`
   (closed form in `core.integer_parity_excess_closed`). This *sharpens* Fact 2
   and corrects an overstatement in v0.1.0: the integer "boundary term" is tiny
   everywhere, so the count-vs-value choice **cannot** inflate a large bias
   through the ambient integers — any real bias is a property of the primes.
4. **The base-10 bias lives in the high digits** — every prime `p > 5` ends in an
   odd digit `{1,3,7,9}`, so `s_10(p) mod 2 = 1 XOR s_10(⌊p/10⌋) mod 2`. The
   entire base-10 parity bias equals (minus) that of the "beheaded" primes
   `⌊p/10⌋`.

The open piece — whether the normalized base-10 prime bias `z(x) = (odd−even)/√π(x)`
stays `O(1)` at value cutoffs `x = 10^k` (artefact) or grows one-signed (real) —
is computed by the job below. **Results are pending that run; no numeric claim
about primes is made here before it is computed.**

## Layout

| file | what |
|------|------|
| `core.py` | stdlib-only reference: digit sums, sieve, brute-force + closed-form excess |
| `fast.py` | numpy segmented sieve + vectorised parity + checkpoints |
| `job.py` | runs the study to `10^9`, writes `results.json` / `summary.json` |
| `test_core.py` | correctness + the four proved facts (no numpy) |
| `NOTES.md` | the research note, prior art, and the hypothesis |

## Run

```bash
python -m unittest -v           # correctness (no dependencies)
pip install -r requirements.txt
python job.py --out-dir out     # the full computation to 10^9
```

`job.py` records, for bases `{2,4,6,8,10,16}`, the odd−even excess and
`z = excess/√n` at value cutoffs `10^k` and at prime-count cutoffs every `10^6`
primes. Base 2 (the Shevelev/Chaffin odious-prime race, verified to `10^19`) is a
control that validates the pipeline's sign.

## Provenance and license

Work by **Munim, Inc.** This project was **produced with AI** (autonomous
research agents) and is published openly under that disclosure. Nothing here is
claimed to be written or verified by a human unless a human says so.

- Code: **Apache-2.0** (`LICENSE`).
- Text and data (`NOTES.md`, result tables): **CC BY 4.0**.

Corrections and prior-art pointers are welcome via issues.
