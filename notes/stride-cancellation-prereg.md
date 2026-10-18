# Pre-registration: does base-2 parity signal come from base-4 stream cancellation?

*Pre-registered by MV2 of Munim, Inc. Produced with AI. Written before the
measurement is run. Filed, not yet queued (the shared compute FIFO holds two
earlier pre-registered jobs; this one is queued only after one of them is graded,
to avoid clogging the queue).*

## Background and the exact identity being exploited

From `STRIDE-DECOMP.md` (proved and exhaustively tested):

```
P_2(g) = P_4(g) XOR P_4(floor(g/2)).
```

Define, for a sieve's consecutive gaps g_1, g_2, ..., the two base-4 strided
streams
- a_n = P_4(g_n)            (even-position bits of the gap)
- b_n = P_4(floor(g_n/2))   (odd-position bits of the gap)

so that the base-2 parity stream is x_n = a_n XOR b_n (equivalently, in +/-1 form,
chi2_n = chi4a_n * chi4b_n).

The study's lag-1 serial statistic, in +/-1 form, is
`r1(stream) = mean_n s_n s_{n+1}` with the pre-registered permutation null
(mean 0, variance from the stored null-moment formulas in this repo), reported as
a z-score. We will compute r1 and its z-score for four streams: a (base-4 even),
b (base-4 odd), x = a XOR b (base-2), and the full base-4 observable itself.

## The question

The observed cross-sieve signature is: primes and squarefree numbers are base-2
loud; sums of two squares are base-2 quiet but base-4 loud. The decomposition
says base-2 is the XOR of the two base-4 streams. **Is the sums-of-two-squares
base-2 quiet because its two base-4 streams individually carry serial signal that
cancels under XOR, or because the streams carry no serial signal at all?**

## Pre-registered predictions (to be graded exactly as written)

Let z(.) denote the lag-1 z-score magnitude at a cutoff of 1e8 for each sieve.

- **P1 (cancellation for sums of two squares).** For sums of two squares,
  `max(z(a), z(b)) >= 30` while `z(x) = z(base-2) <= 30`. That is, at least one
  base-4 stream is strongly serially correlated but the base-2 XOR is weak. This
  is the "cancellation" outcome and would confirm the mechanism.
  - If instead `max(z(a), z(b)) < 30` as well, the mechanism is **falsified**: the
    base-2 quiet is because neither strided stream carries signal, not cancellation.

- **P2 (primes are not pure cancellation).** For primes, `z(x) = z(base-2) >= 100`
  (already established at ~120 at 1e8). We additionally predict that the base-2
  signal is not dominated by a single stream alone: `min(z(a), z(b)) >= 10`, i.e.
  both strided streams contribute. If one stream has `z < 10` while base-2 is
  loud, the base-2 signal is essentially carried by the other single stream.

- **P3 (squarefree mirrors primes).** For squarefree numbers, same shape as P2:
  `z(x) >= 100` and `min(z(a), z(b)) >= 10`.

## Method (fixed in advance)

1. Sieve each sequence up to 1e8 (primes by Eratosthenes; squarefree by a
   mu-sieve; sums of two squares by marking a^2+b^2). Form consecutive gaps.
2. For each gap compute a_n, b_n via the exact bit formulas in `stride_decomp.py`;
   form x_n = a_n XOR b_n and the full base-4 stream.
3. Compute lag-1 r1 for each stream and convert to a z-score with the repo's
   pre-registered permutation-null variance.
4. Grade P1-P3 against the thresholds above. Report all z-scores regardless of
   outcome, including negatives.

## Self-test baked into the job

The job re-derives, on a tiny slice, each sieve against brute factorisation, the
null variance against brute permutation averaging, and `x_n = a_n XOR b_n`
against a direct base-2 digit-sum. The job refuses to report headline numbers
unless every self-test passes.

## Honesty

Thresholds (30, 100, 10) are set before running and will not be moved. A null or
mixed result is a valid outcome and will be published as such. Nothing here has
been measured or refereed yet.
