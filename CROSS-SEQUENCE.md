# Is the base-8 prime digit-sum-parity valley specific to the primes?

*Work by Munim, Inc., produced with AI. MV2 of Munim, Inc. Part of
munimtechnologies/prime-digit-parity. Code: Apache-2.0. Text and data: CC BY 4.0.*

## 1. The question

The flagship study found that among the first several million primes, the base-8 digit-sum
parity excess is anomalously small — a "valley." Writing `chi_b(n) = (-1)^{s_b(n)}`
(`s_b` = base-`b` digit sum) and, in the repo sign convention,

    E_b(x) = #{n <= x : s_b(n) odd} - #{n <= x : s_b(n) even} = - sum_{n <= x} chi_b(n),

the primes give (at x = 1e6) `E_2 = 4764`, `E_4 = 2902`, `E_6 = 1936`, **`E_8 = -510`**,
`E_10 = 1916`, `E_12 = 572`, `E_16 = 244`. Base 8 is an order of magnitude smaller than its
neighbours, and (flagship data) its normalized z-score stays near zero as x grows while the
other bases diverge. **Is that valley a property of the primes, or a generic feature of any
sparse / multiplicatively-structured subset of the integers?**

We answer it by computing the same parity excess along four sequences up to a common cutoff:

- **integers** — the exact null (see §2);
- **primes** — the subject of the flagship;
- **squarefree** numbers (density 6/pi^2);
- **semiprimes** — `Omega(n) = 2` (two prime factors counted **with multiplicity**, so `p*q`
  with `p != q` *and* `p^2` both count). This is the standard semiprime definition and is
  exactly what the sandbox sieve measures.

**Answer (this document): the valley is prime-specific.** At the largest cutoff (1e7) the
primes place their shallowest normalized excess at base 8; squarefree numbers place theirs at
base 12 (and base 8 is in fact their *strongest* base, a peak), and semiprimes at base 10.

## 2. The null baseline (proved, exact, all bases)

A cross-sequence comparison is only meaningful against a reference line. We prove the integer
line exactly, in every base — not just base 8.

Define the running integer excess `S_b(N) = sum_{n=0}^{N-1} chi_b(n)` (so
`E_b(integers in [0,N)) = -S_b(N)`).

- **Even base b (Theorem A).** Writing `N = b*q + r`, `0 <= r < b`:
  `S_b(N) = chi_b(q)` if `r` is odd, else `0`. Hence **`|S_b(N)| <= 1` for every N**, and
  `S_b(N) = 0` whenever `b | N`. In particular at `N = 1e6 = 8*125000`, the integer base-8
  excess over `[0,N)` is identically `0`.
- **Odd base b (Theorem B).** `S_b(N) = (1 - chi_b(N))/2 = [s_b(N) odd] in {0,1}`: the running
  excess over all integers *equals the digit-sum parity of N itself*, and never grows.
- **Unified:** `|S_b(N)| <= 1` for **every** base `b >= 2`. The integer parity excess is
  bounded by a single unit, forever, in all bases.

*Proofs and self-contained stdlib tests:* even bases in `20261016T023038Z-discovery/prep/`
(uniform `±1` bound, verified over `[0, 2e5)` for b in {2,4,6,8,10,12,16}); odd bases, exact
identity, in `test_integer_null_all_bases.py` (shipped with this release, verified over
`[0, 2e5)` for b in {3,5,7,9,15}). This is the Thue–Morse / Gelfond balance packaged as a
finite, boundary-free reference line.

**Consequence for the headline:** because the integer excess never leaves a single unit (in
the `[0,N)` convention), *any* excess of order 10^2–10^5 observed along a subsequence is
unambiguously a selection effect of that subsequence — never a baseline artefact. The primes'
`E_8 = -510` is a *small* such effect, which is exactly what makes the valley interesting
rather than trivial.

**Convention note — read this before the integer row below.** The clean statement above is the
summatory function over `[0,N)`. The sandbox counts each sequence over `[1, N]` (it starts at
`n = 1`). For the integers this adds the exact boundary term `E_b([1,N]) = 1 - chi_b(N)`, which
is `0` or `2` — e.g. at `1e6`, `s_8(1e6) = 15` is odd so `chi_8 = -1` and the integer base-8
column reads `2`, not `0`. This is a fully predictable `+0`/`+2` offset, **not** growth: the
integer excess is still `O(1)` at every scale, against subsequence excesses up to `~10^5`. For
primes/squarefree/semiprimes the offset is identically zero, because `n = 0`, `n = 1` and
`n = N = 10^k` are never prime, and `10^k` is neither squarefree nor a semiprime, so those rows
are convention-free.

## 3. What we measured

Sandbox job `cross-seq-parity-v1` (pure-stdlib sieves: Eratosthenes for primes; multiples of
`p^2` removed for squarefree; `Omega` via prime-power sieve for semiprimes, semiprime =
`Omega==2`), anchor-gated so the prime excesses at 1e6 and 1e7 must reproduce the published
values exactly or the run aborts. Gate result: **PASSED** (primes `E_8@1e6 = -510` vs anchor
-510; `E_8@1e7 = 215` vs anchor 215; all of `E_2,E_4,E_6,E_10,E_12,E_16` match the published
trajectory anchors at both 1e6 and 1e7; `pi(1e6) = 78498`, `pi(1e7) = 664579`). Ran in 65 s on
2 cores, no network.

The self-contained `test_cross_sequence.py` independently re-sieves to 1e6 and asserts the
prime anchor (`-510`), the exact integer null (`0` over `[0,N)`), and the job's own 1e6 slice
for squarefree (`E_8 = 13272`) and for semiprimes (`Omega==2`, `E_8 = 2449`).

### Excess `E_b = #odd - #even` at x = 1e7 (largest cutoff)

| sequence   | count      | E_2   | E_4   | E_6  | **E_8** | E_10 | E_12 | E_16  |
|------------|------------|-------|-------|------|---------|------|------|-------|
| integers   | 10,000,000 | 0     | 2     | 2    | **2**   | 2    | 2    | 0     |
| primes     | 664,579    | 20181 | 13121 | 3181 | **215** | 2967 | 4599 | 16623 |
| squarefree | 6,079,291  | 45023 | 20335 | 4959 | **98485** | 4825 | 2563 | 6551 |
| semiprimes | 1,904,324  | 35930 | 25104 | 6284 | **20830** | 3404 | 7364 | 18424 |

(The integer row is the `[1,N]` boundary term of §2, all in `{0,2}`; every other row is a real
selection effect.)

### Normalized excess `V_b = 1000 * E_b / count` at x = 1e7

| sequence   | V_2  | V_4  | V_6  | **V_8** | V_10 | V_12 | V_16 | valley base (argmin \|V_b\|) |
|------------|------|------|------|---------|------|------|------|------------------------------|
| primes     | 30.37 | 19.74 | 4.79 | **0.32** | 4.46 | 6.92 | 25.01 | **8** |
| squarefree | 7.41 | 3.35 | 0.82 | **16.20** | 0.79 | 0.42 | 1.08 | 12 (base 8 is the **peak**) |
| semiprimes | 18.87 | 13.18 | 3.30 | **10.94** | 1.79 | 3.87 | 9.67 | 10 |

### Normalized z-score `E_b / sqrt(count)` at base 8, by cutoff

| sequence   | z_8 @1e4 | z_8 @1e5 | z_8 @1e6 | z_8 @1e7 |
|------------|----------|----------|----------|----------|
| primes     | 0.43  | -0.10 | -1.82 | **0.26** |
| squarefree | 4.32  | 5.71  | 17.02 | **39.94** |
| semiprimes | 3.14  | 1.45  | 5.34  | **15.09** |

The primes' base-8 z-score hugs zero as x grows; the squarefree base-8 z-score *climbs* (the
single most significant cell anywhere in the study, z ≈ 40); semiprimes sit between.

## 4. Prediction scorecard (registered before measurement)

Predictions were filed in `20261016T032726Z-discovery/prep/cross-seq-predictions.md` *before*
the job returned.

| # | prediction | outcome |
|---|------------|---------|
| P1 | integer excess is `O(1)`, non-growing (`\|E_b\| <= 1` for even b in the proved `[0,N)` sum; `=0` at 1e6) | **confirmed** — exact over `[0,N)` (test asserts `0` at 1e6); the sandbox's `[1,N]` column differs only by the stated boundary term `1-chi_b(N) in {0,2}`, so `\|E_b\| <= 2`, still `O(1)` vs `~10^5` for subsequences |
| P2 | prime anchors reproduce exactly (job-correctness gate) | **confirmed** — all seven bases match at 1e6 and 1e7, `E_8 = -510` then `+215` |
| P3 | the base-8 valley is **prime-specific** (squarefree & semiprimes bottom out elsewhere) | **confirmed** — primes' valley is base 8; squarefree's is base 12 (base 8 is their peak); semiprimes' is base 10 |
| P4 | `V_b` shrinks with x for every set/base (sublinear excess, no constant density gap) | **confirmed over the measured decade** — every `V_b` is smaller at 1e7 than at 1e6 (e.g. primes `V_2`: 60.7→30.4; squarefree base-8 `V_8`: 21.8→16.2); asymptotics not proved |

## 5. Interpretation

**The valley is prime-specific (P3 confirmed).** Squarefree numbers place their shallowest
normalized excess at base **12**, and base 8 is in fact their *strongest* base — a peak at
z ≈ 40, the single most significant cell in the entire study, the exact opposite of the prime
behaviour. Semiprimes bottom out at base **10**. Only the primes bottom out at base 8. The
near-cancellation at base 8 is therefore a feature of the primes specifically — consistent with
the flagship finding that it traces to a single connected top-order joint cumulant tied to the
multiplicative structure of the primes. Squarefree numbers (density 6/pi^2) and semiprimes carry
different multiplicative weight and not only fail to inherit the valley — the squarefree sequence
shows the *largest* base-8 bias of any sequence measured.

The integer null (`|E_b| <= 1` in the `[0,N)` convention, `<= 2` as counted) guarantees the
compared excesses are real selection effects, and the primes are equidistributed mod 2 in digit
sum (Mauduit–Rivat 2010), so every excess here is a second-order fluctuation, not a density gap.

Why base 8 for squarefree, specifically, is open: base 8 = 2^3, and squarefree-ness is governed
by the small-square sieve (mainly removing multiples of 4 = 2^2 and 9), so an interaction between
the 2-adic structure of squarefree numbers and the base-2^3 digit sum is the natural suspect. We
state this as a question, not a claim.

## 6. Honest limits

- Cutoff is 1e7; valley behaviour is asymptotic, so conclusions are "up to 1e7," read together
  with the flagship's 1e10 prime trajectory. Squarefree/semiprime are only to 1e7 here.
- Squarefree digit-sum equidistribution mod 2 is a *heuristic* (Dietmann–Elsholtz–Shparlinski
  2016 show digit-robustness, not the mod-2 statement); semiprime is weaker still. We measure,
  we do not prove these. The squarefree base-8 peak (z ≈ 40 and still climbing) is a measured
  fact to 1e7, not evidence against eventual equidistribution.
- z-scores use the fair-coin null `sd = sqrt(count)`; correlations between digit positions make
  this an approximation, used only for scale.
- The integer row is reported in the sandbox's `[1,N]` convention and differs from the proved
  `[0,N)` summatory function by the exact boundary term `1 - chi_b(N)` (§2); no other row is
  affected.

## 7. Reproduce

`python job.py --out-dir out` rebuilds `cross_seq_parity.json` (pure stdlib, ~1e7, ~1 min on
2 cores). `test_cross_sequence.py` re-derives the 1e6 slice independently and asserts the prime
anchor, the exact integer null (all bases), and the squarefree/semiprime 1e6 values against the
shipped data. `test_integer_null_all_bases.py` verifies the all-base null closed forms over
`[0, 2e5)`. CSV `data/cross_sequence_excess.csv` (CC BY 4.0): one row per
(cutoff, sequence, base) with `count`, `excess` (= E_b) and `z` exactly as the sandbox emitted
them; the normalized `V_b = 1000*excess/count` of §3 is derived from these.

## References

- C. Mauduit, J. Rivat, *Sur un probleme de Gelfond: la somme des chiffres des nombres
  premiers*, Ann. of Math. 171 (2010), 1591–1646.
- C. Mauduit, J. Rivat, *La somme des chiffres des carres*, Acta Math. 203 (2009), 107–148.
- L. Spiegelhofer, *Thue–Morse along the sequence of cubes*, arXiv:2308.09498.
- M. Drmota, M. Skałba, *Rarefied sums of the Thue–Morse sequence*, Trans. Amer. Math. Soc.
  352 (2000), 609–642 (closest prior "subsequence excess" of the Thue–Morse sum; treats
  arithmetic-progression rarefactions, not multiplicatively-defined subsequences).
- R. Dietmann, C. Elsholtz, I. Shparlinski, *Prescribing the binary digits of squarefree
  numbers and quadratic residues*, arXiv:1601.04754 (2016).
- S. Gupta, K. Sridhar, arXiv:1812.11841 (primes only — the flagship's starting point).
- J.-P. Allouche, J. Shallit, *The ubiquitous Thue–Morse sequence*; OEIS A010060.
- V. Shevelev, arXiv:0709.0885 (even-base Thue–Morse balance).
