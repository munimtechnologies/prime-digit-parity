"""gap_corr_job.py -- lag-1 serial correlation of prime-gap digit-sum parity.

The sealed-sandbox measurement script behind GAP-CORR.md (prime-digit-parity
v0.7.0). Follow-up to v0.6.2 (GAP-PARITY.md): v0.6.2 showed the TOTAL
digit-sum-parity excess E_b^gap of the prime gaps g_n = p_{n+1}-p_n is large in
every base -- but that total is fully fixed by the gap-size histogram. This job
asks the open question it leaves: are consecutive gap parities chi_b(g_n),
chi_b(g_{n+1}) serially correlated BEYOND the fixed parity multiset? (the Lemke
Oliver-Soundararajan analogue for gap parity).

Observable, per base b in {2,4,6,8,10,16} and cutoff in {1e6,1e7,1e8}:
    S_b = sum_{n=1}^{N-1} chi_b(g_n) chi_b(g_{n+1})      (same-parity - different),
with chi_b(m) = (-1)^{s_b(m)}, N = n_gaps. The marginal-preserving (random
permutation) null fixes the multiset of chi values; its exact mean is
    E_perm[S_b] = (E_b^2 - N)/N        (hand-derived, see notes/prereg),
and its exact variance is Var[S_b] = (N-1) - 2 E_b^2/N + E_b^4/(N^2(N-1))
    - 1/(N-1)  (NULL-MOMENTS.md / null_moments.perm_var, v0.6.3).
We confirm the mean empirically with K fixed-seed shuffles and report the
connected signal z_b = (S_b - mean_shuffle)/std_shuffle.

HARD GATE (job exits non-zero unless all pass -- a failure is a bug, not a result):
  * pi(x) = 78498 / 664579 / 5761455 exactly;
  * telescoping sum of gaps == p_M - 2 at each cutoff;
  * exactly ONE odd gap (g_1 = 1);
  * half-gap Thue-Morse identity chi_2(g) == chi_2(g/2) for even gaps;
  * two-routine digit-sum guard for the first 200000 gaps, every base;
  * REPRODUCE v0.6.2 exactly: E_b at each cutoff equals the published values;
  * NULL-MEAN identity: empirical shuffle mean within 5*std/sqrt(K) of (E_b^2-N)/N
    at 1e6 and 1e7.

Hypotheses (prereg, decided before numbers): C0 null (|z_b|=O(1)); C1 real lag-1
structure in some base; C2 only base 8 special.

Work by Munim, Inc., produced with AI. MV2 of Munim, Inc.
Part of munimtechnologies/prime-digit-parity. Code Apache-2.0; data CC BY 4.0.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

BASES = (2, 4, 6, 8, 10, 16)
CUTOFFS = (1_000_000, 10_000_000, 100_000_000)
PI_ANCHORS = {1_000_000: 78498, 10_000_000: 664579, 100_000_000: 5761455}
GUARD_GAPS = 200_000
# v0.6.2 published excesses E_b = #odd - #even of s_b(g_n); reproduce EXACTLY.
EXCESS_ANCHORS = {
    1_000_000: {2: -7563, 4: 16607, 6: 3031, 8: -11497, 10: -15033, 16: -38231},
    10_000_000: {2: -77744, 4: 132234, 6: 31814, 8: -69764, 10: -96440, 16: -264548},
    100_000_000: {2: -756388, 4: 1051924, 6: 341536, 8: -437046, 10: -665788, 16: -1901898},
}
SHUFFLES = {1_000_000: 2000, 10_000_000: 800, 100_000_000: 200}
SEED = 20261016


def sieve_primes(limit):
    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i::i] = b"\x00" * len(range(i * i, limit + 1, i))
        i += 1
    return s


def parity_div(m, b):
    """Parity of the base-b digit sum of m, by repeated division."""
    s = 0
    while m:
        s ^= (m % b) & 1
        m //= b
    return s


def parity_str(m, b):
    """Independent base-b digit-sum parity via explicit digit list."""
    if m == 0:
        return 0
    digits = []
    while m:
        digits.append(m % b)
        m //= b
    return sum(digits) & 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--limit", type=int, default=100_000_000)
    args = ap.parse_args()
    limit = args.limit

    print(f"sieving primes up to {limit} ...", flush=True)
    sieve = sieve_primes(limit)

    # Per-base parity bit (0 = even digit sum/chi=+1, 1 = odd/chi=-1) for each gap.
    par_seq = {b: bytearray() for b in BASES}
    excess = {b: 0 for b in BASES}
    halfgap_chi2 = 0
    gap_sum = 0
    odd_gaps = 0
    n_gaps = 0
    prime_count = 0
    prev = None
    cur_prime = None
    guard_ok = True
    snapshots = {}
    snap_order = sorted(CUTOFFS)
    cut_iter = iter(snap_order)
    next_cut = next(cut_iter, None)

    for p in range(2, limit + 1):
        if not sieve[p]:
            continue
        prime_count += 1
        cur_prime = p
        if prev is not None:
            g = p - prev
            n_gaps += 1
            gap_sum += g
            if g & 1:
                odd_gaps += 1
            for b in BASES:
                par = parity_div(g, b)
                par_seq[b].append(par)
                excess[b] += 1 if par else -1
                if n_gaps <= GUARD_GAPS and par != parity_str(g, b):
                    guard_ok = False
            h = g >> 1 if (g & 1) == 0 else g
            halfgap_chi2 += 1 if parity_div(h, 2) else -1
        prev = p

        if next_cut is not None and prime_count == PI_ANCHORS[next_cut]:
            snapshots[next_cut] = {
                "cutoff": next_cut,
                "p_M": cur_prime,
                "prime_count": prime_count,
                "n_gaps": n_gaps,
                "gap_sum": gap_sum,
                "odd_gaps": odd_gaps,
                "excess": dict(excess),
                "halfgap_chi2_excess": halfgap_chi2,
            }
            next_cut = next(cut_iter, None)

    # ---- method gate (identical to v0.6.2) ----
    errors = []
    if not guard_ok:
        errors.append("two-routine digit-sum guard failed")
    for cut in snap_order:
        if cut > limit:
            continue
        if cut not in snapshots:
            errors.append(f"no snapshot at pi({cut})")
            continue
        snap = snapshots[cut]
        if snap["prime_count"] != PI_ANCHORS[cut]:
            errors.append(f"pi({cut})={snap['prime_count']} != {PI_ANCHORS[cut]}")
        if snap["gap_sum"] != snap["p_M"] - 2:
            errors.append(f"telescoping at {cut}: {snap['gap_sum']} != {snap['p_M']-2}")
        if snap["odd_gaps"] != 1:
            errors.append(f"odd gaps at {cut} = {snap['odd_gaps']} != 1")
        if snap["halfgap_chi2_excess"] != snap["excess"][2]:
            errors.append(f"half-gap chi2 identity at {cut}")
        for b in BASES:
            want = EXCESS_ANCHORS[cut][b]
            if snap["excess"][b] != want:
                errors.append(f"excess reproduce at {cut} base {b}: "
                              f"{snap['excess'][b]} != v0.6.2 {want}")
    if errors:
        print("GATE FAILED:", flush=True)
        for e in errors:
            print("  " + e, flush=True)
        sys.exit(1)
    print("method gate PASSED (pi, telescoping, single odd gap, half-gap identity, "
          "two-routine guard, v0.6.2 excess reproduced).", flush=True)

    # ---- lag-1 correlation + marginal-preserving shuffle null ----
    results = {}
    null_mean_errors = []
    for cut in snap_order:
        if cut not in snapshots:
            continue
        snap = snapshots[cut]
        N = snap["n_gaps"]
        K = SHUFFLES[cut]
        rng = np.random.default_rng(SEED + cut)
        per_base = {}
        for b in BASES:
            par = np.frombuffer(bytes(par_seq[b][:N]), dtype=np.uint8).astype(np.int8)
            chi = (1 - 2 * par).astype(np.int8)  # 0->+1, 1->-1
            # sanity: sum(chi) == -excess (chi=+1 for par 0 = even digit sum)
            assert int(chi.sum()) == -snap["excess"][b], (cut, b)
            s_actual = int(np.dot(chi[:-1].astype(np.int64), chi[1:].astype(np.int64)))
            e_perm = (snap["excess"][b] ** 2 - N) / N  # closed-form permutation mean
            # empirical shuffle null (fixed seed)
            samples = np.empty(K, dtype=np.float64)
            work = chi.astype(np.int64)
            for k in range(K):
                rng.shuffle(work)
                samples[k] = np.dot(work[:-1], work[1:])
            mean_sh = float(samples.mean())
            std_sh = float(samples.std(ddof=1))
            z = (s_actual - mean_sh) / std_sh if std_sh > 0 else 0.0
            per_base[b] = {
                "excess": snap["excess"][b],
                "S_actual": s_actual,
                "E_perm_mean": e_perm,
                "shuffle_mean": mean_sh,
                "shuffle_std": std_sh,
                "K": K,
                "z": z,
            }
            # null-mean identity gate at the two affordable cutoffs
            if cut in (1_000_000, 10_000_000):
                tol = 5 * std_sh / math.sqrt(K)
                if abs(mean_sh - e_perm) > max(tol, 1.0):
                    null_mean_errors.append(
                        f"null-mean mismatch cut {cut} base {b}: emp {mean_sh:.2f} "
                        f"vs closed form {e_perm:.2f} (tol {tol:.2f})")
        results[str(cut)] = {"cutoff": cut, "n_gaps": N, "bases": per_base}

    if null_mean_errors:
        print("GATE FAILED (null-mean identity):", flush=True)
        for e in null_mean_errors:
            print("  " + e, flush=True)
        sys.exit(1)
    print("null-mean gate PASSED (empirical shuffle mean == (E_b^2-N)/N).", flush=True)

    os.makedirs(args.out_dir, exist_ok=True)
    out = {
        "limit": limit,
        "bases": list(BASES),
        "seed": SEED,
        "results": results,
        "note": ("S_b = sum chi_b(g_n) chi_b(g_{n+1}) over prime gaps; "
                 "z = (S_b - shuffle_mean)/shuffle_std, marginal-preserving null. "
                 "E_perm_mean = (E_b^2 - N)/N. See notes/gap-parity-correlation-prereg.md."),
    }
    with open(os.path.join(args.out_dir, "gap_parity_corr.json"), "w") as f:
        json.dump(out, f, indent=2)

    print("\nlag-1 gap-parity correlation (z = connected signal in shuffle-std units):")
    for cut in snap_order:
        if str(cut) not in results:
            continue
        r = results[str(cut)]
        print(f"  cutoff {cut:>12}  n_gaps={r['n_gaps']:>9}")
        for b in BASES:
            d = r["bases"][b]
            print(f"      base {b:>2}  S={d['S_actual']:>10}  "
                  f"E_perm={d['E_perm_mean']:>12.1f}  z={d['z']:+8.3f}  K={d['K']}")
    print("\nwrote", os.path.join(args.out_dir, "gap_parity_corr.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
