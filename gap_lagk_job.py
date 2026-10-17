"""gap-parity-lagk-v1 -- lag-k serial correlation of prime-gap digit-sum parity.

Follow-on to prime-digit-parity v0.7.0 (GAP-CORR.md, the lag-1 result). v0.7.0
measured S_b^(1) = sum chi_b(g_n) chi_b(g_{n+1}) of the prime gaps g_n against a
marginal-preserving null and found consecutive gap parities ARE serially
correlated beyond their fixed multiset (bases 2/8/16 anti-correlated to -120 at
1e8; base 4 positively correlated). This job asks the open question it leaves:
does the gap-parity sequence carry serial memory BEYOND one step?

Observable, per base b in {2,4,6,8,10,16}, cutoff in {1e6,1e7,1e8}, lag k in 1..5:
    S_b^(k) = sum_{n=1}^{N-k} chi_b(g_n) chi_b(g_{n+k}),  chi_b(m) = (-1)^{s_b(m)}.
The marginal-preserving (random permutation) null fixes the multiset of chi values;
its exact mean for lag k is (hand-derived, see notes/gap-parity-lagk-prereg.md)
    E_perm[S_b^(k)] = (N - k) * (E_b^2 - N) / (N * (N - 1)),
which for k=1 reduces to (E_b^2 - N)/N, the lag-1 form shipped in v0.7.0. We
confirm it with K fixed-seed shuffles and report z_b^(k) = (S - mean)/std.

HARD GATE (job exits non-zero unless all pass -- a failure is a bug, not a result):
  * pi(x) = 78498 / 664579 / 5761455 exactly;
  * telescoping sum of gaps == p_M - 2 at each cutoff;
  * exactly ONE odd gap (g_1 = 1);
  * half-gap Thue-Morse identity chi_2(g) == chi_2(g/2) for even gaps;
  * two-routine digit-sum guard for the first 200000 gaps, every base;
  * REPRODUCE v0.6.2 excess E_b at each cutoff exactly;
  * REPRODUCE v0.7.0 lag-1 S_b^(1) exactly (regression; skipped only if the anchor
    dict below is left empty, with a printed warning);
  * NULL-MEAN identity: empirical shuffle mean of S_b^(k) within 5*std/sqrt(K) of
    (N-k)(E_b^2-N)/(N(N-1)) for every measured k at 1e6 and 1e7.

Hypotheses (prereg, decided before numbers): D0 null (|z|=O(1) all b,k); D1 decaying
memory; D2 periodic/base-specific memory. See the prereg.

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
LAGS = (1, 2, 3, 4, 5)
# v0.6.2 published excesses E_b = #odd - #even of s_b(g_n); reproduce EXACTLY.
EXCESS_ANCHORS = {
    1_000_000: {2: -7563, 4: 16607, 6: 3031, 8: -11497, 10: -15033, 16: -38231},
    10_000_000: {2: -77744, 4: 132234, 6: 31814, 8: -69764, 10: -96440, 16: -264548},
    100_000_000: {2: -756388, 4: 1051924, 6: 341536, 8: -437046, 10: -665788, 16: -1901898},
}
# v0.7.0 lag-1 S_b^(1) regression anchors, from corr-v1's gap_parity_corr.json
# (results[cut].bases[b].S_actual). Ties this job to the published v0.7.0 result.
LAG1_S_ANCHORS = {
    1_000_000: {2: -4572, 4: 3180, 6: -562, 8: 522, 10: 2050, 16: 17118},
    10_000_000: {2: -27941, 4: 29385, 6: -2459, 8: -2549, 10: 9321, 16: 89817},
    100_000_000: {2: -166721, 4: 233701, 6: 12329, 8: -45597, 10: 67285, 16: 508353},
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

    # ---- method gate (identical to v0.6.2 / v0.7.0) ----
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

    # ---- lag-k correlation + marginal-preserving shuffle null ----
    results = {}
    null_mean_errors = []
    lag1_regress_errors = []
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
            chi = (1 - 2 * par).astype(np.int64)  # 0->+1, 1->-1
            assert int(chi.sum()) == -snap["excess"][b], (cut, b)
            E = snap["excess"][b]
            # actual S for each lag
            s_actual = {}
            for k in LAGS:
                if N - k <= 0:
                    continue
                s_actual[k] = int(np.dot(chi[:-k], chi[k:]))
            # closed-form permutation mean per lag
            e_perm = {k: (N - k) * (E * E - N) / (N * (N - 1)) for k in LAGS if N - k > 0}
            # empirical shuffle null: one shuffle -> all lags (shared across k)
            samples = {k: np.empty(K, dtype=np.float64) for k in LAGS if N - k > 0}
            work = chi.copy()
            for j in range(K):
                rng.shuffle(work)
                for k in LAGS:
                    if N - k > 0:
                        samples[k][j] = np.dot(work[:-k], work[k:])
            lags_out = {}
            for k in LAGS:
                if N - k <= 0:
                    continue
                mean_sh = float(samples[k].mean())
                std_sh = float(samples[k].std(ddof=1))
                z = (s_actual[k] - mean_sh) / std_sh if std_sh > 0 else 0.0
                lags_out[k] = {
                    "S_actual": s_actual[k],
                    "E_perm_mean": e_perm[k],
                    "shuffle_mean": mean_sh,
                    "shuffle_std": std_sh,
                    "K": K,
                    "z": z,
                }
                if cut in (1_000_000, 10_000_000):
                    tol = 5 * std_sh / math.sqrt(K)
                    if abs(mean_sh - e_perm[k]) > max(tol, 1.0):
                        null_mean_errors.append(
                            f"null-mean mismatch cut {cut} base {b} lag {k}: emp "
                            f"{mean_sh:.2f} vs closed form {e_perm[k]:.2f} (tol {tol:.2f})")
            # v0.7.0 lag-1 regression
            if LAG1_S_ANCHORS:
                want = LAG1_S_ANCHORS.get(cut, {}).get(b)
                if want is not None and s_actual.get(1) != want:
                    lag1_regress_errors.append(
                        f"lag-1 regress cut {cut} base {b}: {s_actual.get(1)} != v0.7.0 {want}")
            per_base[b] = {"excess": E, "lags": lags_out}
        results[str(cut)] = {"cutoff": cut, "n_gaps": N, "bases": per_base}

    if null_mean_errors:
        print("GATE FAILED (null-mean identity):", flush=True)
        for e in null_mean_errors:
            print("  " + e, flush=True)
        sys.exit(1)
    if lag1_regress_errors:
        print("GATE FAILED (v0.7.0 lag-1 regression):", flush=True)
        for e in lag1_regress_errors:
            print("  " + e, flush=True)
        sys.exit(1)
    if not LAG1_S_ANCHORS:
        print("WARNING: LAG1_S_ANCHORS empty -- v0.7.0 lag-1 regression gate SKIPPED. "
              "Fill it from corr-v1 before filing.", flush=True)
    print("null-mean gate PASSED (empirical shuffle mean == (N-k)(E^2-N)/(N(N-1)) "
          "for every lag).", flush=True)
    if LAG1_S_ANCHORS:
        print("v0.7.0 lag-1 regression PASSED (S_b^(1) reproduces GAP-CORR exactly).",
              flush=True)

    os.makedirs(args.out_dir, exist_ok=True)
    out = {
        "limit": limit,
        "bases": list(BASES),
        "lags": list(LAGS),
        "seed": SEED,
        "results": results,
        "note": ("S_b^(k) = sum chi_b(g_n) chi_b(g_{n+k}) over prime gaps; "
                 "z = (S - shuffle_mean)/shuffle_std, marginal-preserving null; "
                 "E_perm_mean = (N-k)(E_b^2-N)/(N(N-1)). "
                 "See notes/gap-parity-lagk-prereg.md."),
    }
    with open(os.path.join(args.out_dir, "gap_parity_lagk.json"), "w") as f:
        json.dump(out, f, indent=2)

    print("\nlag-k gap-parity correlation (z = connected signal in shuffle-std units):")
    for cut in snap_order:
        if str(cut) not in results:
            continue
        r = results[str(cut)]
        print(f"  cutoff {cut:>12}  n_gaps={r['n_gaps']:>9}")
        for b in BASES:
            zrow = "  ".join(f"k{k}:{r['bases'][b]['lags'][k]['z']:+6.2f}"
                             for k in LAGS if k in r['bases'][b]['lags'])
            print(f"      base {b:>2}  {zrow}")
    print("\nwrote", os.path.join(args.out_dir, "gap_parity_lagk.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
