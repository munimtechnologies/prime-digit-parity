#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""Does the base-2 gap-parity "lag-2 unique-quiet" signature track ADMITTING EVEN INTEGERS,
or does it need the FULL squarefree congruence structure?

Companion to `odd_squarefree_lag_signature_job.py`, which removes the prime 2 (odd squarefree /
coprime-6 admit no even integers, forcing every gap even). This job fills the complementary
2x2 cell: lift an ODD prime from the squarefree sieve while KEEPING all the even integers.
Squarefree coprime to 3 (and to 5, and to 15) still admit even numbers, so odd gaps and the
half-shifted base-4 stride phase (floor(g/2)=(g-1)/2 for odd g) survive, yet the sieve is no
longer the full squarefree set. If the signature survives the odd-prime lift but dies under the
prime-2 lift (odd squarefree, included here as an anchor), the driver is ADMITTING EVEN
INTEGERS, not full mu^2 structure and not multiplicativity alone.

Pre-registered in prime-digit-parity/notes/admit-even-lag-signature-prereg.md (committed before
this job runs). The observable (base-2 gap digit-sum parity = Thue-Morse), the exact
marginal-preserving permutation null `null_z`, the base-4 stride split
(P_2(g)=P_4(g) XOR P_4(g>>1)), the stream analysis, the scale-free lag-2 signature detector and
the segmented squarefree-coprime sieve are reused VERBATIM (and re-self-tested) from
odd_squarefree_lag_signature_job.py. The self-test additionally asserts the mechanism contrast:
every odd-prime-lift sieve contains even members and at least one ODD gap, whereas the prime-2
lift has all even gaps. Self-tests run first; the job exits non-zero if any fails, so only
self-validated numbers are emitted.
Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
"""
import argparse
import json
import math
import os
import statistics
import sys
from array import array
from fractions import Fraction
from itertools import permutations

LAGS = [1, 2, 3, 4, 5, 6, 7, 8]


# ------------------------------------------------- observable + null (verbatim from arc)
def digit_sum_parity(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


def lag_stats(pb):
    n_eff = len(pb)
    ones = pb.count(1)
    s = n_eff - 2 * ones
    out = {}
    for k in LAGS:
        if k >= n_eff:
            out[k] = None
            continue
        a = int.from_bytes(pb[: n_eff - k], "big")
        c = int.from_bytes(pb[k:], "big")
        d_k = (a ^ c).bit_count()
        m = n_eff - k
        out[k] = m - 2 * d_k
    return n_eff, s, out


def null_z(n_eff, s, S, k):
    N_ = n_eff
    M = N_ - k
    n_shared = max(0, N_ - 2 * k)
    n_disjoint = M * (M - 1) // 2 - n_shared
    rho2 = Fraction(s * s - N_, N_ * (N_ - 1))
    rho4 = Fraction(
        s ** 4 - 6 * N_ * s * s + 8 * s * s + 3 * N_ * N_ - 6 * N_,
        N_ * (N_ - 1) * (N_ - 2) * (N_ - 3),
    )
    eS = M * rho2
    eS2 = M + 2 * n_shared * rho2 + 2 * n_disjoint * rho4
    var = eS2 - eS * eS
    varf = float(var)
    z = (S - float(eS)) / math.sqrt(varf) if varf > 0 else 0.0
    return float(eS), varf, z


# ------------------------------------------------- sieves
def _small_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    return [i for i in range(2, limit + 1) if sieve[i]]


def collect_prime_gaps(N):
    sieve = bytearray(b"\x01") * (N + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    gaps = array("H")
    prev = None
    maxg = 0
    cnt = 0
    for n in range(2, N + 1):
        if sieve[n]:
            cnt += 1
            if prev is not None:
                g = n - prev
                if g > maxg:
                    maxg = g
                gaps.append(g)
            prev = n
    del sieve
    return gaps, cnt, maxg


def collect_squarefree_coprime_gaps(N, forbidden=(), seg_size=1 << 22):
    """Segmented squarefree sieve (n squarefree iff no p^2 | n), keeping only v with
    v % p != 0 for every p in `forbidden`. forbidden=() -> all squarefree (A005117);
    forbidden=(2,) -> odd squarefree (A056911); forbidden=(3,) -> squarefree coprime to 3
    (keeps even integers); forbidden=(2,3) -> squarefree coprime to 6."""
    root = int(round(N ** 0.5)) + 2
    primes = _small_primes(root)
    pk = [p * p for p in primes if p * p <= N]  # ascending
    forb = tuple(forbidden)
    gaps = array("H")
    cnt = 0
    prev = None
    maxg = 0
    lo = 1
    while lo <= N:
        hi = min(lo + seg_size, N + 1)
        length = hi - lo
        seg = bytearray(b"\x01") * length
        for q in pk:
            if q > hi - 1:
                break
            start = ((lo + q - 1) // q) * q
            for m in range(start, hi, q):
                seg[m - lo] = 0
        for off in range(length):
            if seg[off]:
                v = lo + off
                if forb:
                    bad = False
                    for p in forb:
                        if v % p == 0:
                            bad = True
                            break
                    if bad:
                        continue
                cnt += 1
                if prev is not None:
                    g = v - prev
                    if g > maxg:
                        maxg = g
                    gaps.append(g)
                prev = v
        lo = hi
    return gaps, cnt, maxg


# ------------------------------------------------- per-sieve stream analysis (verbatim)
def stream_bit(v, which):
    if which == "a_base4_even":
        return digit_sum_parity(v, 4)
    if which == "b_base4_odd":
        return digit_sum_parity(v >> 1, 4)
    if which == "x_base2":
        return digit_sum_parity(v, 2)
    raise ValueError(which)


STREAMS = ["a_base4_even", "b_base4_odd", "x_base2"]


def analyze_sieve(gaps, max_gap):
    out = {}
    for name in STREAMS:
        tbl = bytes(stream_bit(v, name) for v in range(max_gap + 1))
        pb = bytes(tbl[g] for g in gaps)
        n_eff, s, S = lag_stats(pb)
        frac_odd = (n_eff - s) / 2 / n_eff if n_eff else 0.0
        per_lag = {}
        for k in LAGS:
            if S[k] is None:
                continue
            eS, var, z = null_z(n_eff, s, S[k], k)
            per_lag[str(k)] = {"S": S[k], "E_null": eS, "var_null": var, "z": z}
        out[name] = {"n_eff": n_eff, "s": s, "marginal_frac_odd": frac_odd,
                     "lags": per_lag}
    return out


# ------------------------------------------------- scale-free signature detector (verbatim)
def lag2_signature(x_lags):
    az = {k: abs(x_lags[str(k)]["z"]) for k in LAGS}
    quiet_lag = min(az, key=az.get)
    others = [az[k] for k in LAGS if k != 2]
    med_others = statistics.median(others)
    ratio = (az[2] / med_others) if med_others > 0 else float("inf")
    present = bool(quiet_lag == 2 and az[2] <= 0.1 * med_others)
    return {
        "abs_z_by_lag": {str(k): az[k] for k in LAGS},
        "quiet_lag": quiet_lag,
        "abs_z_lag2": az[2],
        "median_abs_z_other_lags": med_others,
        "ratio_lag2_over_median_others": ratio,
        "signature_present": present,
    }


# ------------------------------------------------- self tests
def digit_parity_tests():
    assert digit_sum_parity(19, 10) == 0
    assert digit_sum_parity(7, 2) == 1
    assert digit_sum_parity(8, 8) == 1
    assert digit_sum_parity(255, 16) == 0


def test_stride_identity():
    for g in range(0, 4096):
        a = digit_sum_parity(g, 4)
        b = digit_sum_parity(g >> 1, 4)
        assert (a ^ b) == digit_sum_parity(g, 2), g
        assert digit_sum_parity(g, 2) == bin(g).count("1") % 2, g


def test_lag_stats():
    pb = bytes([0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1])
    n_eff, s, S = lag_stats(pb)
    assert s == (pb.count(0) - pb.count(1))
    for k in LAGS:
        if k >= n_eff:
            continue
        disag = sum(1 for i in range(len(pb) - k) if pb[i] != pb[i + k])
        assert S[k] == (len(pb) - k) - 2 * disag, (k, S[k])


def test_null_formula():
    for seq in ([1, 1, 1, -1, -1, -1, 1, -1], [1, 1, 1, 1, -1, -1, 1]):
        N_ = len(seq)
        s = sum(seq)
        for k in (1, 2, 3, 4):
            tot = tot2 = cnt = 0
            for perm in permutations(seq):
                Sk = sum(perm[i] * perm[i + k] for i in range(N_ - k))
                tot += Sk
                tot2 += Sk * Sk
                cnt += 1
            emp_mean = Fraction(tot, cnt)
            emp_m2 = Fraction(tot2, cnt)
            M = N_ - k
            n_shared = max(0, N_ - 2 * k)
            n_disjoint = M * (M - 1) // 2 - n_shared
            rho2 = Fraction(s * s - N_, N_ * (N_ - 1))
            rho4 = Fraction(
                s ** 4 - 6 * N_ * s * s + 8 * s * s + 3 * N_ * N_ - 6 * N_,
                N_ * (N_ - 1) * (N_ - 2) * (N_ - 3),
            )
            assert emp_mean == M * rho2, (k, emp_mean, M * rho2)
            assert emp_m2 == M + 2 * n_shared * rho2 + 2 * n_disjoint * rho4, k


def _brute_squarefree(n):
    """Obviously-correct squarefree test by trial division."""
    if n <= 0:
        return False
    d = 2
    m = n
    while d * d <= m:
        if m % (d * d) == 0:
            return False
        while m % d == 0:
            m //= d
        d += 1
    return True


def _ref_members(limit, forbidden):
    forb = tuple(forbidden)
    out = []
    for v in range(1, limit + 1):
        if not _brute_squarefree(v):
            continue
        if any(v % p == 0 for p in forb):
            continue
        out.append(v)
    return out


def test_filtered_sieve():
    # Brute-force squarefree checker sanity.
    assert _brute_squarefree(1) and _brute_squarefree(2) and _brute_squarefree(30)
    assert not _brute_squarefree(4) and not _brute_squarefree(12)
    assert not _brute_squarefree(18) and not _brute_squarefree(50)  # 2*3^2, 2*5^2
    # Exact reference match (count, full gap stream, max gap) over [1, 20000] for each
    # variant, with a deliberately tiny segment so many segment boundaries are crossed.
    # Covers the complementary cell: odd-prime lifts (3,), (5,), (3,5) that KEEP even integers.
    LIMIT = 20000
    for forbidden in ((), (2,), (2, 3), (3,), (5,), (3, 5)):
        members = _ref_members(LIMIT, forbidden)
        ref_gaps = [members[i + 1] - members[i] for i in range(len(members) - 1)]
        gaps, cnt, maxg = collect_squarefree_coprime_gaps(LIMIT, forbidden, seg_size=97)
        assert cnt == len(members), (forbidden, cnt, len(members))
        assert list(gaps) == ref_gaps, forbidden
        assert maxg == (max(ref_gaps) if ref_gaps else 0), (forbidden, maxg)
    # MECHANISM CONTRAST (the whole point of this job):
    # - odd-prime lifts keep even members AND produce at least one ODD gap (half-shifted phase
    #   survives);
    # - the prime-2 lift (odd squarefree) forces ALL gaps even.
    for forbidden in ((3,), (5,), (3, 5)):
        members = _ref_members(20000, forbidden)
        gaps, _, _ = collect_squarefree_coprime_gaps(20000, forbidden, seg_size=97)
        assert any(v % 2 == 0 for v in members), ("no even member", forbidden)
        assert any(g % 2 == 1 for g in gaps), ("no odd gap", forbidden)
    g_odd, _, _ = collect_squarefree_coprime_gaps(20000, (2,), seg_size=97)
    assert all(x % 2 == 0 for x in g_odd), "odd squarefree must have all even gaps"
    g6, _, _ = collect_squarefree_coprime_gaps(20000, (2, 3), seg_size=97)
    assert all(x % 2 == 0 for x in g6), "coprime-6 must have all even gaps"
    # A005117 count Q(10^6) = 607926 (Gegenbauer); re-check the full squarefree count exactly.
    _, qsf, _ = collect_squarefree_coprime_gaps(10 ** 6, (), seg_size=1 << 20)
    assert qsf == 607926, qsf


def test_primes():
    _, pc, _ = collect_prime_gaps(10 ** 6)
    assert pc == 78498, pc


def test_signature_detector():
    def mk(zs):
        return {str(k): {"z": zs[k - 1]} for k in range(1, 9)}
    sq_like = mk([-502, -7, -496, -426, -441, -439, -427, -334])
    assert lag2_signature(sq_like)["signature_present"] is True
    pr_like = mk([112.8, 11.4, 0.62, 0.95, 0.19, 2.6, 2.2, 0.64])
    assert lag2_signature(pr_like)["signature_present"] is False
    weak = mk([10, 3, 12, 11, 9, 13, 8, 10])
    assert lag2_signature(weak)["signature_present"] is False


def run_selftests():
    digit_parity_tests()
    test_stride_identity()
    test_lag_stats()
    test_null_formula()
    test_filtered_sieve()
    test_primes()
    test_signature_detector()


# ------------------------------------------------- grading
def grade(res):
    sieves = list(res.keys())
    sig = {s: lag2_signature(res[s]["x_base2"]["lags"]) for s in sieves}

    # P-AE-sf (control): full squarefree shows the signature (replicates v0.18.0 at this N).
    p_sf = {"squarefree_signature_present": sig["squarefree"]["signature_present"],
            "prediction": "present (control replication)",
            "holds": bool(sig["squarefree"]["signature_present"]),
            "detail": sig["squarefree"]}

    # P-AE-c3 (PRIMARY): squarefree coprime to 3 (keeps even integers) shows the signature.
    p_c3 = {"squarefree_coprime3_signature_present": sig["squarefree_coprime3"]["signature_present"],
            "prediction": "present (keeps even integers -> odd gaps + half-shifted stride phase "
                          "survive, driving the lag-2 XOR near-null)",
            "holds": bool(sig["squarefree_coprime3"]["signature_present"]),
            "detail": sig["squarefree_coprime3"]}

    # P-AE-c5 (replication): squarefree coprime to 5 (keeps even integers) shows the signature.
    p_c5 = {"squarefree_coprime5_signature_present": sig["squarefree_coprime5"]["signature_present"],
            "prediction": "present (second odd-prime lift, still even-admitting)",
            "holds": bool(sig["squarefree_coprime5"]["signature_present"]),
            "detail": sig["squarefree_coprime5"]}

    # P-AE-c15 (replication, double lift): squarefree coprime to 15 shows the signature.
    p_c15 = {"squarefree_coprime15_signature_present": sig["squarefree_coprime15"]["signature_present"],
             "prediction": "present (two odd primes lifted, still even-admitting)",
             "holds": bool(sig["squarefree_coprime15"]["signature_present"]),
             "detail": sig["squarefree_coprime15"]}

    # P-AE-odd (anchor): odd squarefree (removes even integers) does NOT show the signature.
    p_odd = {"odd_squarefree_signature_present": sig["odd_squarefree"]["signature_present"],
             "prediction": "absent (removes even integers -> all gaps even, anchor)",
             "holds": bool(not sig["odd_squarefree"]["signature_present"]),
             "detail": sig["odd_squarefree"]}

    # Decisive read: the signature tracks admitting even integers iff all odd-prime lifts keep
    # it while the prime-2 lift loses it.
    tracks_even_admission = bool(
        sig["squarefree_coprime3"]["signature_present"]
        and sig["squarefree_coprime5"]["signature_present"]
        and sig["squarefree_coprime15"]["signature_present"]
        and not sig["odd_squarefree"]["signature_present"]
    )

    return {
        "P_AE_sf_control": p_sf,
        "P_AE_c3_primary": p_c3,
        "P_AE_c5_replication": p_c5,
        "P_AE_c15_replication": p_c15,
        "P_AE_odd_anchor": p_odd,
        "signature_tracks_even_admission_not_full_mu2": tracks_even_admission,
        "signature_all_sieves": sig,
    }


# ------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=20_000_000)
    ap.add_argument("--out-dir", default="/out")
    args = ap.parse_args()
    N = args.n_max
    os.makedirs(args.out_dir, exist_ok=True)

    print("[self-test] running ...", flush=True)
    try:
        run_selftests()
    except AssertionError as e:
        json.dump({"self_tests": {"all_passed": False, "error": str(e)}},
                  open(os.path.join(args.out_dir, "results.json"), "w"), indent=2)
        print("SELF-TEST FAILED:", e, file=sys.stderr)
        sys.exit(1)
    print("[self-test] passed", flush=True)

    res = {}
    counts = {}
    collectors = (
        ("squarefree", lambda n: collect_squarefree_coprime_gaps(n, ())),
        ("squarefree_coprime3", lambda n: collect_squarefree_coprime_gaps(n, (3,))),
        ("squarefree_coprime5", lambda n: collect_squarefree_coprime_gaps(n, (5,))),
        ("squarefree_coprime15", lambda n: collect_squarefree_coprime_gaps(n, (3, 5))),
        ("odd_squarefree", lambda n: collect_squarefree_coprime_gaps(n, (2,))),
    )
    for key, fn in collectors:
        print(f"[sieve] {key} up to {N} ...", flush=True)
        gaps, cnt, maxg = fn(N)
        print(f"[sieve] {key}: count={cnt} gaps={len(gaps)} max_gap={maxg}", flush=True)
        res[key] = analyze_sieve(gaps, maxg)
        counts[key] = {"count": cnt, "n_gaps": len(gaps), "max_gap": maxg}
        for name in STREAMS:
            zs = " ".join(
                f"z{k}={res[key][name]['lags'].get(str(k), {}).get('z', float('nan')):+9.2f}"
                for k in LAGS)
            print(f"    {key:22s} {name:14s} {zs}", flush=True)
        del gaps

    grading = grade(res)
    out = {
        "title": "Does the base-2 gap-parity lag-2 unique-quiet signature track admitting even "
                 "integers (lift an odd prime, keep the even ones) or the full squarefree "
                 "congruence structure?",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "sieves": {
            "squarefree": "A005117 (control, expect present)",
            "squarefree_coprime3": "squarefree & coprime to 3 (PRIMARY: keeps even, expect present)",
            "squarefree_coprime5": "squarefree & coprime to 5 (keeps even, expect present)",
            "squarefree_coprime15": "squarefree & coprime to 15 (keeps even, expect present)",
            "odd_squarefree": "A056911 (anchor: removes even, expect absent)",
        },
        "null": "exact marginal-preserving permutation null of S^(k), reported as z",
        "signature": "lag-2 unique-quiet: argmin_k|z_x(k)| == 2 AND |z_x(2)| <= 0.1 x "
                     "median|z_x(k!=2)| (scale-free ratio test)",
        "counts": counts,
        "per_sieve": res,
        "grading": grading,
        "self_tests": {"all_passed": True},
        "produced_by": "AI agents working as MV2, Munim, Inc.",
    }
    with open(os.path.join(args.out_dir, "results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("[grading]", json.dumps(grading, indent=2), flush=True)
    print("[done] wrote results.json", flush=True)


if __name__ == "__main__":
    main()
