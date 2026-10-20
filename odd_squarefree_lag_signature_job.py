#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Munim, Inc. Produced with AI by MV2 of Munim, Inc.
"""Does the base-2 gap-parity "lag-2 unique-quiet" signature need the squarefree sieve to
ADMIT EVEN NUMBERS, or is it robust to removing them?

Follow-up to the graded k-free result (prime-digit-parity v0.18.0), which showed the lag-2
unique-quiet XOR-cancellation signature is SPECIFIC to squarefree (k=2 / mu^2) and absent
for cubefree, 4th-power-free and primes. The mechanism question this job isolates: the full
squarefree set admits both even and odd integers, so its gaps include ODD gaps (e.g. 2->3,
gap 1). Restricting to ODD squarefree numbers (A056911) removes every even element, forcing
EVERY gap to be even; the stride identity P_2(g)=P_4(g) XOR P_4(floor(g/2)) then degenerates
(for even g, floor(g/2)=g/2 exactly), changing the phase relationship of the two base-4
stride streams that produced the lag-2 near-null. This job asks whether that restriction
KILLS the signature, and corroborates with squarefree numbers coprime to 6 (another all-even-
gap multiplicative variant that also drops multiples of 3).

Pre-registered in prime-digit-parity/notes/odd-squarefree-lag-signature-prereg.md (committed
before this job runs). The observable (base-2 gap digit-sum parity = Thue-Morse), the exact
marginal-preserving permutation null `null_z`, the base-4 stride split
(P_2(g)=P_4(g) XOR P_4(g>>1)), the stream analysis and the scale-free lag-2 signature
detector are reused VERBATIM (and re-self-tested) from kfree_lag_signature_job.py. The only
new code is a segmented squarefree sieve with a coprimality filter, re-self-tested against a
brute-force trial-division reference (full count AND the complete gap stream over [1, 20000],
across many forced segment boundaries). Self-tests run first; the job exits non-zero if any
fails, so only self-validated numbers are emitted.
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
    forbidden=(2,) -> odd squarefree (A056911); forbidden=(2,3) -> squarefree coprime to 6."""
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
    LIMIT = 20000
    for forbidden in ((), (2,), (2, 3)):
        members = _ref_members(LIMIT, forbidden)
        ref_gaps = [members[i + 1] - members[i] for i in range(len(members) - 1)]
        gaps, cnt, maxg = collect_squarefree_coprime_gaps(LIMIT, forbidden, seg_size=97)
        assert cnt == len(members), (forbidden, cnt, len(members))
        assert list(gaps) == ref_gaps, forbidden
        assert maxg == (max(ref_gaps) if ref_gaps else 0), (forbidden, maxg)
    # Anchor facts: all odd-squarefree gaps are even; so are coprime-6 gaps.
    g_odd, _, _ = collect_squarefree_coprime_gaps(20000, (2,), seg_size=97)
    assert all(x % 2 == 0 for x in g_odd)
    g6, _, _ = collect_squarefree_coprime_gaps(20000, (2, 3), seg_size=97)
    assert all(x % 2 == 0 for x in g6)
    # A005117 count Q(10^6) = 607926 (Gegenbauer), and odd squarefree are 2/3 of squarefree
    # in density; here just re-check the full squarefree count exactly.
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

    # P-OS-sf (control): full squarefree shows the signature (replicates v0.18.0 at this N).
    p_sf = {"squarefree_signature_present": sig["squarefree"]["signature_present"],
            "prediction": "present (control replication)",
            "holds": bool(sig["squarefree"]["signature_present"]),
            "detail": sig["squarefree"]}

    # P-OS-odd (PRIMARY): odd squarefree (all gaps even) does NOT show the signature.
    p_odd = {"odd_squarefree_signature_present": sig["odd_squarefree"]["signature_present"],
             "prediction": "absent (removing even numbers forces all gaps even and is "
                           "hypothesized to kill the lag-2 XOR near-null)",
             "holds": bool(not sig["odd_squarefree"]["signature_present"]),
             "detail": sig["odd_squarefree"]}

    # P-OS-coprime6 (corroborating variant): squarefree coprime to 6 does NOT show it.
    p_c6 = {"squarefree_coprime6_signature_present": sig["squarefree_coprime6"]["signature_present"],
            "prediction": "absent (second all-even-gap multiplicative variant)",
            "holds": bool(not sig["squarefree_coprime6"]["signature_present"]),
            "detail": sig["squarefree_coprime6"]}

    # P-OS-primes (negative control): primes do NOT show it.
    p_pr = {"primes_signature_present": sig["primes"]["signature_present"],
            "prediction": "absent (negative control)",
            "holds": bool(not sig["primes"]["signature_present"]),
            "detail": sig["primes"]}

    return {
        "P_OS_sf_control": p_sf,
        "P_OS_odd_primary": p_odd,
        "P_OS_coprime6_variant": p_c6,
        "P_OS_primes_negative_control": p_pr,
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
        ("odd_squarefree", lambda n: collect_squarefree_coprime_gaps(n, (2,))),
        ("squarefree_coprime6", lambda n: collect_squarefree_coprime_gaps(n, (2, 3))),
        ("primes", collect_prime_gaps),
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
        "title": "Does the base-2 gap-parity lag-2 unique-quiet signature need the squarefree "
                 "sieve to admit even numbers?",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "sieves": {
            "squarefree": "A005117 (control, expect present)",
            "odd_squarefree": "A056911 (PRIMARY: all gaps even, expect absent)",
            "squarefree_coprime6": "squarefree & coprime to 6 (variant, expect absent)",
            "primes": "A000040 (negative control, expect absent)",
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
