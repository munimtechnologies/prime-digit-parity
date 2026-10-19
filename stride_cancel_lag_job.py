#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Is the base-2 stride-cancellation a lag-1 artefact, or does it persist at lags 2-3?

Pre-registered in prime-digit-parity/notes/stride-cancel-lag-prereg.md. Builds on the
graded lag-1 result (STRIDE-CANCEL-RESULTS.md, v0.13.0): for sums of two squares the base-2
gap-parity stream is quiet because its two base-4 strided streams are loud but cancel under
XOR. This job repeats the measurement at lags 1, 2 and 3 to test whether the cancellation
persists (general structural feature) or is specific to adjacent gaps.

Streams, for a sieve's consecutive gaps:
    a_n = P_4(g_n)            (even-position bits)  == the base-4 observable
    b_n = P_4(floor(g_n/2))   (odd-position bits)
    x_n = a_n XOR b_n = P_2(g_n)   (the base-2 observable)

For each sieve (primes, squarefree, sums of two squares) up to N_MAX and each lag k in
{1,2,3}, compute S^(k) = sum_n chi_n chi_{n+k} (chi = +/-1) with the EXACT
marginal-preserving permutation null, reported as a z-score. Grade the pre-registered
P-L2 / P-L3 / P-PERSIST predictions against their fixed thresholds.

Self-tests run first; the job exits non-zero if any fails, so only self-validated numbers
are emitted. Pure standard library, no network. Produced with AI by MV2, Munim, Inc.
Code Apache-2.0; result data CC BY 4.0.
"""
import argparse
import json
import math
import os
import sys
from array import array
from fractions import Fraction
from itertools import permutations

LAGS = [1, 2, 3]


# ------------------------------------------------- observable + null (verbatim from arc)
def digit_sum_parity(n, b):
    s = 0
    while n:
        s += n % b
        n //= b
    return s & 1


def lag_stats(pb):
    """From a 0/1 byte string pb, compute N_eff, s, and S^(k)=M-2*D_k for k in LAGS."""
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
    """Exact marginal-preserving permutation-null mean, variance and z of S^(k)."""
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


# ------------------------------------------------- the three sieves
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


def _small_primes(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    return [i for i in range(2, limit + 1) if sieve[i]]


def collect_squarefree_gaps(N):
    """Segmented mu^2 sieve: n is squarefree iff no p^2 divides n."""
    SEG = 1 << 22
    root = int(N ** 0.5)
    primes = _small_primes(root)
    psq = [p * p for p in primes if p * p <= N]
    gaps = array("H")
    cnt = 0
    prev = None
    maxg = 0
    lo = 1
    while lo <= N:
        hi = min(lo + SEG, N + 1)
        length = hi - lo
        seg = bytearray(b"\x01") * length
        for p2 in psq:
            if p2 > hi - 1:
                break
            start = ((lo + p2 - 1) // p2) * p2
            for m in range(start, hi, p2):
                seg[m - lo] = 0
        for off in range(length):
            if seg[off]:
                q = lo + off
                cnt += 1
                if prev is not None:
                    g = q - prev
                    if g > maxg:
                        maxg = g
                    gaps.append(g)
                prev = q
        lo = hi
    return gaps, cnt, maxg


def build_bad(limit):
    """bad[n] == 1 iff n is NOT a sum of two squares.
    n in S iff every prime p = 3 (mod 4) divides n to an even power."""
    root = int(limit ** 0.5)
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, root + 1):
        if sieve[i]:
            sieve[i * i::i] = b"\x00" * len(sieve[i * i::i])
    bad = bytearray(limit + 1)
    scratch = bytearray(limit + 1)
    for p in range(3, limit + 1):
        if not sieve[p] or p % 4 != 3:
            continue
        if p > root:
            bad[p::p] = b"\x01" * (limit // p)
        else:
            pk = p
            while pk <= limit:
                for n in range(pk, limit + 1, pk):
                    scratch[n] ^= 1
                pk *= p
            for n in range(p, limit + 1, p):
                if scratch[n]:
                    bad[n] = 1
                    scratch[n] = 0
    del sieve, scratch
    return bad


def collect_sumsq_gaps(N):
    bad = build_bad(N)
    gaps = array("H")
    cnt = 0
    prev = None
    maxg = 0
    for n in range(1, N + 1):
        if bad[n] == 0:
            cnt += 1
            if prev is not None:
                g = n - prev
                if g > maxg:
                    maxg = g
                gaps.append(g)
            prev = n
    del bad
    return gaps, cnt, maxg


# ------------------------------------------------- per-sieve stream analysis
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
            eS, var, z = null_z(n_eff, s, S[k], k)
            per_lag[str(k)] = {"S": S[k], "E_null": eS, "var_null": var, "z": z}
        out[name] = {"n_eff": n_eff, "s": s, "marginal_frac_odd": frac_odd,
                     "lags": per_lag}
    return out


# ------------------------------------------------- self tests
def test_digit_parity():
    assert digit_sum_parity(19, 10) == 0
    assert digit_sum_parity(7, 2) == 1
    assert digit_sum_parity(8, 8) == 1
    assert digit_sum_parity(255, 16) == 0
    assert digit_sum_parity(0, 10) == 0


def test_stride_identity():
    for g in range(0, 4096):
        a = digit_sum_parity(g, 4)
        b = digit_sum_parity(g >> 1, 4)
        assert (a ^ b) == digit_sum_parity(g, 2), g
        assert digit_sum_parity(g, 2) == bin(g).count("1") % 2, g


def test_lag_stats():
    """S^(k) = M - 2*disagreements_at_k, checked directly for k=1,2,3."""
    pb = bytes([0, 1, 1, 0, 0, 1, 0, 1, 0, 0])
    n_eff, s, S = lag_stats(pb)
    assert n_eff == len(pb)
    assert s == (pb.count(0) - pb.count(1))
    for k in LAGS:
        disag = sum(1 for i in range(len(pb) - k) if pb[i] != pb[i + k])
        assert S[k] == (len(pb) - k) - 2 * disag, (k, S[k])


def test_null_formula():
    """Brute-force permutation average vs formula, for k=1,2,3."""
    for seq in ([1, 1, 1, -1, -1, -1, 1, -1], [1, 1, 1, 1, -1, -1, 1]):
        N_ = len(seq)
        s = sum(seq)
        for k in (1, 2, 3):
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


def test_small_sieves():
    _, pc, _ = collect_prime_gaps(10 ** 6)
    assert pc == 78498, pc
    _, sc, _ = collect_squarefree_gaps(100)
    assert sc == 61, sc
    bad = build_bad(60)
    first = [n for n in range(1, 30) if bad[n] == 0]
    assert first[:12] == [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20], first[:12]


def run_selftests():
    test_digit_parity()
    test_stride_identity()
    test_lag_stats()
    test_null_formula()
    test_small_sieves()


# ------------------------------------------------- grading
def grade(res):
    def az(sieve, stream, k):
        return abs(res[sieve][stream]["lags"][str(k)]["z"])

    def classify(sieve, k):
        a, b, x = (az(sieve, s, k) for s in STREAMS)
        cancel = bool(max(a, b) >= 30 and x <= 30)
        falsified = bool(x > 30)
        inconclusive = bool(max(a, b) < 30 and x <= 30)
        return {"cancellation_persists": cancel, "falsified_lag_specific": falsified,
                "inconclusive_no_stream_signal": inconclusive,
                "abs_z": {"a": a, "b": b, "x_base2": x}}

    return {
        "P_L2_sumsq": classify("sums_of_two_squares", 2),
        "P_L3_sumsq": classify("sums_of_two_squares", 3),
        "P_PERSIST_primes_l2": bool(az("primes", "x_base2", 2) >= 10),
        "P_PERSIST_squarefree_l2": bool(az("squarefree", "x_base2", 2) >= 10),
        "thresholds": {"stream_loud": 30, "xor_quiet": 30, "base2_persist": 10},
    }


# ------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=100_000_000)
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
    for key, fn in (("primes", collect_prime_gaps),
                    ("squarefree", collect_squarefree_gaps),
                    ("sums_of_two_squares", collect_sumsq_gaps)):
        print(f"[sieve] {key} up to {N} ...", flush=True)
        gaps, cnt, maxg = fn(N)
        print(f"[sieve] {key}: count={cnt} gaps={len(gaps)} max_gap={maxg}", flush=True)
        res[key] = analyze_sieve(gaps, maxg)
        counts[key] = {"count": cnt, "n_gaps": len(gaps), "max_gap": maxg}
        for name in STREAMS:
            zs = " ".join(f"z{k}={res[key][name]['lags'][str(k)]['z']:+9.2f}" for k in LAGS)
            print(f"    {key:22s} {name:14s} {zs}", flush=True)
        del gaps

    grading = grade(res)
    out = {
        "title": "Stride-cancellation at lags 1-3: does base-2 XOR cancellation persist?",
        "n_max": N,
        "lags": LAGS,
        "identity": "P_2(g) = P_4(g) XOR P_4(floor(g/2))",
        "streams": {
            "a_base4_even": "P_4(g) (even-position bits) == the base-4 observable",
            "b_base4_odd": "P_4(floor(g/2)) (odd-position bits)",
            "x_base2": "P_2(g) = a XOR b (the base-2 observable)",
        },
        "null": "exact marginal-preserving permutation null of S^(k), reported as z",
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
