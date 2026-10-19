"""Cross-sequence Markov-order test: is higher-order gap-parity memory prime-specific?

Run 97 established (prime-digit-parity v0.14.0) that the base-b digit-sum parity of
consecutive PRIME gaps is genuinely higher-order Markov, not first-order, in all six
bases {2,4,8,10,12,16}: the triple conditional-independence test H0 (s_{n+1} _||_
s_{n-1} | s_n) is rejected globally and survives a blocked drift control. That test
has only ever been run on primes.

This job repeats the *exact* same triple test on three sequences at N = 1e8:
  * primes                 (A000040)  -- re-run as an anchor / internal replication,
  * squarefree numbers     (A005117),
  * sums of two squares    (A001481, taken from 1).
The arc already showed the lag-1 serial memory itself is NOT prime-specific (it is
present and often stronger in these sieves) and that each sequence has a different
dominant base (no universal default base). The open question graded here is whether
the MARKOV ORDER is also universal: do all three reject first-order, and is base 2
still the loudest, or is the higher-order memory's base profile sequence-specific too?

Pre-registered in prime-digit-parity/notes/markov-order-cross-prereg.md (committed
before this job ran). Numpy provides the fast paths; every fast path is cross-checked
inside run_selftests() against the proven pure-Python markov_order.py reference and,
for the two sieves, against proven stdlib reference sieves (the squarefree-lag arc's
build_bad and a brute squarefree check). Self-tests run first and the job exits
non-zero if any fails, so only self-validated numbers are ever emitted.

Produced with AI by MV2, Munim, Inc. Code Apache-2.0; result data CC BY 4.0.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

import markov_order as mo

BASES = [2, 4, 8, 10, 12, 16]


# --------------------------------------------------------------- numpy sieves
def sieve_np(n):
    """All primes <= n as an int64 array (Sieve of Eratosthenes)."""
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].astype(np.int64)


def squarefree_np(n):
    """Squarefree numbers in 1..n (A005117): n is squarefree iff no p^2 | n."""
    sf = np.ones(n + 1, dtype=bool)
    sf[0] = False
    r = int(n ** 0.5)
    for d in range(2, r + 1):
        d2 = d * d
        if d2 > n:
            break
        sf[d2::d2] = False  # marking by every d (not just primes) is redundant but exact
    return np.nonzero(sf)[0].astype(np.int64)


def sumsq_bad_np(n, primes):
    """bad[m]==True iff m is NOT a sum of two squares.

    m in A001481 iff every prime p = 3 (mod 4) divides m to an even power.
    For primes p > sqrt(n) every multiple carries p to the first (odd) power, so all
    their multiples are bad; for small p we accumulate the exponent parity exactly.
    """
    bad = np.zeros(n + 1, dtype=bool)
    r = int(n ** 0.5)
    p3 = primes[(primes % 4) == 3]
    small = p3[p3 <= r]
    large = p3[p3 > r]
    for p in large.tolist():
        bad[p::p] = True
    par = np.zeros(n + 1, dtype=np.int8)
    for p in small.tolist():
        par.fill(0)
        pk = p
        while pk <= n:
            par[pk::pk] ^= 1
            pk *= p
        bad |= (par == 1)
    return bad


def sumsq_np(n, primes):
    """Sums of two squares in 1..n (A001481 from 1)."""
    bad = sumsq_bad_np(n, primes)
    ss = np.nonzero(~bad)[0]
    return ss[ss >= 1].astype(np.int64)


# --------------------------------------------------------------- markov analysis (numpy)
def np_triple_counts(bits):
    b = bits.astype(np.int64)
    code = 4 * b[:-2] + 2 * b[1:-1] + b[2:]
    return np.bincount(code, minlength=8)


def np_lag_k(bits, k):
    u = 1.0 - 2.0 * bits.astype(np.float64)
    u -= u.mean()
    denom = float(np.dot(u, u))
    if denom == 0.0:
        return 0.0
    return float(np.dot(u[:-k], u[k:]) / denom)


def blocked_g2(bits, blocks):
    """Sum of the per-block triple G2 statistics; df = 2 * (#usable blocks)."""
    n = bits.size
    edges = np.linspace(0, n, blocks + 1, dtype=np.int64)
    g2_sum = 0.0
    used = 0
    for i in range(blocks):
        seg = bits[edges[i]:edges[i + 1]]
        if seg.size < 3:
            continue
        res = mo.g_test_triples(np_triple_counts(seg).tolist())
        g2_sum += res["G2"]
        used += 1
    df = 2 * used
    p = mo.chi2_sf(g2_sum, df) if df > 0 else 1.0
    # how many sigma the blocked G2 sits above its null mean (= df, sd = sqrt(2 df))
    sigma = (g2_sum - df) / ((2.0 * df) ** 0.5) if df > 0 else 0.0
    return {"G2": g2_sum, "df": df, "p": p, "blocks_used": used, "sigma_over_df": sigma}


def analyze_sequence(gaps, blocks):
    max_gap = int(gaps.max())
    per_base = {}
    for base in BASES:
        tbl = np.array(mo.parity_table(max_gap, base), dtype=np.int8)
        bits = tbl[gaps]
        counts = np_triple_counts(bits)
        glob = mo.g_test_triples(counts.tolist())
        blk = blocked_g2(bits, blocks)
        r1 = np_lag_k(bits, 1)
        r2 = np_lag_k(bits, 2)
        r1sq = r1 * r1
        per_base[str(base)] = {
            "triple_counts": counts.tolist(),
            "n_triples": int(bits.size - 2),
            "global": glob,
            "blocked": blk,
            "r1": r1,
            "r2": r2,
            "r1_squared": r1sq,
            "r2_minus_r1sq": r2 - r1sq,
            "frac_ones": float(bits.mean()),
        }
    g2_by_base = {b: per_base[b]["global"]["G2"] for b in per_base}
    argmax_base = max(g2_by_base, key=g2_by_base.get)
    return {
        "max_gap": max_gap,
        "per_base": per_base,
        "global_G2_by_base": g2_by_base,
        "base_with_largest_global_G2": argmax_base,
    }


# --------------------------------------------------------------- stdlib references (self-test only)
def _stdlib_build_bad(limit):
    """Proven sums-of-two-squares reference from the squarefree-lag arc (stdlib)."""
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
                for m in range(pk, limit + 1, pk):
                    scratch[m] ^= 1
                pk *= p
            for m in range(p, limit + 1, p):
                if scratch[m]:
                    bad[m] = 1
                    scratch[m] = 0
    return bad


def _brute_squarefree(m):
    """True iff m (>=1) has no repeated prime factor; trial division reference."""
    if m < 1:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            m //= d
            if m % d == 0:   # d divides m at least twice -> d^2 | original m
                return False
        d += 1
    return True


# --------------------------------------------------------------- self tests
def run_selftests():
    checks = {}
    # 1. prime counts
    checks["pi_1e6"] = int(sieve_np(1_000_000).size) == 78498
    checks["pi_1e4"] = int(sieve_np(10_000).size) == 1229
    # 2. squarefree count to 100 (known: 61) and brute agreement to 2000
    checks["sf_count_100"] = int(squarefree_np(100).size) == 61
    sf_small = set(squarefree_np(2000).tolist())
    checks["sf_brute_2000"] = all(
        (m in sf_small) == _brute_squarefree(m) for m in range(1, 2001))
    # 3. sums of two squares: first 12 and exact bad-array agreement vs stdlib to 20000
    first12 = sumsq_np(60, sieve_np(60)).tolist()[:12]
    checks["sumsq_first12"] = first12 == [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20]
    LIM = 20000
    bad_np = sumsq_bad_np(LIM, sieve_np(LIM))
    bad_ref = _stdlib_build_bad(LIM)
    checks["sumsq_bad_matches_stdlib"] = all(
        bool(bad_np[m]) == bool(bad_ref[m]) for m in range(LIM + 1))
    # 4. numpy triple counts / autocorr cross-check vs pure-Python reference
    rng = np.random.default_rng(12345)
    rb = rng.integers(0, 2, size=5000).astype(np.int8)
    checks["triple_crosscheck"] = (
        np_triple_counts(rb).tolist() == mo.triple_counts(rb.tolist()))
    checks["lag1_crosscheck"] = abs(
        np_lag_k(rb, 1) - mo.lag_k_autocorr(rb.tolist(), 1)) < 1e-9
    checks["lag2_crosscheck"] = abs(
        np_lag_k(rb, 2) - mo.lag_k_autocorr(rb.tolist(), 2)) < 1e-9
    # exact autocorr value on a hand-computed case: bits=[0,0,1,1] -> r1 = 0.25
    checks["lag_exact"] = abs(mo.lag_k_autocorr([0, 0, 1, 1], 1) - 0.25) < 1e-12
    # 5. chi2 survival closed forms
    checks["chi2_df2"] = all(
        abs(mo.chi2_sf(x, 2) - np.exp(-x / 2.0)) < 1e-9 for x in (1.0, 5.0, 30.0))
    checks["chi2_df200"] = 0.0 < mo.chi2_sf(220.0, 200) < 1.0
    # 6. G-test gates: a conditionally-independent table gives G2~0; a dependent one is large.
    #    counts indexed 4a+2b+c. Independent given b: each 2x2 (a,c) table rank-1.
    indep = [0] * 8
    for b in (0, 1):
        ra = (100, 300) if b == 0 else (50, 100)   # row weights for a=0,1
        cc = (1, 2) if b == 0 else (1, 3)           # col weights for c=0,1
        for a in (0, 1):
            for c in (0, 1):
                indep[4 * a + 2 * b + c] = ra[a] * cc[c]  # rank-1 => exactly independent
    g_indep = mo.g_test_triples(indep)
    checks["g_independent_small"] = g_indep["G2"] < 1e-6
    dep = [0] * 8
    for b in (0, 1):
        # strong dependence: mass only on a==c
        dep[4 * 0 + 2 * b + 0] = 1000
        dep[4 * 1 + 2 * b + 1] = 1000
    g_dep = mo.g_test_triples(dep)
    checks["g_dependent_large"] = g_dep["G2"] > 100.0
    checks["g_independent_G2"] = g_indep["G2"]
    checks["g_dependent_G2"] = g_dep["G2"]
    passed = all(v is True for k, v in checks.items() if not k.endswith("_G2"))
    checks["all_passed"] = passed
    return checks


# --------------------------------------------------------------- grading
def grade(per_seq):
    seqs = ["primes", "squarefree", "sums_of_two_squares"]

    def base2_global(s):
        return per_seq[s]["per_base"]["2"]["global"]["G2"]

    def base2_sigma(s):
        return per_seq[s]["per_base"]["2"]["blocked"]["sigma_over_df"]

    def all_bases_reject(s):
        return all(per_seq[s]["per_base"][str(b)]["global"]["G2"] >= 40.0 for b in BASES)

    p1 = {s: base2_global(s) >= 100.0 for s in seqs}
    p2 = {s: base2_sigma(s) >= 5.0 for s in seqs}
    p3 = {s: all_bases_reject(s) for s in seqs}
    argmax = {s: per_seq[s]["base_with_largest_global_G2"] for s in seqs}
    distinct_argmax = len(set(argmax.values())) >= 2
    return {
        "P_MC_1_base2_rejects_firstorder_global": {
            "per_sequence": p1, "holds": all(p1.values()),
            "detail_G2": {s: base2_global(s) for s in seqs}},
        "P_MC_2_base2_survives_blocking": {
            "per_sequence": p2, "holds": all(p2.values()),
            "detail_sigma_over_df": {s: base2_sigma(s) for s in seqs}},
        "P_MC_3_all_six_bases_reject": {
            "per_sequence": p3, "holds": all(p3.values())},
        "P_MC_4_base_profile_sequence_specific": {
            "argmax_base_per_sequence": argmax,
            "holds": distinct_argmax,
            "note": "holds if the loudest base is NOT the same across all three sequences"},
        "thresholds": {"reject_global_G2": 40.0, "base2_global_G2": 100.0,
                       "blocking_sigma": 5.0},
    }


# --------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=100_000_000)
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--blocks", type=int, default=100)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, "results.json")

    t0 = time.time()
    print("[self-test] running ...", flush=True)
    selftests = run_selftests()
    if not selftests["all_passed"]:
        json.dump({"self_tests": selftests, "error": "self-tests failed"},
                  open(out_path, "w"), indent=2)
        print("SELF-TESTS FAILED:", selftests, file=sys.stderr)
        sys.exit(1)
    print("[self-test] passed", flush=True)

    N = args.n_max
    primes = sieve_np(N)

    per_seq = {}
    counts = {}

    def record(key, nums):
        gaps = np.diff(nums).astype(np.int64)
        per_seq[key] = analyze_sequence(gaps, args.blocks)
        counts[key] = {"count": int(nums.size), "n_gaps": int(gaps.size),
                       "max_gap": int(gaps.max())}
        print(f"[seq] {key}: count={nums.size} gaps={gaps.size} "
              f"max_gap={int(gaps.max())} "
              f"G2_base2={per_seq[key]['per_base']['2']['global']['G2']:.1f} "
              f"argmax_base={per_seq[key]['base_with_largest_global_G2']}", flush=True)
        del gaps

    print("[seq] primes ...", flush=True)
    record("primes", primes)
    print("[seq] squarefree ...", flush=True)
    record("squarefree", squarefree_np(N))
    print("[seq] sums_of_two_squares ...", flush=True)
    record("sums_of_two_squares", sumsq_np(N, primes))

    grading = grade(per_seq)
    results = {
        "title": "Cross-sequence Markov-order test of gap digit-sum parity "
                 "(primes vs squarefree vs sums of two squares)",
        "question": "Is the higher-order (non-first-order) Markov memory in base-b gap "
                    "digit-sum parity prime-specific, or universal across sieves -- and "
                    "is base 2 always the loudest?",
        "n_max": N,
        "bases": BASES,
        "blocks": args.blocks,
        "H0": "first-order Markov: s_{n+1} independent of s_{n-1} given s_n (triple test, df=2)",
        "counts": counts,
        "per_sequence": per_seq,
        "grading": grading,
        "self_tests": selftests,
        "elapsed_sec": round(time.time() - t0, 2),
        "produced_by": "AI agents working as MV2, Munim, Inc.",
    }
    json.dump(results, open(out_path, "w"), indent=2)
    print("[grading]", json.dumps(grading, indent=2), flush=True)
    print(f"[done] wrote results.json ({results['elapsed_sec']}s)", flush=True)


if __name__ == "__main__":
    main()
