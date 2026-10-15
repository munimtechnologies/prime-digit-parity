# Work by Munim, Inc., produced with AI. Apache-2.0.
#
# Independent re-computation of the base-10 / base-8 prime digit-sum parity
# excess around the Gupta-Sridhar region, using a PURE-STDLIB segmented sieve
# (no numpy), to confirm the dip reported in RESULTS.md with a different
# algorithm than the numpy-based job. Records base-10 and base-8 excess at
# value cutoffs 1e8 and 1e9, and base-10 excess at count cutoffs 3e7..7e7.
# See CONFIRM.md for the results and cross-checks.

import json
import math
import os

LIMIT = 1_450_000_000          # enough to pass the first 7e7 primes
SEG = 1 << 24                  # ~16.7M per segment (bytearray, memory-light)

VALUE_CUTOFFS = [10 ** 8, 10 ** 9]
COUNT_CUTOFFS = [30_000_000, 40_000_000, 50_000_000, 60_000_000, 70_000_000]


def digit_parity_base10(n):
    s = 0
    while n:
        s += n % 10
        n //= 10
    return s & 1


def digit_parity_base8(n):
    s = 0
    while n:
        s += n & 7
        n >>= 3
    return s & 1


def main():
    out_dir = os.environ.get("OUT_DIR", "/out")
    os.makedirs(out_dir, exist_ok=True)

    root = int(math.isqrt(LIMIT)) + 1
    base_sieve = bytearray([1]) * (root + 1)
    base_sieve[0] = base_sieve[1] = 0
    for i in range(2, int(math.isqrt(root)) + 1):
        if base_sieve[i]:
            base_sieve[i * i::i] = b"\x00" * len(base_sieve[i * i::i])
    base_primes = [i for i in range(2, root + 1) if base_sieve[i]]

    n = 0                       # running prime count
    e10 = 0                     # running odd-even excess, base 10
    e8 = 0                      # running odd-even excess, base 8
    fact4_ok = True

    value_pending = sorted(VALUE_CUTOFFS)
    count_pending = sorted(COUNT_CUTOFFS)
    value_snaps = {}
    count_snaps = {}

    def handle_prime(p):
        nonlocal n, e10, e8, fact4_ok
        # value-cutoff snapshot: record state of all primes <= cutoff, before
        # counting any prime that exceeds it.
        while value_pending and p > value_pending[0]:
            c = value_pending.pop(0)
            value_snaps[c] = {"n": n, "excess_base10": e10, "excess_base8": e8}
        if p > 5 and (p % 10) % 2 == 0:
            fact4_ok = False
        e10 += 1 if digit_parity_base10(p) else -1
        e8 += 1 if digit_parity_base8(p) else -1
        n += 1
        if count_pending and n == count_pending[0]:
            c = count_pending.pop(0)
            count_snaps[c] = {"largest_prime": p, "excess_base10": e10}

    low = 0
    while low <= LIMIT and (value_pending or count_pending):
        high = min(low + SEG, LIMIT + 1)
        seg = bytearray([1]) * (high - low)
        if low == 0:
            seg[0] = 0
            if high > 1:
                seg[1] = 0
        for q in base_primes:
            if q * q >= high:
                break
            start = max(q * q, ((low + q - 1) // q) * q)
            seg[start - low::q] = b"\x00" * len(seg[start - low::q])
        for i in range(high - low):
            if seg[i]:
                handle_prime(low + i)
                if not value_pending and not count_pending:
                    break
        low = high

    result = {
        "limit_reached": low,
        "algorithm": "pure-stdlib segmented sieve (independent of numpy v3)",
        "fact4_last_digit_identity_holds": fact4_ok,
        "value_cutoffs": value_snaps,
        "count_cutoffs": count_snaps,
        "note": "Work by Munim, Inc., produced with AI. Apache-2.0.",
    }
    # add z scores for convenience
    for c, s in value_snaps.items():
        s["z_base10"] = s["excess_base10"] / math.sqrt(s["n"]) if s["n"] else 0.0
        s["z_base8"] = s["excess_base8"] / math.sqrt(s["n"]) if s["n"] else 0.0
    for c, s in count_snaps.items():
        s["z_base10"] = s["excess_base10"] / math.sqrt(c)

    with open(os.path.join(out_dir, "confirm.json"), "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
