"""Vectorised computation of digit-sum parity excess over primes (numpy).

A segmented sieve streams primes up to `limit` in blocks; for each even base we
compute the base-b digit-sum parity of every prime with a vectorised digit peel,
and record the running odd-minus-even excess at two kinds of checkpoints:

  * VALUE cutoffs x = 10**k (and a few linear points): primes p <= x,
  * COUNT cutoffs N = m * step: the first N primes.

The distinction is the whole point of the study (see NOTES.md): integer
digit-sum parity is exactly balanced at value cutoffs x = b**k but not at an
arbitrary prime-count cutoff, so the two can disagree.

Correctness of parity_base / the sieve is cross-checked against core.py in
test_core.py on small ranges.

Work by Munim, Inc.  Produced with AI (see README.md).  Apache-2.0.
"""

from __future__ import annotations

import numpy as np


def _base_primes(n: int) -> np.ndarray:
    sieve = np.ones(n + 1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if sieve[i]:
            sieve[i * i :: i] = False
    return np.nonzero(sieve)[0].astype(np.int64)


def segmented_primes(limit: int, seg: int = 10_000_000):
    """Yield numpy int64 arrays of primes, in order, covering 2..limit."""
    root = int(limit**0.5) + 1
    bp = _base_primes(root)
    yield bp[bp <= limit]
    low = root + 1
    while low <= limit:
        high = min(low + seg - 1, limit)
        size = high - low + 1
        block = np.ones(size, dtype=bool)
        for p in bp:
            if p * p > high:
                break
            start = max(p * p, ((low + p - 1) // p) * p)
            block[start - low :: p] = False
        yield np.nonzero(block)[0].astype(np.int64) + low
        low = high + 1


def parity_base(arr: np.ndarray, base: int) -> np.ndarray:
    """Digit-sum parity (0/1, int8) of every entry of `arr` in the given base."""
    x = arr.copy()
    par = np.zeros(len(arr), dtype=np.int8)
    while x.any():
        par ^= (x % base % 2).astype(np.int8)
        x //= base
    return par


def run(limit: int, bases, value_checkpoints, count_step: int, seg: int = 10_000_000):
    """Return a results dict: checkpoints with per-base odd counts and excess."""
    bases = list(bases)
    vpts = sorted(int(v) for v in value_checkpoints if v <= limit)
    vi = 0  # index into vpts of next pending value checkpoint

    cum_count = 0
    cum_odd = {b: 0 for b in bases}
    value_rows = []
    count_rows = []
    seg_rows = []
    next_count = count_step

    for pr in segmented_primes(limit, seg):
        if len(pr) == 0:
            continue
        L = len(pr)
        # per-base parity and within-segment prefix sums of "odd"
        par = {b: parity_base(pr, b) for b in bases}
        pref = {b: np.cumsum(par[b], dtype=np.int64) for b in bases}
        seg_max = int(pr[-1])

        # ---- value checkpoints that fall at or below this segment's max ----
        while vi < len(vpts) and vpts[vi] <= seg_max:
            v = vpts[vi]
            idx = int(np.searchsorted(pr, v, side="right"))  # primes in seg <= v
            n = cum_count + idx
            row = {"x": v, "n": n, "odd": {}}
            for b in bases:
                odd_in = int(pref[b][idx - 1]) if idx > 0 else 0
                odd = cum_odd[b] + odd_in
                row["odd"][b] = odd
                row.setdefault("excess", {})[b] = 2 * odd - n
            value_rows.append(row)
            vi += 1

        # ---- count checkpoints that fall inside this segment ----
        while next_count <= cum_count + L:
            idx = next_count - cum_count  # 1..L
            p_val = int(pr[idx - 1])
            row = {"n": next_count, "x": p_val, "odd": {}, "excess": {}}
            for b in bases:
                odd = cum_odd[b] + int(pref[b][idx - 1])
                row["odd"][b] = odd
                row["excess"][b] = 2 * odd - next_count
            count_rows.append(row)
            next_count += count_step

        # advance running totals
        cum_count += L
        for b in bases:
            cum_odd[b] += int(pref[b][-1])

        # trajectory: running excess at the end of every segment
        seg_rows.append(
            {
                "x": seg_max,
                "n": cum_count,
                "excess": {b: 2 * cum_odd[b] - cum_count for b in bases},
            }
        )

    # value checkpoints beyond the largest prime (e.g. x = limit = 10**k, which is
    # composite, so no segment's last prime reaches it): totals are final
    while vi < len(vpts):
        v = vpts[vi]
        row = {"x": v, "n": cum_count, "odd": {}, "excess": {}}
        for b in bases:
            row["odd"][b] = cum_odd[b]
            row["excess"][b] = 2 * cum_odd[b] - cum_count
        value_rows.append(row)
        vi += 1

    # final totals as one more row
    final = {"x": limit, "n": cum_count, "odd": {}, "excess": {}}
    for b in bases:
        final["odd"][b] = cum_odd[b]
        final["excess"][b] = 2 * cum_odd[b] - cum_count

    return {
        "limit": limit,
        "bases": bases,
        "value_checkpoints": value_rows,
        "count_checkpoints": count_rows,
        "segment_checkpoints": seg_rows,
        "final": final,
    }
