"""choose() must pick the same tiling however its dynamic programme is organised.

``reference_choose`` is the straightforward layer-by-layer version; the shipped one filters the
candidate amplicons once per threshold and pushes values forward, and has to agree with it on
any set of candidates, ties included.
"""

from __future__ import annotations

import bisect
import random

import pytest

from vhirprimerlab.candidates import end_loss
from vhirprimerlab.models import PrimerCandidate, TilingParams
from vhirprimerlab.tiling import _popcount, choose, geometry

PARAMS = TilingParams(
    n_amplicons=5, overlap=50, tm_min=57, tm_opt=60, tm_max=63, max_tm_diff=2,
    len_min=18, len_max=25, gc_min=40, gc_max=60, max_mismatches=2, window=5,
    degenerate=True, max_degeneracy=4, min_coverage=95,
)


def reference_choose(fwds, revs, n, geom, params, span):
    fwds = sorted(fwds, key=lambda c: c.start)
    revs = sorted(revs, key=lambda c: c.start)
    f_starts = [f.start for f in fwds]
    r_starts = [r.start for r in revs]
    by_rev = []
    for r in revs:
        lo = bisect.bisect_left(f_starts, r.end - geom.len_max + 1)
        hi = bisect.bisect_right(f_starts, r.end - geom.len_min + 1)
        amps = []
        for fi in range(lo, hi):
            f = fwds[fi]
            if abs(f.tm - r.tm) <= params.max_tm_diff:
                amps.append((fi, _popcount(f.mask & r.mask), _popcount(f.mask) + _popcount(r.mask),
                             -(f.penalty + r.penalty)))
        by_rev.append(amps)
    preds = []
    for f in fwds:
        lo = bisect.bisect_left(r_starts, f.end + 1 + params.overlap)
        hi = bisect.bisect_right(r_starts, f.start + geom.max_outer_overlap - 1)
        preds.append([ri for ri in range(lo, hi)
                      if revs[ri].end - f.start + 1 <= geom.max_outer_overlap])
    first_ok = [end_loss(f, span) <= params.max_start_loss for f in fwds]
    last_ok = [end_loss(r, span) <= params.max_end_loss for r in revs]

    def run(t):
        value = {}
        back = [{}]
        for ri, amps in enumerate(by_rev):
            for fi, cov, pm, neg in amps:
                v = (cov, pm, -fwds[fi].start, neg)
                if cov >= t and first_ok[fi] and (ri not in value or v > value[ri]):
                    value[ri] = v
                    back[0][ri] = (fi, None)
        for _ in range(n - 1):
            best_in = {}
            for fi, rs in enumerate(preds):
                cands = [(value[ri], ri) for ri in rs if ri in value]
                if cands:
                    best_in[fi] = max(cands)
            new = {}
            layer = {}
            for ri, amps in enumerate(by_rev):
                for fi, cov, pm, neg in amps:
                    if cov < t or fi not in best_in:
                        continue
                    (s, m, start, p), prev = best_in[fi]
                    v = (s + cov, m + pm, start, p + neg)
                    if ri not in new or v > new[ri]:
                        new[ri] = v
                        layer[ri] = (fi, prev)
            value = new
            back.append(layer)
            if not value:
                return None
        ends = [((v[0], v[1], v[2] - (span[1] - revs[ri].end), v[3]), ri)
                for ri, v in value.items() if last_ok[ri]]
        if not ends:
            return None
        ri = max(ends)[1]
        path = []
        for layer in reversed(back):
            fi, prev = layer[ri]
            path.append((fwds[fi], revs[ri]))
            ri = prev
        return path[::-1]

    levels = sorted({cov for amps in by_rev for _, cov, _, _ in amps})
    best = None
    lo, hi = 0, len(levels) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        path = run(levels[mid])
        if path is None:
            hi = mid - 1
        else:
            best, lo = path, mid + 1
    return best


def candidates(rng: random.Random, genome: int, density: float, bits: int):
    """Random primers on both strands; masks, Tm and penalty vary so that ties are common."""
    def make(strand: str, start: int) -> PrimerCandidate:
        length = rng.randint(18, 25)
        mask = rng.getrandbits(bits) | rng.getrandbits(bits)           # mostly set bits, some gaps
        return PrimerCandidate(strand, start, length, "A" * length, "A" * length, "A" * length,
                               round(rng.uniform(58.0, 62.0), 1), None, 50.0,
                               rng.choice([0.0, 0.5, 1.0, 2.0]), mask, [])
    fwds = [make("+", s) for s in range(0, genome - 25) if rng.random() < density]
    revs = [make("-", s) for s in range(0, genome - 25) if rng.random() < density]
    return fwds, revs


def draw(seed: int):
    rng = random.Random(seed)
    genome = rng.choice([1500, 2400, 3600])
    n = rng.randint(2, 6)
    geom = geometry(genome, n, PARAMS.overlap, 0.2, PARAMS.len_max)
    if geom is None:
        return None
    fwds, revs = candidates(rng, genome, rng.choice([0.04, 0.1, 0.25]), rng.choice([6, 10, 16]))
    params = TilingParams(**{**PARAMS.__dict__, "n_amplicons": n,
                             "max_start_loss": rng.choice([30, 60, 400]),
                             "max_end_loss": rng.choice([30, 60, 400]),
                             "max_tm_diff": rng.choice([1.0, 2.0, 4.0])})
    return fwds, revs, n, geom, params, (0, genome - 1)


@pytest.mark.parametrize("seed", range(40))
def test_same_tiling_as_the_plain_version(seed: int) -> None:
    case = draw(seed)
    if case is None:
        pytest.skip("amplicons too short for this draw")
    assert choose(*case) == reference_choose(*case)


def test_the_comparison_is_not_just_none_against_none() -> None:
    cases = [c for c in (draw(seed) for seed in range(40)) if c is not None]
    assert sum(choose(*c) is not None for c in cases) >= 8
