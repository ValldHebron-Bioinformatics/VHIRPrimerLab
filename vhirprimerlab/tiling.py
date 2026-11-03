"""Tiling the whole genome with N overlapping amplicons in two pools."""

from __future__ import annotations

import bisect
import logging
import math
from typing import Callable

from vhirprimerlab.candidates import end_loss, enumerate_primers, pick_conserved, score_primers
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.models import (
    Alignment,
    Amplicon,
    ColumnStats,
    Consensus,
    Geometry,
    PipelineError,
    PrimerCandidate,
    TilingParams,
    TilingResult,
)
from vhirprimerlab.scoring import score_pair
from vhirprimerlab.structure import DIMER_DG, dimer_dg

log = logging.getLogger("vhirprimerlab")

Progress = Callable[[str, float], None]
Path_ = list[tuple[PrimerCandidate, PrimerCandidate]]

TOLERANCES = (0.20, 0.35, 0.50)  # amplicon length tolerances tried in order
LONG_AMPLICON = 3000  # warn above this target length
MAX_DIMER_SWAPS = 20  # primers dropped at most while trying to remove dimers


def target_length(genome_len: int, n: int, overlap: int) -> int:
    """Amplicon length that tiles genome_len with n amplicons overlapping by overlap."""
    return round((genome_len + (n - 1) * overlap) / n)


def geometry(
    genome_len: int, n: int, overlap: int, tolerance: float, primer_len_max: int
) -> Geometry | None:
    """Size limits for n amplicons, or None when amplicons would be too short.

    Neighbours may overlap (outer primer edges) by at most half the shortest
    amplicon, which keeps amplicons k and k+2 (same pool) apart. That overlap must
    still hold two primers plus the required insert overlap.
    """
    target = target_length(genome_len, n, overlap)
    len_min = math.floor(target * (1 - tolerance))
    max_outer = (len_min - 1) // 2
    if max_outer < overlap + 2 * primer_len_max:
        return None
    return Geometry(
        target_len=target,
        len_min=len_min,
        len_max=math.ceil(target * (1 + tolerance)),
        tolerance=tolerance,
        max_outer_overlap=max_outer,
    )


def max_amplicons(genome_len: int, overlap: int, primer_len_max: int) -> int:
    """Largest n whose amplicons are still long enough at the default tolerance (0 if none)."""
    best = 0
    for n in range(1, genome_len + 1):
        if geometry(genome_len, n, overlap, TOLERANCES[0], primer_len_max) is None:
            break
        best = n
    return best


def _popcount(x: int) -> int:
    return bin(x).count("1")


def choose(
    fwds: list[PrimerCandidate],
    revs: list[PrimerCandidate],
    n: int,
    geom: Geometry,
    params: TilingParams,
    genome_len: int,
) -> Path_ | None:
    """Exactly n amplicons in genome order maximising the worst coverage, then the sum.

    Ties are broken by how many sequences each primer matches on its own (a primer
    beside one that misses sequences, e.g. ones truncated at the genome end, should
    still match them), then by how close the tiling reaches both genome ends, then by
    the summed primer penalties.

    Dynamic programming over layers k = 1..n whose state is the reverse primer that
    ends amplicon k; feasibility for a coverage threshold is monotone, so the
    threshold is binary-searched over the distinct amplicon coverages.
    """
    fwds = sorted(fwds, key=lambda c: c.start)
    revs = sorted(revs, key=lambda c: c.start)
    f_starts = [f.start for f in fwds]
    r_starts = [r.start for r in revs]

    # candidate amplicons, grouped by reverse primer: r -> [(f, cov, primer matches, -penalty)]
    by_rev: list[list[tuple[int, int, float]]] = []
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

    # reverse primers that may end the amplicon before one starting with fwds[fi]
    preds: list[list[int]] = []
    for f in fwds:
        lo = bisect.bisect_left(r_starts, f.end + 1 + params.overlap)
        hi = bisect.bisect_right(r_starts, f.start + geom.max_outer_overlap - 1)
        preds.append([ri for ri in range(lo, hi)
                      if revs[ri].end - f.start + 1 <= geom.max_outer_overlap])

    first_ok = [end_loss(f, genome_len) <= params.max_end_loss for f in fwds]
    last_ok = [end_loss(r, genome_len) <= params.max_end_loss for r in revs]

    def run(t: int) -> Path_ | None:
        # layer value: (sum of coverage, sum of primer matches, -bases before the first
        # primer, -sum of penalty); back[k][ri] = (fi, prev ri)
        value: dict[int, tuple[int, int, int, float]] = {}
        back: list[dict[int, tuple[int, int | None]]] = [{}]
        for ri, amps in enumerate(by_rev):
            for fi, cov, pm, neg in amps:
                v = (cov, pm, -fwds[fi].start, neg)
                if cov >= t and first_ok[fi] and (ri not in value or v > value[ri]):
                    value[ri] = v
                    back[0][ri] = (fi, None)
        for _ in range(n - 1):
            best_in: dict[int, tuple[tuple[int, int, int, float], int]] = {}
            for fi, rs in enumerate(preds):
                cands = [(value[ri], ri) for ri in rs if ri in value]
                if cands:
                    best_in[fi] = max(cands)
            new: dict[int, tuple[int, int, int, float]] = {}
            layer: dict[int, tuple[int, int | None]] = {}
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
        ends = [((v[0], v[1], v[2] - (genome_len - 1 - revs[ri].end), v[3]), ri)
                for ri, v in value.items() if last_ok[ri]]
        if not ends:
            return None
        ri = max(ends)[1]
        path: Path_ = []
        for layer in reversed(back):
            fi, prev = layer[ri]
            path.append((fwds[fi], revs[ri]))
            ri = prev
        return path[::-1]

    levels = sorted({cov for amps in by_rev for _, cov, _, _ in amps})
    best: Path_ | None = None
    lo, hi = 0, len(levels) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        path = run(levels[mid])
        if path is None:
            hi = mid - 1
        else:
            best, lo = path, mid + 1
    return best


def _primer_names(index: int) -> tuple[str, str]:
    return f"VHIRPrimerLab_{index}_LEFT", f"VHIRPrimerLab_{index}_RIGHT"


Key = tuple[str, int, int]  # (strand, start, length) identifies a candidate primer


def _key(c: PrimerCandidate) -> Key:
    return c.strand, c.start, c.length


def _dimer_pairs(path: Path_) -> list[tuple[int, str, str, float, PrimerCandidate, PrimerCandidate]]:
    """Primer pairs in one pool whose worst variant pair has ΔG below DIMER_DG.

    Amplicon k (1-based) is in pool 1 when k is odd, pool 2 when even. Every variant
    of a degenerate primer is checked against every variant of the other primer.
    """
    found = []
    for pool in (1, 2):
        oligos: list[tuple[str, PrimerCandidate]] = []
        for k, (f, r) in enumerate(path, start=1):
            if (1 if k % 2 else 2) == pool:
                left, right = _primer_names(k)
                oligos += [(left, f), (right, r)]
        for i, (name_a, a) in enumerate(oligos):
            for name_b, b in oligos[i + 1:]:
                dg = dimer_dg(a.seq, b.seq)
                if dg < DIMER_DG:
                    found.append((pool, name_a, name_b, dg, a, b))
    return found


def _dimers(path: Path_) -> list[str]:
    return [f"pool {pool}: {na} × {nb} ΔG {dg / 1000:.1f} kcal/mol"
            for pool, na, nb, dg, _, _ in _dimer_pairs(path)]


def _worst_coverage(path: Path_) -> int:
    return min(_popcount(f.mask & r.mask) for f, r in path)


def _choose_any(fwds, revs, n, geoms, params, genome_len, banned: set[Key]):
    """choose() at the first length tolerance that works, skipping banned primers."""
    fwds = [c for c in fwds if _key(c) not in banned]
    revs = [c for c in revs if _key(c) not in banned]
    for geom in geoms:
        path = choose(fwds, revs, n, geom, params, genome_len)
        if path:
            return path, geom
    return None, geoms[-1]


def _avoid_dimers(path: Path_, geom: Geometry, fwds, revs, n, geoms, params, genome_len):
    """Swap out primers that form dimers within a pool, when that costs no coverage.

    For the first dimer pair not yet given up on, the tiling is chosen again without
    one of its primers (the one in most dimer pairs first). The new tiling is kept when
    it has fewer dimer pairs and its weakest amplicon covers at least as many
    sequences as before; if neither primer can go, the pair is kept and reported.
    """
    banned: set[Key] = set()
    kept: set[tuple[Key, Key]] = set()
    pairs = _dimer_pairs(path)
    swaps = 0
    while swaps < MAX_DIMER_SWAPS:
        todo = [p for p in pairs if (_key(p[4]), _key(p[5])) not in kept]
        if not todo:
            break
        counts: dict[Key, int] = {}
        for p in pairs:
            for c in (p[4], p[5]):
                counts[_key(c)] = counts.get(_key(c), 0) + 1
        a, b = todo[0][4], todo[0][5]
        order = sorted((a, b), key=lambda c: (-counts[_key(c)], -c.penalty))
        for c in order:
            trial, trial_geom = _choose_any(fwds, revs, n, geoms, params, genome_len,
                                            banned | {_key(c)})
            if trial is None or _worst_coverage(trial) < _worst_coverage(path):
                continue
            trial_pairs = _dimer_pairs(trial)
            if len(trial_pairs) < len(pairs):
                banned.add(_key(c))
                path, geom, pairs = trial, trial_geom, trial_pairs
                swaps += 1
                log.info("Dropped a primer that formed a dimer in its pool (%d dimer pairs left)",
                         len(pairs))
                break
        else:
            kept.add((_key(a), _key(b)))
            swaps += 1
    return path, geom


def _ref(cmap: CoordMap, cons: Consensus, pos: int, side: str) -> int | None:
    try:
        return cmap.col_to_ref(cons.cons_to_col[pos], side)  # type: ignore[arg-type]
    except PipelineError:
        return None


def _no_fit(
    fwds, revs, n: int, genome_len: int, params: TilingParams, cons: Consensus, cmap: CoordMap
) -> PipelineError:
    """Explain why no tiling exists: a genome end out of reach, or the longest primer-free
    stretch and N values that work."""
    for side, cands in (("start", fwds), ("end", revs)):
        closest = min((end_loss(c, genome_len) for c in cands), default=genome_len)
        if closest > params.max_end_loss:
            return PipelineError(
                f"No tiling with {n} amplicons: no acceptable primer within "
                f"{params.max_end_loss} nt of the genome {side} (the closest leaves {closest:,} nt "
                f"out); raise --max-end-loss"
            )
    sites = sorted((c.start, c.end) for c in fwds + revs)
    gap, prev_end = (0, -1, -1), -1
    for start, end in sites + [(genome_len, genome_len)]:
        if start - prev_end - 1 > gap[0]:
            gap = (start - prev_end - 1, prev_end + 1, start - 1)
        prev_end = max(prev_end, end)
    a = _ref(cmap, cons, gap[1], "start") or gap[1] + 1
    b = _ref(cmap, cons, gap[2], "end") or gap[2] + 1
    limit = max_amplicons(genome_len, params.overlap, params.len_max)
    # a primer-free stretch is bridged by fewer, longer amplicons: look below n first
    tries = [n - 1, n - 2, n - 3, n // 2, n // 3, n // 4, n + 1, n + 2, n + 3]
    works: list[str] = []
    for m in dict.fromkeys(m for m in tries if 1 <= m <= limit and m != n):
        g = geometry(genome_len, m, params.overlap, TOLERANCES[-1], params.len_max)
        if g and choose(fwds, revs, m, g, params, genome_len):
            works.append(str(m))
            if len(works) == 3:
                break
    hint = ", ".join(works) if works else "a different value"
    return PipelineError(
        f"No tiling with {n} amplicons: no acceptable primers between {a:,} and {b:,}; try -n {hint}"
    )


def tile(
    aln: Alignment,
    cons: Consensus,
    stats: list[ColumnStats],
    cmap: CoordMap,
    params: TilingParams,
    progress: Progress | None = None,
) -> TilingResult:
    """Choose params.n_amplicons overlapping amplicons covering the whole consensus."""
    report = progress or (lambda stage, fraction: None)
    n, genome_len = params.n_amplicons, len(cons.seq)
    geoms = [g for tol in TOLERANCES
             if (g := geometry(genome_len, n, params.overlap, tol, params.len_max))]
    if not geoms:
        raise PipelineError(
            f"N too large: amplicons would be ~{target_length(genome_len, n, params.overlap)} bp; "
            f"use at most {max_amplicons(genome_len, params.overlap, params.len_max)}"
        )
    if geoms[0].target_len > LONG_AMPLICON:
        log.warning("Amplicons will be ~%d bp; long PCRs are harder (try a larger -n)",
                    geoms[0].target_len)

    report("finding primers", 0.10)
    raw_f, raw_r = enumerate_primers(cons.seq, params)
    picked = pick_conserved(raw_f, stats, cons, params) + pick_conserved(raw_r, stats, cons, params)
    scored = score_primers(picked, aln, cons, params,
                           progress=lambda stage, f: report(stage, 0.15 + 0.65 * f))
    fwds = [c for c in scored if c.strand == "+"]
    revs = [c for c in scored if c.strand == "-"]
    log.info("Candidate primers: %d forward, %d reverse", len(fwds), len(revs))

    report("choosing amplicons", 0.85)
    path, geom = _choose_any(fwds, revs, n, geoms, params, genome_len, set())
    if not path:
        raise _no_fit(fwds, revs, n, genome_len, params, cons, cmap)
    report("checking primer dimers", 0.90)
    path, geom = _avoid_dimers(path, geom, fwds, revs, n, geoms, params, genome_len)
    if geom.tolerance != TOLERANCES[0]:
        log.warning("Used ±%d %% amplicon length tolerance (no tiling within ±%d %%)",
                    round(geom.tolerance * 100), round(TOLERANCES[0] * 100))

    amplicons: list[Amplicon] = []
    for k, (f, r) in enumerate(path, start=1):
        cov = score_pair(f.site, f.cols, r.site, r.cols, aln, params.scoring)
        nxt = path[k][0] if k < len(path) else None
        amplicons.append(Amplicon(
            index=k, pool=1 if k % 2 else 2, fwd=f, rev=r, coverage=cov,
            below_threshold=cov.coverage_pct < params.min_coverage,
            overlap_next=None if nxt is None else r.start - nxt.end - 1,
        ))
    first = _ref(cmap, cons, path[0][0].end + 1, "start")
    last = _ref(cmap, cons, path[-1][1].start - 1, "end")
    ends = (first - 1 if first else 0, cmap.ref_length - last if last else 0)
    return TilingResult(amplicons, geom, len(fwds), len(revs), _dimers(path), ends)
