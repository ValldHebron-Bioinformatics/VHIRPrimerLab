"""Tiling the whole genome with N overlapping amplicons in two pools."""

from __future__ import annotations

import bisect
import logging
import math
from typing import Callable

import primer3

from vhirprimerlab.candidates import enumerate_primers, pick_conserved, score_primers
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

log = logging.getLogger("vhirprimerlab")

Progress = Callable[[str, float], None]
Path_ = list[tuple[PrimerCandidate, PrimerCandidate]]

TOLERANCES = (0.20, 0.35, 0.50)  # amplicon length tolerances tried in order
LONG_AMPLICON = 3000  # warn above this target length
DIMER_DG = -9000.0  # cal/mol; primer pairs in one pool below this are reported


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
        end_slack=target // 4,
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

    Ties are broken by how close the tiling reaches both genome ends, then by the
    summed primer penalties.

    Dynamic programming over layers k = 1..n whose state is the reverse primer that
    ends amplicon k; feasibility for a coverage threshold is monotone, so the
    threshold is binary-searched over the distinct amplicon coverages.
    """
    fwds = sorted(fwds, key=lambda c: c.start)
    revs = sorted(revs, key=lambda c: c.start)
    f_starts = [f.start for f in fwds]
    r_starts = [r.start for r in revs]

    # candidate amplicons, grouped by reverse primer: r -> [(f, cov, -penalty)]
    by_rev: list[list[tuple[int, int, float]]] = []
    for r in revs:
        lo = bisect.bisect_left(f_starts, r.end - geom.len_max + 1)
        hi = bisect.bisect_right(f_starts, r.end - geom.len_min + 1)
        amps = []
        for fi in range(lo, hi):
            f = fwds[fi]
            if abs(f.tm - r.tm) <= params.max_tm_diff:
                amps.append((fi, _popcount(f.mask & r.mask), -(f.penalty + r.penalty)))
        by_rev.append(amps)

    # reverse primers that may end the amplicon before one starting with fwds[fi]
    preds: list[list[int]] = []
    for f in fwds:
        lo = bisect.bisect_left(r_starts, f.end + 1 + params.overlap)
        hi = bisect.bisect_right(r_starts, f.start + geom.max_outer_overlap - 1)
        preds.append([ri for ri in range(lo, hi)
                      if revs[ri].end - f.start + 1 <= geom.max_outer_overlap])

    first_ok = [f.start <= geom.end_slack for f in fwds]
    last_ok = [r.end >= genome_len - 1 - geom.end_slack for r in revs]

    def run(t: int) -> Path_ | None:
        # layer value: (sum of coverage, -bases before the first primer, -sum of penalty);
        # back[k][ri] = (fi, prev ri)
        value: dict[int, tuple[int, int, float]] = {}
        back: list[dict[int, tuple[int, int | None]]] = [{}]
        for ri, amps in enumerate(by_rev):
            for fi, cov, neg in amps:
                v = (cov, -fwds[fi].start, neg)
                if cov >= t and first_ok[fi] and (ri not in value or v > value[ri]):
                    value[ri] = v
                    back[0][ri] = (fi, None)
        for _ in range(n - 1):
            best_in: dict[int, tuple[tuple[int, int, float], int]] = {}
            for fi, rs in enumerate(preds):
                cands = [(value[ri], ri) for ri in rs if ri in value]
                if cands:
                    best_in[fi] = max(cands)
            new: dict[int, tuple[int, int, float]] = {}
            layer: dict[int, tuple[int, int | None]] = {}
            for ri, amps in enumerate(by_rev):
                for fi, cov, neg in amps:
                    if cov < t or fi not in best_in:
                        continue
                    (s, start, p), prev = best_in[fi]
                    v = (s + cov, start, p + neg)
                    if ri not in new or v > new[ri]:
                        new[ri] = v
                        layer[ri] = (fi, prev)
            value = new
            back.append(layer)
            if not value:
                return None
        ends = [((v[0], v[1] - (genome_len - 1 - revs[ri].end), v[2]), ri)
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

    levels = sorted({cov for amps in by_rev for _, cov, _ in amps})
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


def _dimers(amplicons: list[Amplicon]) -> list[str]:
    warnings: list[str] = []
    for pool in (1, 2):
        oligos: list[tuple[str, str]] = []
        for a in amplicons:
            if a.pool == pool:
                left, right = _primer_names(a.index)
                oligos += [(left, a.fwd.base_seq), (right, a.rev.base_seq)]
        for i, (name_a, seq_a) in enumerate(oligos):
            for name_b, seq_b in oligos[i + 1:]:
                dg = primer3.calc_heterodimer(seq_a, seq_b).dg
                if dg < DIMER_DG:
                    warnings.append(f"pool {pool}: {name_a} × {name_b} ΔG {dg / 1000:.1f} kcal/mol")
    return warnings


def _ref(cmap: CoordMap, cons: Consensus, pos: int, side: str) -> int | None:
    try:
        return cmap.col_to_ref(cons.cons_to_col[pos], side)  # type: ignore[arg-type]
    except PipelineError:
        return None


def _no_fit(
    fwds, revs, n: int, genome_len: int, params: TilingParams, cons: Consensus, cmap: CoordMap
) -> PipelineError:
    """Explain why no tiling exists: the longest primer-free stretch and N values that work."""
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
    path, geom = None, geoms[0]
    for geom in geoms:
        path = choose(fwds, revs, n, geom, params, genome_len)
        if path:
            break
    if not path:
        raise _no_fit(fwds, revs, n, genome_len, params, cons, cmap)
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
    return TilingResult(amplicons, geom, len(fwds), len(revs), _dimers(amplicons), ends)
