"""Tiling the whole genome with N overlapping amplicons in two pools."""

from __future__ import annotations

import bisect
import dataclasses
import logging
import math
from typing import Callable

from vhirprimerlab.candidates import (
    CandidateMemo,
    end_loss,
    enumerate_primers,
    loss_limit,
    pick_conserved,
    score_primers,
    structure_valid,
)
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.models import (
    Alignment,
    Amplicon,
    ColumnStats,
    Consensus,
    EdgeAdvice,
    Geometry,
    PipelineError,
    PrimerCandidate,
    Region,
    TilingParams,
    TilingResult,
    primer_names,
)
from vhirprimerlab.scoring import score_pair
from vhirprimerlab.structure import dimer_dg
from vhirprimerlab.thermo import DEFAULT_THERMO, Thermo

log = logging.getLogger("vhirprimerlab")

Progress = Callable[[str, float], None]
Path_ = list[tuple[PrimerCandidate, PrimerCandidate]]

TOLERANCES = (0.20, 0.35, 0.50)  # amplicon length tolerances tried in order
LONG_AMPLICON = 3000  # warn above this target length
MAX_DIMER_SWAPS = 20  # primers dropped at most while trying to remove dimers


def target_length(genome_len: int, n: int, overlap: int) -> int:
    """Amplicon length that tiles genome_len with n amplicons overlapping by overlap."""
    return round((genome_len + (n - 1) * overlap) / n)


def amplicons_for(span_len: int, amplicon_len: int, overlap: int) -> int:
    """Fewest amplicons of amplicon_len, overlapping by overlap, that cover span_len."""
    if amplicon_len <= overlap:
        raise PipelineError("--amplicon-len must be larger than --overlap")
    return max(1, math.ceil((span_len - overlap) / (amplicon_len - overlap)))


def geometry(
    genome_len: int, n: int, overlap: int, tolerance: float, primer_len_max: int
) -> Geometry | None:
    """Size limits for n amplicons, or None when amplicons would be too short.

    Neighbours may overlap (outer primer edges) by at most half the shortest
    amplicon, which keeps amplicons k and k+2 (same pool) apart. That overlap must
    still hold two primers plus the required insert overlap. A single amplicon has
    no neighbour, so it only needs to hold two primers and at least one base of insert.
    """
    target = target_length(genome_len, n, overlap)
    len_min = math.floor(target * (1 - tolerance))
    max_outer = (len_min - 1) // 2
    if n > 1 and max_outer < overlap + 2 * primer_len_max:
        return None
    if n == 1 and len_min < 2 * primer_len_max + 1:
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


def min_amplicon_len(span_len: int, overlap: int, primer_len_max: int) -> int:
    """Shortest --amplicon-len that still tiles span_len, or 0 when none does.

    Tiling a span derives N from --amplicon-len, so this is the number to name when the
    amplicons come out too short: fewer amplicons are longer, so the most the geometry
    allows is the fewest the span can be cut into, and this length asks for that many.
    """
    most = max_amplicons(span_len, overlap, primer_len_max)
    return overlap + max(1, math.ceil((span_len - overlap) / most)) if most else 0


def _popcount(x: int) -> int:
    return bin(x).count("1")


def choose(
    fwds: list[PrimerCandidate],
    revs: list[PrimerCandidate],
    n: int,
    geom: Geometry,
    params: TilingParams,
    span: tuple[int, int],
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

    # the forward primers each reverse primer may follow, the other way round
    succ: list[list[int]] = [[] for _ in revs]
    for fi, rs in enumerate(preds):
        for ri in rs:
            succ[ri].append(fi)

    first_ok = [end_loss(f, span) <= params.max_start_loss for f in fwds]
    last_ok = [end_loss(r, span) <= params.max_end_loss for r in revs]

    def run(t: int) -> Path_ | None:
        # Only amplicons covering at least t sequences can be used, whatever the layer, so they
        # are picked out once. A layer's value is (sum of coverage, sum of primer matches, -bases
        # before the first primer, -sum of penalty); back[k][ri] = (fi, prev ri).
        viable = [[a for a in amps if a[1] >= t] for amps in by_rev]
        value: dict[int, tuple[int, int, int, float]] = {}
        back: list[dict[int, tuple[int, int | None]]] = [{}]
        for ri, amps in enumerate(viable):
            for fi, cov, pm, neg in amps:
                v = (cov, pm, -fwds[fi].start, neg)
                if first_ok[fi] and (ri not in value or v > value[ri]):
                    value[ri] = v
                    back[0][ri] = (fi, None)
        for _ in range(n - 1):
            # the best way to arrive at each forward primer: the largest (value, reverse primer)
            best_in: dict[int, tuple[tuple[int, int, int, float], int]] = {}
            for ri, v in value.items():
                for fi in succ[ri]:
                    cand = (v, ri)
                    if fi not in best_in or cand > best_in[fi]:
                        best_in[fi] = cand
            new: dict[int, tuple[int, int, int, float]] = {}
            layer: dict[int, tuple[int, int | None]] = {}
            for ri, amps in enumerate(viable):
                for fi, cov, pm, neg in amps:
                    if fi not in best_in:
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


def _pair_name(index: int) -> str:
    """Default name of the amplicon at a 1-based whole-genome index."""
    return f"VHIRPrimerLab_{index}"


Key = tuple[str, int, int]  # (strand, start, length) identifies a candidate primer


def _key(c: PrimerCandidate) -> Key:
    return c.strand, c.start, c.length


def _dimer_pairs(path: Path_, names: list[str], pools: list[int] | None = None,
                 thermo: Thermo = DEFAULT_THERMO) -> list[
        tuple[int, str, str, float, PrimerCandidate, PrimerCandidate]]:
    """Primer pairs in one pool whose worst variant pair has ΔG below thermo.max_dimer_dg.

    Amplicon k (1-based) is named names[k - 1] and is in pools[k - 1]; without pools,
    pool 1 when k is odd and pool 2 when even. Every variant of a degenerate primer is
    checked against every variant of the other primer.
    """
    pools = pools or [1 if k % 2 else 2 for k in range(1, len(path) + 1)]
    found = []
    for pool in (1, 2):
        oligos: list[tuple[str, PrimerCandidate]] = []
        for (f, r), name, amp_pool in zip(path, names, pools):
            if amp_pool == pool:
                left, right = primer_names(name)
                oligos += [(left, f), (right, r)]
        for i, (name_a, a) in enumerate(oligos):
            for name_b, b in oligos[i + 1:]:
                dg = dimer_dg(a.seq, b.seq, thermo)
                if dg < thermo.max_dimer_dg:
                    found.append((pool, name_a, name_b, dg, a, b))
    return found


def _dimers(path: Path_, names: list[str], pools: list[int] | None = None,
            thermo: Thermo = DEFAULT_THERMO) -> list[str]:
    return [f"pool {pool}: {na} × {nb} ΔG {dg / 1000:.1f} kcal/mol"
            for pool, na, nb, dg, _, _ in _dimer_pairs(path, names, pools, thermo)]


def assign_pools(spans: list[list[tuple[str, int, int]]]) -> tuple[list[int], list[str]]:
    """Pools (1 or 2) of the amplicons of several spans, and a warning for each clash left.

    Amplicons are numbered across the scheme and alternate pools by that number, which keeps
    neighbours of a span apart. Spans that overlap one another (neighbouring genes, say) can
    put overlapping amplicons in one tube, where each one's primers prime the other's
    product. A whole span is therefore flipped to the other parity when that overlaps
    fewer same-pool amplicons; a clash no choice avoids (three amplicons that overlap
    one another cannot be split over two pools) is reported. Each span is (name, start, end).
    """
    placed: list[tuple[str, int, int, int]] = []  # name, start, end, pool

    def clashes(candidate: list[tuple[str, int, int, int]]) -> list[tuple[str, str, int]]:
        found = []
        for name, start, end, pool in candidate:
            for other, o_start, o_end, o_pool in placed:
                shared = min(end, o_end) - max(start, o_start) + 1
                if o_pool == pool and shared > 0:
                    found.append((other, name, shared))
        return found

    warnings: list[str] = []
    for span in spans:
        options = []
        for flip in (0, 1):
            candidate = [(name, start, end, 1 + ((len(placed) + j + flip) % 2))
                         for j, (name, start, end) in enumerate(span)]
            found = clashes(candidate)
            options.append(((len(found), sum(c[2] for c in found)), candidate, found))
        _, candidate, found = min(options, key=lambda o: o[0])  # the numbering's own pools win ties
        warnings += [f"{a} and {b} are in the same pool and their amplicons overlap by {n} bp, so "
                     f"the primers of one can prime the other's product; tile overlapping "
                     f"regions in separate runs" for a, b, n in found]
        placed += candidate
    return [pool for *_, pool in placed], warnings


def _worst_coverage(path: Path_) -> int:
    return min(_popcount(f.mask & r.mask) for f, r in path)


def _choose_any(fwds, revs, n, geoms, params, span, banned: set[Key]):
    """choose() at the first length tolerance that works, skipping banned primers."""
    fwds = [c for c in fwds if _key(c) not in banned]
    revs = [c for c in revs if _key(c) not in banned]
    for geom in geoms:
        path = choose(fwds, revs, n, geom, params, span)
        if path:
            return path, geom
    return None, geoms[-1]


def _avoid_dimers(path: Path_, names: list[str], geom: Geometry, fwds, revs, n, geoms, params, span):
    """Swap out primers that form dimers within a pool, when that costs no coverage.

    For the first dimer pair not yet given up on, the tiling is chosen again without
    one of its primers (the one in most dimer pairs first). The new tiling is kept when
    it has fewer dimer pairs and its weakest amplicon covers at least as many
    sequences as before; if neither primer can go, the pair is kept and reported.
    """
    banned: set[Key] = set()
    kept: set[tuple[Key, Key]] = set()
    pairs = _dimer_pairs(path, names, thermo=params.thermo)
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
            trial, trial_geom = _choose_any(fwds, revs, n, geoms, params, span,
                                            banned | {_key(c)})
            if trial is None or _worst_coverage(trial) < _worst_coverage(path):
                continue
            trial_pairs = _dimer_pairs(trial, names, thermo=params.thermo)
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
    fwds, revs, n: int, span: tuple[int, int], params: TilingParams, cons: Consensus, cmap: CoordMap
) -> PipelineError:
    """Explain why no tiling exists: a genome end out of reach, or the longest primer-free
    stretch and N values that work."""
    genome_len = span[1] - span[0] + 1
    what = "genome" if span == (0, len(cons.seq) - 1) else "span"
    for side, strand, cands in (("start", "+", fwds), ("end", "-", revs)):
        limit = loss_limit(params, strand)
        closest = min((end_loss(c, span) for c in cands), default=genome_len)
        if closest > limit:
            return PipelineError(
                f"No tiling with {n} amplicons: no acceptable primer within {limit} nt of the "
                f"{what} {side} (the closest leaves {closest:,} nt out); "
                f"set --max-{side}-loss to {closest:,} or more"
            )
    sites = sorted((c.start, c.end) for c in fwds + revs)
    gap, prev_end = (0, span[0] - 1, span[0] - 1), span[0] - 1
    for start, end in sites + [(span[1] + 1, span[1] + 1)]:
        if start - prev_end - 1 > gap[0]:
            gap = (start - prev_end - 1, prev_end + 1, start - 1)
        prev_end = max(prev_end, end)
    a = _ref(cmap, cons, gap[1], "start") or gap[1] + 1
    b = _ref(cmap, cons, gap[2], "end") or gap[2] + 1
    limit = max_amplicons(genome_len, params.overlap, params.len_max)
    # a primer-free stretch is bridged by fewer, longer amplicons: look below n first
    tries = [n - 1, n - 2, n - 3, n // 2, n // 3, n // 4, n + 1, n + 2, n + 3]
    works: list[int] = []
    for m in dict.fromkeys(m for m in tries if 1 <= m <= limit and m != n):
        g = geometry(genome_len, m, params.overlap, TOLERANCES[-1], params.len_max)
        if g and choose(fwds, revs, m, g, params, span):
            works.append(m)
            if len(works) == 3:
                break
    if what == "span":
        # Over a span N is derived from --amplicon-len, so naming N would name a number the
        # caller never chose. Each N that works is offered as the shortest amplicon length that
        # asks for exactly that many, and only when it really does.
        lens = [str(length) for m in works
                if (length := params.overlap + math.ceil((genome_len - params.overlap) / m))
                > params.overlap and amplicons_for(genome_len, length, params.overlap) == m]
        advice = f"try --amplicon-len {', '.join(lens) if lens else 'a different value'}"
    else:
        advice = f"try -n {', '.join(str(m) for m in works) if works else 'a different value'}"
    return PipelineError(
        f"No tiling with {n} amplicons: no acceptable primers between {a:,} and {b:,}; {advice}"
    )


def _candidates(aln: Alignment, cons: Consensus, stats: list[ColumnStats], params: TilingParams,
                span: tuple[int, int], progress: Progress | None = None,
                memo: CandidateMemo | None = None):
    """Forward and reverse primers for span, scored against every sequence.

    With a memo (made for these parameters, span and alignment) primers found and scored before
    are reused: a re-run that only moves the end allowances repeats the work for the ends alone.
    """
    if memo is not None and not memo.fits(params, span, aln):
        memo = None
    # Most primers are never kept, and checking a primer's hairpins and dimers is what takes the
    # time: that check is left until pick_conserved is about to keep one (the same primers result).
    raw_f, raw_r = enumerate_primers(cons.seq, params, span, memo, check_structure=False)
    valid = structure_valid(params)
    picked = (pick_conserved(raw_f, stats, cons, params, valid=valid)
              + pick_conserved(raw_r, stats, cons, params, valid=valid))
    scored = score_primers(picked, aln, cons, params, progress=progress, span=span, memo=memo)
    return [c for c in scored if c.strand == "+"], [c for c in scored if c.strand == "-"]


def _coverages(path: Path_, aln: Alignment, params: TilingParams) -> list[float]:
    """Each amplicon's coverage %, scored exactly as Amplicon.coverage is."""
    return [score_pair(f.site, f.cols, r.site, r.cols, aln, params.scoring).coverage_pct
            for f, r in path]


def _edge_advice(amplicons: list[Amplicon], fwds, revs, n, geoms, geom: Geometry,
                 params: TilingParams, span, aln: Alignment, cons: Consensus,
                 stats: list[ColumnStats], memo: CandidateMemo | None = None) -> list[EdgeAdvice]:
    """Where giving up more bases at one span edge would lift that edge amplicon's coverage.

    Judged on the edge amplicon itself; the weakest amplicon is only a guard, so advice
    never makes it worse. Each end is explored alone, with the other end's limit
    untouched, because a generous allowance at both ends also spends bases at the end that
    gains nothing. Exploring reuses the scored primers, which is cheap, but the primers
    depend on the limit (the end rules reach as far as it does), so every answer is
    checked by choosing again, dimers avoided, under the limit it names, and only what
    that check shows is advised. Nothing here changes the tiling chosen.
    """
    worst = min(a.coverage.coverage_pct for a in amplicons)
    names = [_pair_name(k) for k in range(1, n + 1)]
    advice: list[EdgeAdvice] = []
    for end, edge, flag, at in (("start", amplicons[0], "max_start_loss", 0),
                                ("end", amplicons[-1], "max_end_loss", -1)):
        asked = getattr(params, flag)
        if not edge.below_threshold or asked >= geom.target_len:
            continue  # one amplicon is as far in as it is worth looking
        before = edge.coverage.coverage_pct
        cur, best, found = params, before, None
        cur_f, cur_r = fwds, revs
        while getattr(cur, flag) < geom.target_len:
            path, _ = _choose_any(cur_f, cur_r, n, geoms, dataclasses.replace(
                cur, **{flag: geom.target_len}), span, set())
            if not path:
                break
            pcts = _coverages(path, aln, params)
            if pcts[at] <= best or min(pcts) < worst:
                break  # nothing better to prove, so no costly check
            # the better tiling's own loss, not the cap explored: re-running with it keeps this
            # end from spending more than the tiling needs, and the other end is as it was
            limit = end_loss(path[0][0] if end == "start" else path[-1][1], span)
            if limit <= getattr(cur, flag):
                break
            cur = dataclasses.replace(params, **{flag: limit})
            cur_f, cur_r = _candidates(aln, cons, stats, cur, span, memo=memo)
            path, g = _choose_any(cur_f, cur_r, n, geoms, cur, span, set())
            if not path:
                break
            # a re-run avoids dimers, which can swap primers: promise what it produces
            path, _ = _avoid_dimers(path, names, g, cur_f, cur_r, n, geoms, cur, span)
            pcts = _coverages(path, aln, params)
            if pcts[at] <= best or min(pcts) < worst:
                break
            best, found = pcts[at], limit
        if found is not None:
            advice.append(EdgeAdvice(None, end, found, found - asked, before, best))
    return advice


def tile(
    aln: Alignment,
    cons: Consensus,
    stats: list[ColumnStats],
    cmap: CoordMap,
    params: TilingParams,
    progress: Progress | None = None,
    span: tuple[int, int] | None = None,
) -> TilingResult:
    """Choose params.n_amplicons overlapping amplicons covering span (default: the whole consensus)."""
    report = progress or (lambda stage, fraction: None)
    span = span or (0, len(cons.seq) - 1)
    n, genome_len = params.n_amplicons, span[1] - span[0] + 1
    geoms = [g for tol in TOLERANCES
             if (g := geometry(genome_len, n, params.overlap, tol, params.len_max))]
    if not geoms:
        target = target_length(genome_len, n, params.overlap)
        if span == (0, len(cons.seq) - 1):
            raise PipelineError(
                f"N too large: amplicons would be ~{target} bp; "
                f"use at most {max_amplicons(genome_len, params.overlap, params.len_max)}"
            )
        # Over a span N came from --amplicon-len, so that is the flag to name.
        least = min_amplicon_len(genome_len, params.overlap, params.len_max)
        if not least:
            raise PipelineError(
                f"this {genome_len:,} bp span is too short for one amplicon with primers of up "
                f"to {params.len_max} bases"
            )
        raise PipelineError(
            f"--amplicon-len {params.amplicon_len} gives {n} amplicons of ~{target} bp here, too "
            f"short to hold two primers and a {params.overlap} bp overlap; use "
            f"--amplicon-len {least} or more"
        )
    if geoms[0].target_len > LONG_AMPLICON:
        log.warning("Amplicons will be ~%d bp; long PCRs are harder (try a larger -n)",
                    geoms[0].target_len)

    report("finding primers", 0.10)
    memo = CandidateMemo(params, span, aln)  # lets the edge advice reuse the primers away from the ends
    fwds, revs = _candidates(aln, cons, stats, params, span,
                             lambda stage, f: report(stage, 0.15 + 0.65 * f), memo)
    log.info("Candidate primers: %d forward, %d reverse", len(fwds), len(revs))

    report("choosing amplicons", 0.85)
    path, geom = _choose_any(fwds, revs, n, geoms, params, span, set())
    if not path:
        raise _no_fit(fwds, revs, n, span, params, cons, cmap)
    names = [_pair_name(k) for k in range(1, len(path) + 1)]
    report("checking primer dimers", 0.90)
    path, geom = _avoid_dimers(path, names, geom, fwds, revs, n, geoms, params, span)
    if geom.tolerance != TOLERANCES[0]:
        log.warning("Used ±%d %% amplicon length tolerance (no tiling within ±%d %%)",
                    round(geom.tolerance * 100), round(TOLERANCES[0] * 100))

    amplicons: list[Amplicon] = []
    for k, (f, r) in enumerate(path, start=1):
        cov = score_pair(f.site, f.cols, r.site, r.cols, aln, params.scoring)
        nxt = path[k][0] if k < len(path) else None
        amplicons.append(Amplicon(
            index=k, name=names[k - 1], pool=1 if k % 2 else 2, fwd=f, rev=r, coverage=cov,
            below_threshold=cov.coverage_pct < params.min_coverage,
            overlap_next=None if nxt is None else r.start - nxt.end - 1,
        ))
    first = _ref(cmap, cons, path[0][0].end + 1, "start")
    last = _ref(cmap, cons, path[-1][1].start - 1, "end")
    # the genome's own ends bound a span that reaches them
    span_first = 1 if span[0] == 0 else _ref(cmap, cons, span[0], "start")
    span_last = cmap.ref_length if span[1] == len(cons.seq) - 1 else _ref(cmap, cons, span[1], "end")
    ends = (first - span_first if first and span_first else 0,
            span_last - last if last and span_last else 0)
    advice = _edge_advice(amplicons, fwds, revs, n, geoms, geom, params, span, aln, cons, stats, memo)
    return TilingResult(amplicons, geom, len(fwds), len(revs), _dimers(path, names, thermo=params.thermo), ends,
                        edge_advice=advice)


def tile_spans(
    aln: Alignment,
    cons: Consensus,
    stats: list[ColumnStats],
    cmap: CoordMap,
    params: TilingParams,
    spans: list[Region],
    progress: Progress | None = None,
) -> TilingResult:
    """Tile each span on its own and merge the results into one multiplex scheme.

    Amplicons are numbered across the scheme, pools alternate over that numbering, and
    the dimer check runs once over the merged scheme because every primer of a pool
    shares a tube whatever span it came from.
    """
    if not spans:
        raise PipelineError("No regions to tile")
    report = progress or (lambda stage, fraction: None)
    parts: list[TilingResult] = []
    amplicons: list[Amplicon] = []
    path: Path_ = []
    names: list[str] = []
    extents: list[list[tuple[str, int, int]]] = []
    for i, span in enumerate(spans):
        lo, length = cmap.region_to_included(span)
        n = amplicons_for(length, params.amplicon_len, params.overlap)
        try:
            part = tile(aln, cons, stats, cmap, dataclasses.replace(params, n_amplicons=n),
                        progress=lambda stage, f, i=i: report(stage, (i + f) / len(spans)),
                        span=(lo, lo + length - 1))
        except PipelineError as e:
            raise PipelineError(f"Region '{span.name}': {e}") from e
        parts.append(part)
        extents.append([])
        for k, a in enumerate(part.amplicons, start=1):
            name = span.name if n == 1 else f"{span.name}_{k}"
            index = len(amplicons) + 1
            last = k == len(part.amplicons)
            amplicons.append(dataclasses.replace(
                a, index=index, name=name, overlap_next=None if last else a.overlap_next))
            extents[-1].append((name, a.fwd.start, a.rev.end))
            path.append((a.fwd, a.rev))
            names.append(name)
    # A span of several amplicons names them NAME_1, NAME_2, ...; a span of one keeps its own
    # name. A span called F_1 beside a span called F therefore collides, and primers.tsv and
    # scheme.bed would carry two different oligos called F_1_LEFT. The browser can produce it: a
    # gene F taking several amplicons, and a sibling feature with no /gene whose product
    # sanitises to F_1.
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise PipelineError(
            f"Two amplicons would both be called {', '.join(repeated)}, so their primers could "
            "not be told apart in primers.tsv or scheme.bed: rename one of the regions"
        )
    pools, pool_warnings = assign_pools(extents)
    amplicons = [dataclasses.replace(a, pool=pool) for a, pool in zip(amplicons, pools)]
    for warning in pool_warnings:
        log.warning(warning)
    ends = (max(p.uncovered_ends[0] for p in parts), max(p.uncovered_ends[1] for p in parts))
    advice = [dataclasses.replace(a, span=span.name)
              for span, part in zip(spans, parts) for a in part.edge_advice]
    return TilingResult(amplicons, parts[0].geometry,
                        sum(p.n_fwd for p in parts), sum(p.n_rev for p in parts),
                        _dimers(path, names, pools, params.thermo), ends, spanned=True, edge_advice=advice,
                        pool_warnings=pool_warnings)
