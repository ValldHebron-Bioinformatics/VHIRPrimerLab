"""Single primers for tiling: enumerate on the consensus, keep conserved ones, score them."""

from __future__ import annotations

import dataclasses
from collections import defaultdict
from typing import Callable

from Bio.Seq import reverse_complement

from vhirprimerlab.degenerate import degenerate_site, gc_percent, tm_range
from vhirprimerlab.models import Alignment, ColumnStats, Consensus, PrimerCandidate, TilingParams
from vhirprimerlab.scoring import match_primer
from vhirprimerlab.structure import structure_ok
from vhirprimerlab.thermo import Thermo

Progress = Callable[[str, float], None]

# Genome ends are often AT-rich (the hMPV trailer, for one): a primer that reaches them
# may need more bases to reach the Tm, so there the length limit rises and GC % is free.
END_LEN_MAX = 32


def end_loss(c: PrimerCandidate, span: tuple[int, int]) -> int:
    """Span bases left outside the insert if c is the first forward or last reverse primer.

    span is the (first, last) consensus position tiled; the whole genome is (0, length - 1).
    """
    return c.end + 1 - span[0] if c.strand == "+" else span[1] + 1 - c.start


def loss_limit(params: TilingParams, strand: str) -> int:
    """Bases that may stay outside the insert at the end this primer reaches."""
    return params.max_start_loss if strand == "+" else params.max_end_loss


def _end_params(params: TilingParams) -> TilingParams:
    return dataclasses.replace(
        params, len_max=max(params.len_max, END_LEN_MAX), gc_min=0.0, gc_max=100.0
    )


def _rc_structure_ok(site: str, thermo: Thermo) -> bool:
    return structure_ok(reverse_complement(site), thermo)


def _gc(seq: str) -> float:
    return 100.0 * sum(b in "GC" for b in seq) / len(seq)


def _best_length(
    cons_seq: str, anchor: int, strand: str, params: TilingParams, span: tuple[int, int]
) -> PrimerCandidate | None:
    """Best primer whose 3' end sits at anchor (end for '+', start for '-')."""
    len_opt = (params.len_min + params.len_max) // 2
    best: PrimerCandidate | None = None
    for length in range(params.len_min, params.len_max + 1):
        start = anchor - length + 1 if strand == "+" else anchor
        if start < span[0] or start + length > span[1] + 1:
            break
        site = cons_seq[start:start + length]
        if site.strip("ACGT"):
            break  # longer primers would contain the same non-ACGT base
        if not params.gc_min <= _gc(site) <= params.gc_max:
            continue
        seq = site if strand == "+" else reverse_complement(site)
        tm = params.thermo.tm(seq)
        if not params.tm_min <= tm <= params.tm_max:
            continue
        penalty = abs(tm - params.tm_opt) + abs(length - len_opt)
        if best is None or penalty < best.penalty:
            best = PrimerCandidate(strand, start, length, site, seq, seq, tm, None, _gc(site), penalty)
    if best is None:
        return None
    return best if structure_ok(best.seq, params.thermo) else None


def enumerate_primers(
    cons_seq: str, params: TilingParams, span: tuple[int, int] | None = None
) -> tuple[list[PrimerCandidate], list[PrimerCandidate]]:
    """Forward and reverse primers on the consensus, one per 3' anchor position.

    Primers lie inside span (default: the whole consensus). Within max_start_loss of the
    span start, and max_end_loss of its end, the end rules apply (see END_LEN_MAX).
    """
    span = span or (0, len(cons_seq) - 1)
    first, last = span
    ends = _end_params(params)

    def best(anchor: int, strand: str, near_end: bool) -> PrimerCandidate | None:
        # Near an end the ordinary rules come first and the end rules (up to 32 nt, free GC) only
        # fill the anchors the ordinary ones leave bare. A primer then never changes when
        # --max-start-loss or --max-end-loss grows to reach it. Choosing the end-rule primer
        # whenever one existed made it longer and so more variable, so pick_conserved dropped it
        # from its bin, and the advice "set the limit to N" was wrong by the one base that moved
        # the primer inside the allowance (and then again at N + 1, and so on).
        found = _best_length(cons_seq, anchor, strand, params, span)
        return found or (_best_length(cons_seq, anchor, strand, ends, span) if near_end else None)

    fwds = [c for e in range(first, last + 1)
            if (c := best(e, "+", e + 1 - first <= params.max_start_loss))]
    revs = [c for s in range(first, last + 1)
            if (c := best(s, "-", last + 1 - s <= params.max_end_loss))]
    return fwds, revs


def _anchor(c: PrimerCandidate) -> int:
    return c.end if c.strand == "+" else c.start


def _variability(c: PrimerCandidate, freq: list[float], window: int) -> float:
    """Sum of (1 - majority frequency) over the site, 3' window bases weighted 3x."""
    total = 0.0
    for k in range(c.length):
        in_window = k >= c.length - window if c.strand == "+" else k < window
        total += (3.0 if in_window else 1.0) * (1.0 - freq[c.start + k])
    return total


def pick_conserved(
    cands: list[PrimerCandidate],
    stats: list[ColumnStats],
    cons: Consensus,
    params: TilingParams,
    bin_size: int = 40,
    per_bin: int = 3,
) -> list[PrimerCandidate]:
    """The per_bin least variable candidates of every bin of 3'-anchor positions.

    Every bin keeps its best candidates however variable the region is: in a variable
    stretch (an attachment-protein gene, say) the best available primers can still match
    most sequences, and the tiling search, which maximises the weakest amplicon's
    coverage, decides whether to use them.
    """
    freq = [stats[col].majority_freq for col in cons.cons_to_col]
    bins: dict[int, list[tuple[float, PrimerCandidate]]] = defaultdict(list)
    for c in cands:
        bins[_anchor(c) // bin_size].append((_variability(c, freq, params.window), c))
    picked: list[PrimerCandidate] = []
    for b in sorted(bins):
        keep = sorted(bins[b], key=lambda vc: (vc[0], vc[1].penalty))[:per_bin]
        picked += sorted((c for _, c in keep), key=_anchor)
    return picked


def _mask(site: str, cols: list[int], aln: Alignment, side, params) -> int:
    lo, hi = cols[0], cols[-1] + 1
    offset = [c - lo for c in cols]
    memo: dict[str, bool] = {}
    mask = 0
    for i, row in enumerate(aln.rows):
        key = row[lo:hi]
        ok = memo.get(key)
        if ok is None:
            ok = memo[key] = match_primer(site, offset, key, side, params).ok
        if ok:
            mask |= 1 << i
    return mask


def score_primers(
    cands: list[PrimerCandidate],
    aln: Alignment,
    cons: Consensus,
    params: TilingParams,
    progress: Progress | None = None,
    span: tuple[int, int] | None = None,
) -> list[PrimerCandidate]:
    """Match every candidate against every sequence; add degenerate bases if allowed.

    span (default: the whole consensus) sets where GC % is free for degenerate primers.
    """
    span = span or (0, len(cons.seq) - 1)
    full = (1 << len(aln.rows)) - 1
    sp = params.scoring
    out: list[PrimerCandidate] = []
    for i, c in enumerate(cands):
        if progress and i % 50 == 0:
            progress("scoring", i / max(len(cands), 1))
        cols = cons.cons_to_col[c.start:c.start + c.length]
        side = "right" if c.strand == "+" else "left"
        mask = _mask(c.site, cols, aln, side, sp)
        scored = dataclasses.replace(c, cols=cols, mask=mask)
        if params.degenerate and mask != full:
            at_end = end_loss(c, span) <= loss_limit(params, c.strand)
            gc = (0.0, 100.0) if at_end else (params.gc_min, params.gc_max)
            site = degenerate_site(
                c.site, cols, aln, sp, params.max_degeneracy, side, gc,
                accept=(lambda s: structure_ok(s, params.thermo)) if c.strand == "+"
                else (lambda s: _rc_structure_ok(s, params.thermo)),
            )
            if site != c.site:
                seq = site if c.strand == "+" else reverse_complement(site)
                lo, hi = tm_range(seq, params.thermo)
                scored = dataclasses.replace(
                    scored, site=site, seq=seq, tm=(lo + hi) / 2, tm_range=(lo, hi), gc=gc_percent(seq),
                    mask=_mask(site, cols, aln, side, sp),
                )
        out.append(scored)
    return out
