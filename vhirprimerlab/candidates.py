"""Single primers for tiling: enumerate on the consensus, keep conserved ones, score them."""

from __future__ import annotations

import dataclasses
from collections import defaultdict
from typing import Callable

import primer3
from Bio.Seq import reverse_complement

from vhirprimerlab.degenerate import degenerate_site, gc_percent, tm_range
from vhirprimerlab.models import Alignment, ColumnStats, Consensus, PrimerCandidate, TilingParams
from vhirprimerlab.scoring import match_primer
from vhirprimerlab.structure import structure_ok

Progress = Callable[[str, float], None]

# Genome ends are often AT-rich (the hMPV trailer, for one): a primer that reaches them
# may need more bases to reach the Tm, so there the length limit rises and GC % is free.
END_LEN_MAX = 32


def end_loss(c: PrimerCandidate, genome_len: int) -> int:
    """Genome bases left outside the insert if c is the first forward or last reverse primer."""
    return c.end + 1 if c.strand == "+" else genome_len - c.start


def _end_params(params: TilingParams) -> TilingParams:
    return dataclasses.replace(
        params, len_max=max(params.len_max, END_LEN_MAX), gc_min=0.0, gc_max=100.0
    )


def _rc_structure_ok(site: str) -> bool:
    return structure_ok(reverse_complement(site))


def _gc(seq: str) -> float:
    return 100.0 * sum(b in "GC" for b in seq) / len(seq)


def _best_length(
    cons_seq: str, anchor: int, strand: str, params: TilingParams
) -> PrimerCandidate | None:
    """Best primer whose 3' end sits at anchor (end for '+', start for '-')."""
    len_opt = (params.len_min + params.len_max) // 2
    best: PrimerCandidate | None = None
    for length in range(params.len_min, params.len_max + 1):
        start = anchor - length + 1 if strand == "+" else anchor
        if start < 0 or start + length > len(cons_seq):
            break
        site = cons_seq[start:start + length]
        if site.strip("ACGT"):
            break  # longer primers would contain the same non-ACGT base
        if not params.gc_min <= _gc(site) <= params.gc_max:
            continue
        seq = site if strand == "+" else reverse_complement(site)
        tm = primer3.calc_tm(seq)
        if not params.tm_min <= tm <= params.tm_max:
            continue
        penalty = abs(tm - params.tm_opt) + abs(length - len_opt)
        if best is None or penalty < best.penalty:
            best = PrimerCandidate(strand, start, length, site, seq, seq, tm, None, _gc(site), penalty)
    if best is None:
        return None
    return best if structure_ok(best.seq) else None


def enumerate_primers(
    cons_seq: str, params: TilingParams
) -> tuple[list[PrimerCandidate], list[PrimerCandidate]]:
    """Forward and reverse primers on the consensus, one per 3' anchor position.

    Within max_end_loss of either genome end the end rules apply (see END_LEN_MAX).
    """
    n, ends = len(cons_seq), _end_params(params)
    loss = params.max_end_loss
    fwds = [c for e in range(n)
            if (c := _best_length(cons_seq, e, "+", ends if e + 1 <= loss else params))]
    revs = [c for s in range(n)
            if (c := _best_length(cons_seq, s, "-", ends if n - s <= loss else params))]
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
) -> list[PrimerCandidate]:
    """Match every candidate against every sequence; add degenerate bases if allowed."""
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
            at_end = end_loss(c, len(cons.seq)) <= params.max_end_loss
            gc = (0.0, 100.0) if at_end else (params.gc_min, params.gc_max)
            site = degenerate_site(
                c.site, cols, aln, sp, params.max_degeneracy, side, gc,
                accept=structure_ok if c.strand == "+" else _rc_structure_ok,
            )
            if site != c.site:
                seq = site if c.strand == "+" else reverse_complement(site)
                lo, hi = tm_range(seq)
                scored = dataclasses.replace(
                    scored, site=site, seq=seq, tm=(lo + hi) / 2, tm_range=(lo, hi), gc=gc_percent(seq),
                    mask=_mask(site, cols, aln, side, sp),
                )
        out.append(scored)
    return out
