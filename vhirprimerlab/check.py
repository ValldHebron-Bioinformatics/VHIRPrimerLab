"""Checking primer pairs the user already has against every aligned sequence."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path

import primer3
from Bio.Seq import reverse_complement

from vhirprimerlab.degenerate import fold, gc_percent, tm_range
from vhirprimerlab.models import (
    Alignment,
    Assay,
    CheckResult,
    Consensus,
    Location,
    PairCoverage,
    PipelineError,
    ScoringParams,
)
from vhirprimerlab.regions import NAME_RE
from vhirprimerlab.scoring import IUPAC, score_pair

log = logging.getLogger("vhirprimerlab")

TM_FOLD_LIMIT = 256  # above this, enumerating every expansion for a Tm range is too slow
NOT_DOWNSTREAM = "reverse site is not downstream of forward site"


def _clean(field: str) -> str:
    return field.strip().upper().replace("U", "T")


def _looks_like_header(fwd: str) -> bool:
    """A column title like 'forward' is mostly not A/C/G/T; a primer, even with a typo, is."""
    return sum(ch in "ACGT" for ch in fwd) * 2 < len(fwd)


def _check_oligo(seq: str, what: str, where: str, window: int) -> None:
    for pos, ch in enumerate(seq, start=1):
        if ch not in IUPAC:
            raise PipelineError(f"{where}: {what} primer has an invalid character at position {pos}: '{ch}'")
    if len(seq) < window:
        raise PipelineError(
            f"{where}: {what} primer ({len(seq)} bases) is shorter than --three-prime-window ({window})"
        )


def load_assays(path: Path, window: int) -> list[Assay]:
    """Read a tab-separated name/forward/reverse file (optional header, # comments)."""
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"Primer file not found: {path}")
    assays: list[Assay] = []
    seen_data = False
    for lineno, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = [f.strip() for f in line.split("\t")]
        while fields and not fields[-1]:
            fields.pop()
        where = f"{path.name} line {lineno}"
        if len(fields) != 3:
            raise PipelineError(f"{where}: expected 3 tab-separated columns (name, forward, reverse)")
        name, fwd, rev = fields[0], _clean(fields[1]), _clean(fields[2])
        if not seen_data:
            seen_data = True
            if _looks_like_header(fwd):
                continue  # header line
        if not NAME_RE.fullmatch(name):
            raise PipelineError(f"{where}: invalid name '{name}' (use letters, digits, _ . -)")
        _check_oligo(fwd, "forward", where, window)
        _check_oligo(rev, "reverse", where, window)
        assays.append(Assay(name, fwd, rev))
    if not assays:
        raise PipelineError(f"Primer file contains no primer pairs: {path}")
    dups = sorted(n for n, c in Counter(a.name for a in assays).items() if c > 1)
    if dups:
        raise PipelineError(f"{path.name}: duplicate pair names: {', '.join(dups)}")
    return assays


def _best_sites(site: str, consensus_seq: str) -> tuple[int, list[int]]:
    """Fewest mismatches of site against the consensus, and every start that reaches it."""
    n = len(site)
    allowed = [IUPAC[b] for b in site]
    best, best_starts = n + 1, []
    for start in range(len(consensus_seq) - n + 1):
        window = consensus_seq[start:start + n]
        mm = sum(c not in ok for c, ok in zip(window, allowed))
        if mm < best:
            best, best_starts = mm, [start]
        elif mm == best:
            best_starts.append(start)
    return best, best_starts


def locate(site: str, consensus_seq: str) -> Location | None:
    """Start on the consensus with the fewest mismatches; None if tied or the site is too long."""
    if len(site) > len(consensus_seq):
        return None
    best, starts = _best_sites(site, consensus_seq)
    return Location(starts[0], best) if len(starts) == 1 else None


def _tm(primer: str, label: str, notes: list[str]) -> tuple[float | None, tuple[float, float] | None]:
    n = fold(primer)
    if n == 1:
        return primer3.calc_tm(primer), None
    if n > TM_FOLD_LIMIT:
        notes.append(f"Tm not computed: {label} fold {n} > {TM_FOLD_LIMIT}")
        return None, None
    lo, hi = tm_range(primer)
    return (lo + hi) / 2, (lo, hi)


def check_assay(
    assay: Assay, aln: Alignment, cons: Consensus, params: ScoringParams, min_coverage: float
) -> CheckResult:
    """Locate both primers on the consensus and score the pair on every sequence."""
    notes: list[str] = []
    fwd_site, rev_site = assay.forward, reverse_complement(assay.reverse)
    locs: dict[str, Location | None] = {}
    for label, site in (("fwd", fwd_site), ("rev", rev_site)):
        if len(site) > len(cons.seq):
            notes.append(f"{label} primer is longer than the consensus")
            locs[label] = None
            continue
        best, starts = _best_sites(site, cons.seq)
        loc = Location(starts[0], best) if len(starts) == 1 else None
        if best > params.max_mismatches:
            # a unique poor site is still scored (it usually covers nothing); ties are not located
            found = [f"no binding site for {label} primer within {params.max_mismatches} "
                     f"mismatches (best {best})"]
            flipped, _ = _best_sites(reverse_complement(site), cons.seq)
            if flipped <= params.max_mismatches:
                found.append(f"{label} primer may be written in the wrong orientation")
            notes += found
            log.warning("Pair %s: %s", assay.name, "; ".join(found))
        elif loc is None:
            notes.append(f"ambiguous binding site for {label} primer")
        locs[label] = loc
    fwd, rev = locs["fwd"], locs["rev"]
    located = fwd is not None and rev is not None
    if located and rev.start <= fwd.start + len(fwd_site) - 1:
        notes.append(NOT_DOWNSTREAM)
        located = False

    fwd_cols = rev_cols = []
    if located:
        fwd_cols = cons.cons_to_col[fwd.start:fwd.start + len(fwd_site)]
        rev_cols = cons.cons_to_col[rev.start:rev.start + len(rev_site)]
        coverage = score_pair(fwd_site, fwd_cols, rev_site, rev_cols, aln, params)
    else:
        coverage = PairCoverage([], {sid: "not located" for sid in aln.ids}, 0.0)

    fwd_tm, fwd_range = _tm(assay.forward, "fwd", notes)
    rev_tm, rev_range = _tm(assay.reverse, "rev", notes)
    return CheckResult(
        assay=assay, located=located, note="; ".join(notes), fwd=fwd, rev=rev,
        fwd_cols=list(fwd_cols), rev_cols=list(rev_cols), coverage=coverage,
        fwd_tm=fwd_tm, rev_tm=rev_tm, fwd_tm_range=fwd_range, rev_tm_range=rev_range,
        fwd_gc=gc_percent(assay.forward), rev_gc=gc_percent(assay.reverse),
        degenerate=fold(assay.forward) > 1 or fold(assay.reverse) > 1,
        below_threshold=not located or coverage.coverage_pct < min_coverage,
    )
