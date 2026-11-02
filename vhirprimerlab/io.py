"""File input/output: FASTA loading/validation and output writers."""

from __future__ import annotations

import csv
import logging
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

from Bio import SeqIO

from vhirprimerlab.degenerate import degenerate_bases, fold, gc_percent
from vhirprimerlab.models import (
    Alignment,
    CheckReport,
    CheckResult,
    PipelineError,
    Record,
    RegionResult,
    RunReport,
    TilingReport,
    TilingResult,
)
from vhirprimerlab.ranking import count_reasons

PER_PRIMER_DEGENERACY = [
    "forward_degenerate_bases", "reverse_degenerate_bases", "forward_variants", "reverse_variants",
]


def _degeneracy(fwd: str, rev: str) -> list[int]:
    """Values for PER_PRIMER_DEGENERACY."""
    return [degenerate_bases(fwd), degenerate_bases(rev), fold(fwd), fold(rev)]

if TYPE_CHECKING:
    from vhirprimerlab.coords import CoordMap

log = logging.getLogger(__name__)

IUPAC_NT = set("ACGTRYSWKMBDHVN")


def load_fasta(path: Path, allow_gaps: bool = False) -> list[Record]:
    """Load a multi-FASTA, upper-case it, convert U->T and validate it."""
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"Input file not found: {path}")
    text = path.read_text().replace("\r", "")
    if not text.strip():
        raise PipelineError(f"Input file is empty: {path}")
    if not text.lstrip().startswith(">"):
        raise PipelineError(f"Input file is not valid FASTA (no '>' header first): {path}")

    allowed = IUPAC_NT | {"-"} if allow_gaps else IUPAC_NT
    records: list[Record] = []
    with path.open() as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            seq = "".join(str(rec.seq).split()).upper().replace("U", "T")
            if not seq:
                raise PipelineError(f"Sequence '{rec.id}' is empty")
            for i, ch in enumerate(seq, start=1):
                if ch not in allowed:
                    raise PipelineError(
                        f"Sequence '{rec.id}' has an invalid character at position {i}: '{ch}'"
                    )
            records.append(Record(rec.id, seq))

    if not records:
        raise PipelineError(f"Input file is empty: {path}")
    dupes = sorted(rid for rid, n in Counter(r.id for r in records).items() if n > 1)
    if dupes:
        raise PipelineError(f"Duplicate IDs: {', '.join(dupes)}")
    log.info("Loaded %d sequences from %s", len(records), path)
    return records


def alignment_from_records(records: list[Record]) -> Alignment:
    """Build an Alignment from equal-length gapped records, dropping all-gap columns."""
    if not records:
        raise PipelineError("Alignment is empty")
    width = len(records[0].seq)
    for r in records:
        if len(r.seq) != width:
            raise PipelineError(
                f"Aligned sequences differ in length: '{r.id}' has {len(r.seq)} columns, "
                f"'{records[0].id}' has {width}"
            )
    keep = [c for c in range(width) if any(r.seq[c] != "-" for r in records)]
    rows = ["".join(r.seq[c] for c in keep) for r in records]
    return Alignment(ids=[r.id for r in records], rows=rows)


def load_alignment(path: Path) -> Alignment:
    """Load and validate a pre-aligned FASTA (gap character '-')."""
    return alignment_from_records(load_fasta(path, allow_gaps=True))


def write_alignment(aln: Alignment, path: Path) -> None:
    Path(path).write_text("".join(f">{i}\n{row}\n" for i, row in zip(aln.ids, aln.rows)))


TSV_COLUMNS = [
    "region", "rank", "forward_seq", "reverse_seq", "forward_tm", "reverse_tm", "forward_gc",
    "reverse_gc", "forward_len", "reverse_len", "amplicon_start", "amplicon_end", "amplicon_len",
    "coverage_pct", "n_covered", "n_total", "degenerate", *PER_PRIMER_DEGENERACY,
    "forward_tm_range", "reverse_tm_range", "penalty", "below_threshold", "not_covered",
]


def _fmt_range(r: tuple[float, float] | None) -> str:
    return f"{r[0]:.1f}-{r[1]:.1f}" if r else ""


def write_tsv(results: list[RegionResult], cmap: CoordMap, path: Path) -> None:
    """Write ranked pairs of every region; rank restarts at 1 per region."""
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_COLUMNS)
        for res in results:
            for rank, sp in enumerate(res.rank.pairs, start=1):
                start = cmap.col_to_ref(sp.fwd_cols[0], "start")
                end = cmap.col_to_ref(sp.rev_cols[-1], "end")
                cov = sp.coverage
                writer.writerow([
                    res.region.name, rank, sp.fwd_seq, sp.rev_seq,
                    f"{sp.pair.forward.tm:.1f}", f"{sp.pair.reverse.tm:.1f}",
                    f"{gc_percent(sp.fwd_seq):.1f}", f"{gc_percent(sp.rev_seq):.1f}",
                    sp.pair.forward.length, sp.pair.reverse.length,
                    start, end, end - start + 1,
                    f"{cov.coverage_pct:.1f}", len(cov.covered), len(cov.covered) + len(cov.not_covered),
                    sp.degenerate, *_degeneracy(sp.fwd_seq, sp.rev_seq), _fmt_range(sp.fwd_tm_range), _fmt_range(sp.rev_tm_range),
                    f"{sp.pair.penalty:.3f}", res.rank.below_threshold,
                    ",".join(f"{sid}({why})" for sid, why in cov.not_covered.items()),
                ])


def _region_section(res: RegionResult, min_coverage: float) -> list[str]:
    r = res.region
    title = f"Region {r.name} ({r.start}-{r.end})"
    start, length = res.included
    lines = [title, "-" * len(title),
             f"Consensus positions (primer3 included region): {start + 1}-{start + length}",
             f"primer3 candidate pairs: {res.n_candidates}",
             f"Pairs meeting --min-coverage ({min_coverage:g}%): {res.n_passing}"]
    if res.n_candidates == 0:
        lines.append("primer3 found no candidate pairs in this region. primer3 explanation:")
        lines += [f"  {k}: {v}" for k, v in res.explain.items()]
        return lines
    best = res.rank.pairs[0].coverage.coverage_pct
    lines.append(f"Best coverage: {best:.1f}%")
    if res.rank.below_threshold:
        lines.append(
            f"WARNING: no pair reached {min_coverage:g}% coverage; the best available pairs are "
            "reported with below_threshold=True."
        )
        if res.failures:
            lines.append("Why the best pair misses sequences (count of uncovered sequences):")
            lines += [f"  {reason}: {n}" for reason, n in res.failures]
    return lines


def _header_lines(report: RunReport | CheckReport | TilingReport, title: str) -> list[str]:
    """Title, command, parameters, input and alignment blocks shared by both reports."""
    lf = report.length
    out = [title, "=" * 24,
           f"Date: {datetime.now().isoformat(timespec='seconds')}",
           f"Command: {report.command}", "",
           "Parameters", "----------"]
    out += [f"{k}: {v}" for k, v in report.params.items()]
    out += ["", "Input", "-----",
            f"Sequences loaded: {report.n_input}",
            f"Length median/min/max: {report.median_len:.0f} / {report.min_len} / {report.max_len} bp"]
    if lf.outliers:
        verb = "dropped" if lf.action == "drop" else "kept (warn)"
        out.append(f"Length outliers ({verb}):")
        out += [f"  {sid}: {n} bp" for sid, n in lf.outliers]
    else:
        out.append("Length outliers: none")
    out += [f"Sequences used: {len(lf.kept)}", f"Reference: {report.reference_id}", "",
            "Alignment", "---------",
            f"Alignment columns: {report.aln_length}",
            f"Consensus length: {report.cons_length}",
            f"Fully conserved columns: {report.n_conserved}",
            f"Mean majority-base frequency: {report.mean_conservation:.3f}", ""]
    return out


def write_report(report: RunReport, path: Path) -> None:
    out = _header_lines(report, "VHIRPrimerLab report")
    for res in report.regions:
        out += _region_section(res, report.min_coverage) + [""]
    failed = [res.region.name for res in report.regions if res.rank.below_threshold]
    out += ["Summary", "-------"]
    if failed:
        out.append(
            f"Regions without a pair meeting --min-coverage ({report.min_coverage:g}%): "
            + ", ".join(failed)
        )
    else:
        out.append("Every region has at least one pair meeting --min-coverage.")
    out += ["", "Log", "---", *report.log_lines]
    Path(path).write_text("\n".join(out) + "\n")


CHECK_TSV_COLUMNS = [
    "name", "forward_seq", "reverse_seq", "forward_tm", "reverse_tm", "forward_gc",
    "reverse_gc", "forward_len", "reverse_len", "degenerate", *PER_PRIMER_DEGENERACY,
    "forward_tm_range", "reverse_tm_range", "located", "forward_site_mismatches", "reverse_site_mismatches",
    "amplicon_start", "amplicon_end", "amplicon_len", "coverage_pct", "n_covered",
    "n_total", "below_threshold", "not_covered", "note",
]


def _fmt_tm(tm: float | None) -> str:
    return f"{tm:.1f}" if tm is not None else ""


def _amplicon(res: CheckResult, cmap: CoordMap) -> tuple[int, int] | None:
    """Amplicon span on the reference, or None when unlocated or outside the reference."""
    if not res.located:
        return None
    try:
        start = cmap.col_to_ref(res.fwd_cols[0], "start")
        end = cmap.col_to_ref(res.rev_cols[-1], "end")
    except PipelineError:  # the reference has no bases there (e.g. a partial genome)
        return None
    return (start, end) if start <= end else None


def _note(res: CheckResult, cmap: CoordMap) -> str:
    extra = [f"amplicon not on reference {cmap.reference_id}"] if res.located and not _amplicon(res, cmap) else []
    return "; ".join(filter(None, [res.note, *extra]))


def write_check_tsv(results: list[CheckResult], cmap: CoordMap, path: Path) -> None:
    """One row per checked pair, in input order."""
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(CHECK_TSV_COLUMNS)
        for res in results:
            a, cov = res.assay, res.coverage
            amp = _amplicon(res, cmap)
            writer.writerow([
                a.name, a.forward, a.reverse, _fmt_tm(res.fwd_tm), _fmt_tm(res.rev_tm),
                f"{res.fwd_gc:.1f}", f"{res.rev_gc:.1f}", len(a.forward), len(a.reverse),
                res.degenerate, *_degeneracy(a.forward, a.reverse), _fmt_range(res.fwd_tm_range), _fmt_range(res.rev_tm_range),
                res.located,
                res.fwd.mismatches if res.located else "", res.rev.mismatches if res.located else "",
                *(amp + (amp[1] - amp[0] + 1,) if amp else ("", "", "")),
                f"{cov.coverage_pct:.1f}", len(cov.covered), len(cov.covered) + len(cov.not_covered),
                res.below_threshold,
                ",".join(f"{sid}({why})" for sid, why in cov.not_covered.items()),
                _note(res, cmap),
            ])


def _pair_section(res: CheckResult, cmap: CoordMap) -> list[str]:
    a, cov = res.assay, res.coverage
    title = f"Pair {a.name}"
    n_total = len(cov.covered) + len(cov.not_covered)
    lines = [title, "-" * len(title), f"Forward: {a.forward}", f"Reverse: {a.reverse}"]
    if res.located:
        f0, r0 = res.fwd.start, res.rev.start
        amp = _amplicon(res, cmap)
        lines += [
            f"Consensus positions: forward {f0 + 1}-{f0 + len(a.forward)}, "
            f"reverse {r0 + 1}-{r0 + len(a.reverse)}",
            f"Amplicon on the reference: {amp[0]}-{amp[1]} ({amp[1] - amp[0] + 1} bp)" if amp
            else f"Amplicon on the reference: not on reference {cmap.reference_id}",
            f"Consensus mismatches: forward {res.fwd.mismatches}, reverse {res.rev.mismatches}",
        ]
    else:
        lines.append(f"Not located: {res.note}")
    lines.append(f"Coverage: {cov.coverage_pct:.1f}% ({len(cov.covered)} of {n_total})")
    note = _note(res, cmap)
    if res.located and note:
        lines.append(f"Note: {note}")
    if res.below_threshold and cov.not_covered:
        lines.append("Why it misses sequences (count of uncovered sequences):")
        lines += [f"  {reason}: {n}" for reason, n in count_reasons(cov.not_covered)]
    return lines


def write_check_report(report: CheckReport, cmap: CoordMap, path: Path) -> None:
    out = _header_lines(report, "VHIRPrimerLab check report")
    for res in report.results:
        out += _pair_section(res, cmap) + [""]
    low = [res.assay.name for res in report.results if res.below_threshold]
    out += ["Summary", "-------"]
    if low:
        out.append(f"Pairs below --min-coverage ({report.min_coverage:g}%): " + ", ".join(low))
    else:
        out.append("Every pair meets --min-coverage.")
    out += ["", "Log", "---", *report.log_lines]
    Path(path).write_text("\n".join(out) + "\n")


TILING_TSV_COLUMNS = [
    "amplicon", "pool", "forward_name", "forward_seq", "reverse_name", "reverse_seq",
    "forward_tm", "reverse_tm", "forward_gc", "reverse_gc", "forward_len", "reverse_len",
    "degenerate", *PER_PRIMER_DEGENERACY, "forward_tm_range", "reverse_tm_range",
    "amplicon_start", "amplicon_end", "amplicon_len", "overlap_next", "coverage_pct", "n_covered", "n_total",
    "below_threshold", "not_covered",
]


def _names(index: int) -> tuple[str, str]:
    return f"VHIRPrimerLab_{index}_LEFT", f"VHIRPrimerLab_{index}_RIGHT"


def _span(first_col: int, last_col: int, cmap: CoordMap) -> tuple[int, int] | None:
    """1-based inclusive reference span of alignment columns, or None if off the reference."""
    try:
        start, end = cmap.col_to_ref(first_col, "start"), cmap.col_to_ref(last_col, "end")
    except PipelineError:
        return None
    return (start, end) if start <= end else None


def write_tiling_tsv(result: TilingResult, cmap: CoordMap, path: Path) -> None:
    """One row per amplicon, in genome order."""
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TILING_TSV_COLUMNS)
        for a in result.amplicons:
            left, right = _names(a.index)
            span = _span(a.fwd.cols[0], a.rev.cols[-1], cmap)
            cov = a.coverage
            writer.writerow([
                a.index, a.pool, left, a.fwd.seq, right, a.rev.seq,
                f"{a.fwd.tm:.1f}", f"{a.rev.tm:.1f}", f"{a.fwd.gc:.1f}", f"{a.rev.gc:.1f}",
                a.fwd.length, a.rev.length,
                a.fwd.tm_range is not None or a.rev.tm_range is not None,
                *_degeneracy(a.fwd.seq, a.rev.seq),
                _fmt_range(a.fwd.tm_range), _fmt_range(a.rev.tm_range),
                *(span if span else ("", "")), span[1] - span[0] + 1 if span else "",
                "" if a.overlap_next is None else a.overlap_next,
                f"{cov.coverage_pct:.1f}", len(cov.covered), len(cov.covered) + len(cov.not_covered),
                a.below_threshold,
                ",".join(f"{sid}({why})" for sid, why in cov.not_covered.items()),
            ])


def write_scheme_bed(result: TilingResult, cmap: CoordMap, path: Path) -> None:
    """ARTIC-style primer scheme: chrom, start, end (0-based half-open), name, pool, strand, seq."""
    lines: list[str] = []
    for a in result.amplicons:
        for name, primer, strand in zip(_names(a.index), (a.fwd, a.rev), "+-"):
            span = _span(primer.cols[0], primer.cols[-1], cmap)
            if span is None:
                log.warning("Primer %s is not on reference %s; left out of scheme.bed",
                            name, cmap.reference_id)
                continue
            lines.append("\t".join([cmap.reference_id, str(span[0] - 1), str(span[1]), name,
                                    str(a.pool), strand, primer.seq]))
    Path(path).write_text("".join(line + "\n" for line in lines))


def write_tiling_report(report: TilingReport, cmap: CoordMap, path: Path) -> None:
    res, mc = report.result, report.min_coverage
    g = res.geometry
    out = _header_lines(report, "VHIRPrimerLab tiling report")
    out += ["Tiling", "------",
            f"Amplicons: {len(res.amplicons)}",
            f"Target length: {g.target_len} bp (allowed {g.len_min}-{g.len_max}, "
            f"±{round(g.tolerance * 100)} %)",
            f"Candidate primers: {res.n_fwd} forward, {res.n_rev} reverse",
            f"Not covered at the ends: {res.uncovered_ends[0]} bp (start), "
            f"{res.uncovered_ends[1]} bp (end)", ""]
    for a in res.amplicons:
        span = _span(a.fwd.cols[0], a.rev.cols[-1], cmap)
        where = f"{span[0]}-{span[1]}  {span[1] - span[0] + 1} bp" if span else "not on reference"
        out.append(f"  {a.index} pool {a.pool}  {where}  {a.coverage.coverage_pct:.1f} %")
        if a.below_threshold:
            out.append("    Why sequences are missed (count):")
            out += [f"      {reason}: {n}" for reason, n in count_reasons(a.coverage.not_covered)]
    out += ["", "Primer dimers", "-------------"]
    out += [f"  {w}" for w in res.dimer_warnings] or ["none"]
    low = [str(a.index) for a in res.amplicons if a.below_threshold]
    out += ["", "Summary", "-------",
            f"Amplicons below {mc:g} % coverage: {', '.join(low)}" if low
            else f"All amplicons meet {mc:g} % coverage."]
    out += ["", "Log", "---", *report.log_lines]
    Path(path).write_text("\n".join(out) + "\n")
