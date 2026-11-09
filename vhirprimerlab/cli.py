"""Command-line interface and pipeline orchestration."""

from __future__ import annotations

import argparse
import json
import logging
import os
import shlex
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, NoReturn

from Bio.Seq import reverse_complement

from vhirprimerlab.align import run_aligner
from vhirprimerlab.check import check_assay, cross_dimers, load_assays, structure_warnings
from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.degenerate import degenerate_pair, tm_range
from vhirprimerlab.design import run_primer3
from vhirprimerlab.io import (
    alignment_from_records,
    load_fasta,
    region_order,
    tiling_order,
    write_alignment,
    write_check_report,
    write_check_tsv,
    write_order_tsv,
    write_report,
    write_scheme_bed,
    write_tiling_report,
    write_tiling_tsv,
    write_tsv,
)
from vhirprimerlab.lengths import filter_by_length, ungapped_length
from vhirprimerlab.models import (
    Alignment,
    CandidatePair,
    CheckReport,
    ColumnStats,
    Consensus,
    DesignParams,
    LengthFilterResult,
    PipelineError,
    Record,
    Region,
    RegionResult,
    RunReport,
    ScoredPair,
    ScoringParams,
    TilingParams,
    TilingReport,
)
from vhirprimerlab.ranking import failure_summary, rank
from vhirprimerlab.regions import resolve_regions
from vhirprimerlab.scoring import score_pair
from vhirprimerlab.tiling import tile, tile_spans

log = logging.getLogger("vhirprimerlab")

EXIT_OK, EXIT_ERROR, EXIT_BELOW_THRESHOLD = 0, 1, 2


class _ListHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__()
        self.lines: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.lines.append(self.format(record))


class _Parser(argparse.ArgumentParser):
    """ArgumentParser whose usage errors exit 1 (exit 2 means 'below threshold')."""

    def error(self, message: str) -> NoReturn:
        self.print_usage(sys.stderr)
        self.exit(EXIT_ERROR, f"{self.prog}: error: {message}\n")


def build_parser() -> argparse.ArgumentParser:
    p = _Parser(
        prog="vhirprimerlab",
        description="Design consensus PCR primers from a multi-FASTA of related sequences. "
                    "Simplest use: vhirprimerlab -i genomes.fasta -n 10 -o out (tile the whole "
                    "genome with 10 overlapping amplicons in 2 pools).",
        epilog="To test primer pairs you already have: vhirprimerlab check --help. "
               "To run the web app locally: vhirprimerlab web",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-i", "--input", type=Path, help="multi-FASTA of unaligned sequences")
    src.add_argument("--aligned", type=Path, help="pre-aligned FASTA (skips the alignment step)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--aligner", choices=["mafft", "muscle"], default="mafft")
    p.add_argument("--reference-id", help="reference for coordinates (default: first sequence)")
    p.add_argument("-n", "--amplicons", type=int,
                   help="tile the whole genome with this many overlapping amplicons (2 pools)")
    p.add_argument("--tile-region", action="append", metavar="[NAME:]START-END",
                   help="tile this span of the reference (1-based, inclusive) instead of the "
                        "whole genome; repeatable, all spans go into one scheme")
    p.add_argument("--amplicon-len", type=int, default=400,
                   help="tiling with --tile-region: target amplicon length (bp)")
    p.add_argument("--overlap", type=int, default=50,
                   help="tiling: minimum overlap between neighbouring amplicons (bp)")
    p.add_argument("--max-start-loss", type=int, default=50,
                   help="tiling: most genome bases left outside the first amplicon's insert (bp)")
    p.add_argument("--max-end-loss", type=int, default=50,
                   help="tiling: most genome bases left outside the last amplicon's insert (bp)")
    p.add_argument("--no-degenerate", action="store_true",
                   help="tiling: plain primers only (degenerate bases are on by default)")
    p.add_argument("--progress-file", type=Path, help=argparse.SUPPRESS)

    reg = p.add_mutually_exclusive_group()
    reg.add_argument("--region", action="append", metavar="[NAME:]START-END",
                     help="target region, 1-based inclusive reference coordinates; repeatable "
                          "(default: whole reference)")
    reg.add_argument("--regions-file", type=Path, help="TSV with columns name, start, end")

    p.add_argument("--amplicon-min", type=int, help="region design: minimum amplicon length (bp)")
    p.add_argument("--amplicon-max", type=int, help="region design: maximum amplicon length (bp)")
    p.add_argument("--tm-min", type=float, default=57.0)
    p.add_argument("--tm-opt", type=float, default=60.0)
    p.add_argument("--tm-max", type=float, default=63.0)
    p.add_argument("--max-tm-diff", type=float, default=2.0, help="max Tm difference within a pair")
    p.add_argument("--primer-len-min", type=int, default=18)
    p.add_argument("--primer-len-max", type=int, default=25)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)

    p.add_argument("--length-tolerance", type=float, default=0.10,
                   help="allowed deviation from the median length, as a fraction")
    p.add_argument("--outlier-action", choices=["warn", "drop", "abort"], default="drop")

    p.add_argument("--max-mismatches", type=int, default=2, help="max mismatches per primer")
    p.add_argument("--three-prime-window", type=int, default=5,
                   help="3' bases that must match exactly")
    p.add_argument("--min-coverage", type=float, default=95.0, help="minimum coverage (%%)")
    p.add_argument("--allow-degenerate", action="store_true", help="allow IUPAC bases in primers")
    p.add_argument("--max-degeneracy", type=int,
                   help="max fold per degenerate primer (default 4 when tiling, 8 otherwise)")
    p.add_argument("--candidates", type=int, default=200, help="primer3 pairs to re-score")
    p.add_argument("--num-pairs", type=int, default=5, help="pairs to report per region")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging and tracebacks")
    return p


def _design_problems(a: argparse.Namespace) -> list[str]:
    """Checks shared by region design and tiling."""
    problems: list[str] = []
    if not a.tm_min <= a.tm_opt <= a.tm_max:
        problems.append("--tm-min <= --tm-opt <= --tm-max is required")
    if not 1 <= a.primer_len_min <= a.primer_len_max:
        problems.append("1 <= --primer-len-min <= --primer-len-max is required")
    if a.primer_len_max > 36:
        problems.append("--primer-len-max must be at most 36 (primer3's limit)")
    if not 0 <= a.gc_min <= a.gc_max <= 100:
        problems.append("0 <= --gc-min <= --gc-max <= 100 is required")
    if a.length_tolerance <= 0:
        problems.append("--length-tolerance must be > 0")
    if not 1 <= a.three_prime_window <= a.primer_len_min:
        problems.append("--three-prime-window must be between 1 and --primer-len-min")
    if a.max_mismatches < 0:
        problems.append("--max-mismatches must be >= 0")
    if not 0 <= a.min_coverage <= 100:
        problems.append("--min-coverage must be between 0 and 100")
    if a.max_tm_diff < 0:
        problems.append("--max-tm-diff must be >= 0")
    if a.max_degeneracy < 1:
        problems.append("--max-degeneracy must be >= 1")
    return problems


def check_params(a: argparse.Namespace) -> None:
    problems = _design_problems(a)
    if not 1 <= a.amplicon_min <= a.amplicon_max:
        problems.append("1 <= --amplicon-min <= --amplicon-max is required")
    for name in ("candidates", "num_pairs"):
        if getattr(a, name) < 1:
            problems.append(f"--{name.replace('_', '-')} must be >= 1")
    if problems:
        raise PipelineError("Invalid parameters: " + "; ".join(problems))


def _score_candidate(
    cand: CandidatePair,
    cons: Consensus,
    aln: Alignment,
    sparams: ScoringParams,
    allow_degenerate: bool,
    max_fold: int,
    gc_limits: tuple[float, float] | None = None,
) -> ScoredPair:
    f, r = cand.forward, cand.reverse
    fwd_cols = cons.cons_to_col[f.start: f.start + f.length]
    rev_cols = cons.cons_to_col[r.start: r.start + r.length]
    fwd_site, rev_site = f.seq, reverse_complement(r.seq)
    cov = score_pair(fwd_site, fwd_cols, rev_site, rev_cols, aln, sparams)
    fwd_seq, rev_seq, degenerate = f.seq, r.seq, False
    if allow_degenerate and cov.not_covered:
        new_f, new_r, new_cov = degenerate_pair(
            fwd_site, fwd_cols, rev_site, rev_cols, aln, sparams, max_fold, gc_limits
        )
        if (new_f, new_r) != (fwd_site, rev_site):
            fwd_seq, rev_seq, cov, degenerate = new_f, reverse_complement(new_r), new_cov, True
    return ScoredPair(
        pair=cand, fwd_seq=fwd_seq, rev_seq=rev_seq, fwd_cols=fwd_cols, rev_cols=rev_cols,
        coverage=cov, degenerate=degenerate,
        fwd_tm_range=tm_range(fwd_seq) if degenerate else None,
        rev_tm_range=tm_range(rev_seq) if degenerate else None,
    )


def _design_region(
    region: Region,
    included: tuple[int, int],
    cons: Consensus,
    aln: Alignment,
    a: argparse.Namespace,
    dparams: DesignParams,
    sparams: ScoringParams,
) -> RegionResult:
    try:
        design = run_primer3(cons, included, dparams)
    except PipelineError as exc:
        raise PipelineError(f"Region '{region.name}': {exc}") from exc
    scored = [
        _score_candidate(
            c, cons, aln, sparams, a.allow_degenerate, a.max_degeneracy, (a.gc_min, a.gc_max)
        )
        for c in design.pairs
    ]
    ranked = rank(scored, a.min_coverage, a.num_pairs)
    n_passing = sum(1 for s in scored if s.coverage.coverage_pct >= a.min_coverage)
    failures = failure_summary(scored) if ranked.below_threshold else []
    if not design.pairs:
        log.warning("Region %s: primer3 found no candidate pairs", region.name)
    elif ranked.below_threshold:
        log.warning(
            "Region %s: no pair reached %g%% coverage (best %.1f%%)",
            region.name, a.min_coverage, ranked.pairs[0].coverage.coverage_pct,
        )
    else:
        log.info(
            "Region %s: %d of %d candidate pairs meet %g%% coverage",
            region.name, n_passing, len(design.pairs), a.min_coverage,
        )
    return RegionResult(region, included, len(design.pairs), n_passing, ranked, design.explain, failures)


@dataclass
class _Inputs:
    records: list[Record]
    lengths: list[int]
    lf: LengthFilterResult
    reference_id: str
    ref_len: int


@dataclass
class _Aligned:
    aln: Alignment
    stats: list[ColumnStats]
    cons: Consensus
    cmap: CoordMap


def check_span(region: Region, ref_length: int) -> None:
    """Refuse a --tile-region span that is backwards or not on the reference."""
    if region.start < 1 or region.end > ref_length:
        raise PipelineError(
            f"Region '{region.name}' ({region.start}-{region.end}) is outside the reference "
            f"(1-{ref_length})"
        )
    if region.start >= region.end:
        raise PipelineError(f"Region '{region.name}': start must be smaller than end")


def check_tiling_params(a: argparse.Namespace) -> None:
    problems = _design_problems(a)
    if a.amplicons is not None and a.amplicons < 1:
        problems.append("-n must be >= 1")
    if a.tile_region and a.amplicon_len <= a.overlap:
        problems.append("--amplicon-len must be greater than --overlap")
    if a.overlap < 0:
        problems.append("--overlap must be >= 0")
    if a.max_start_loss < 0:
        problems.append("--max-start-loss must be >= 0")
    if a.max_end_loss < 0:
        problems.append("--max-end-loss must be >= 0")
    if problems:
        raise PipelineError("Invalid parameters: " + "; ".join(problems))


_progress_reporter: Callable[[str, float], None] | None = None


def set_progress_reporter(fn: Callable[[str, float], None] | None) -> None:
    """Also send every progress update to fn (a host program, e.g. the browser worker)."""
    global _progress_reporter
    _progress_reporter = fn


def make_progress(path: Path | None) -> Callable[[str, float], None]:
    """Report progress as {"stage", "progress"} JSON in path (atomically) and to the reporter."""
    reporter = _progress_reporter
    if path is None and reporter is None:
        return lambda stage, fraction: None

    def emit(stage: str, fraction: float) -> None:
        if reporter is not None:
            reporter(stage, round(fraction, 3))
        if path is not None:
            tmp = path.with_name(path.name + ".tmp")
            tmp.write_text(json.dumps({"stage": stage, "progress": round(fraction, 3)}))
            os.replace(tmp, path)

    return emit


def _report_params(a: argparse.Namespace) -> dict[str, object]:
    # `amplicons` is -n, which tiling chosen spans derives per span: it is None there, and
    # "amplicons: None" in the report would read as a setting rather than as its absence.
    return {k: (str(v) if isinstance(v, Path) else v)
            for k, v in vars(a).items() if k != "progress_file" and (k, v) != ("amplicons", None)}


def _load_inputs(a: argparse.Namespace) -> _Inputs:
    """Load and validate the input, filter by length and pick the reference."""
    if a.input is not None:
        records = load_fasta(a.input)
    else:
        records = load_fasta(a.aligned, allow_gaps=True)
        alignment_from_records(records)  # validates equal row lengths
    if a.reference_id is not None and a.reference_id not in {r.id for r in records}:
        raise PipelineError(f"--reference-id '{a.reference_id}' is not in the input")
    lengths = [ungapped_length(r.seq) for r in records]
    lf = filter_by_length(records, a.length_tolerance, a.outlier_action, a.reference_id)
    # default reference: first sequence left after length filtering
    reference_id = a.reference_id or lf.kept[0].id
    ref_len = ungapped_length(next(r.seq for r in lf.kept if r.id == reference_id))
    return _Inputs(records, lengths, lf, reference_id, ref_len)


def _align_inputs(a: argparse.Namespace, inp: _Inputs) -> _Aligned:
    """Align (or take the given alignment), write it, and build consensus and coordinates."""
    kept = inp.lf.kept
    aln = run_aligner(kept, a.aligner) if a.input is not None else alignment_from_records(kept)
    a.outdir.mkdir(parents=True, exist_ok=True)
    write_alignment(aln, a.outdir / "alignment.fasta")
    stats = score_columns(aln)
    cons = build_consensus(stats)
    return _Aligned(aln, stats, cons, CoordMap(aln, inp.reference_id, cons))


def run(a: argparse.Namespace, command: str, log_lines: list[str]) -> int:
    check_params(a)
    progress = make_progress(a.progress_file)
    inp = _load_inputs(a)
    regions = resolve_regions(a.region, a.regions_file, inp.ref_len)
    progress("aligning", 0.05)
    al = _align_inputs(a, inp)
    records, lengths, lf, reference_id = inp.records, inp.lengths, inp.lf, inp.reference_id
    aln, stats, cons, cmap = al.aln, al.stats, al.cons, al.cmap
    # validate every region (incl. its consensus length) before any primer design
    included: dict[str, tuple[int, int]] = {}
    for region in regions:
        cmap.validate_region(region, a.amplicon_min)
        included[region.name] = cmap.region_to_included(region)
        if included[region.name][1] < a.amplicon_min:
            raise PipelineError(
                f"Region '{region.name}' spans only {included[region.name][1]} consensus "
                f"positions (the reference has insertions there), shorter than "
                f"--amplicon-min ({a.amplicon_min})"
            )

    dparams = DesignParams(
        a.tm_min, a.tm_opt, a.tm_max, a.max_tm_diff, a.primer_len_min, a.primer_len_max,
        a.gc_min, a.gc_max, a.amplicon_min, a.amplicon_max, a.candidates,
    )
    sparams = ScoringParams(a.max_mismatches, a.three_prime_window)
    progress("scoring", 0.3)
    results = [
        _design_region(r, included[r.name], cons, aln, a, dparams, sparams) for r in regions
    ]

    write_tsv(results, cmap, a.outdir / "pairs.tsv")
    write_order_tsv(region_order(results), a.outdir / "primers.tsv")
    failed = [res.region.name for res in results if res.rank.below_threshold]
    if failed:
        log.warning("Regions without a pair meeting --min-coverage: %s", ", ".join(failed))
    params = _report_params(a)
    report = RunReport(
        command=command, params=params, n_input=len(records),
        median_len=float(statistics.median(lengths)), min_len=min(lengths), max_len=max(lengths),
        length=lf, reference_id=reference_id, aln_length=aln.length, cons_length=len(cons.seq),
        n_conserved=sum(1 for s in stats if s.majority_freq == 1.0),
        mean_conservation=statistics.fmean(s.majority_freq for s in stats),
        min_coverage=a.min_coverage, regions=results, log_lines=log_lines,
    )
    write_report(report, a.outdir / "report.txt")
    log.info("Wrote results to %s", a.outdir)
    progress("done", 1.0)
    return EXIT_BELOW_THRESHOLD if failed else EXIT_OK


def run_tiling(a: argparse.Namespace, command: str, log_lines: list[str]) -> int:
    check_tiling_params(a)
    progress = make_progress(a.progress_file)
    inp = _load_inputs(a)
    spans = resolve_regions(a.tile_region, None, inp.ref_len) if a.tile_region else None
    for span in spans or []:
        check_span(span, inp.ref_len)
    progress("aligning", 0.05)
    al = _align_inputs(a, inp)
    params = TilingParams(
        n_amplicons=a.amplicons if a.amplicons is not None else 1,  # spans work it out
        overlap=a.overlap, tm_min=a.tm_min, tm_opt=a.tm_opt,
        tm_max=a.tm_max, max_tm_diff=a.max_tm_diff, len_min=a.primer_len_min,
        len_max=a.primer_len_max, gc_min=a.gc_min, gc_max=a.gc_max,
        max_mismatches=a.max_mismatches, window=a.three_prime_window,
        degenerate=not a.no_degenerate, max_degeneracy=a.max_degeneracy,
        min_coverage=a.min_coverage, max_start_loss=a.max_start_loss,
        max_end_loss=a.max_end_loss, amplicon_len=a.amplicon_len,
    )
    if spans:
        result = tile_spans(al.aln, al.cons, al.stats, al.cmap, params, spans, progress)
    else:
        result = tile(al.aln, al.cons, al.stats, al.cmap, params, progress)
    progress("writing", 0.95)
    write_tiling_tsv(result, al.cmap, a.outdir / "amplicons.tsv")
    write_order_tsv(tiling_order(result), a.outdir / "primers.tsv")
    write_scheme_bed(result, al.cmap, a.outdir / "scheme.bed")
    low = [str(x.index) for x in result.amplicons if x.below_threshold]
    if low:
        log.warning("Amplicons below --min-coverage: %s", ", ".join(low))
    else:
        log.info("All %d amplicons meet %g%% coverage", len(result.amplicons), a.min_coverage)
    report = TilingReport(
        command=command, params=_report_params(a), n_input=len(inp.records),
        median_len=float(statistics.median(inp.lengths)), min_len=min(inp.lengths),
        max_len=max(inp.lengths), length=inp.lf, reference_id=inp.reference_id,
        aln_length=al.aln.length, cons_length=len(al.cons.seq),
        n_conserved=sum(1 for s in al.stats if s.majority_freq == 1.0),
        mean_conservation=statistics.fmean(s.majority_freq for s in al.stats),
        min_coverage=a.min_coverage, result=result, log_lines=log_lines,
    )
    write_tiling_report(report, al.cmap, a.outdir / "report.txt")
    log.info("Wrote results to %s", a.outdir)
    progress("done", 1.0)
    return EXIT_BELOW_THRESHOLD if low else EXIT_OK


def build_check_parser() -> argparse.ArgumentParser:
    p = _Parser(
        prog="vhirprimerlab check",
        description="Check primer pairs you already have against every sequence of a family.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-i", "--input", type=Path, help="multi-FASTA of unaligned sequences")
    src.add_argument("--aligned", type=Path, help="pre-aligned FASTA (skips the alignment step)")
    p.add_argument("--primers", type=Path, required=True,
                   help="TSV with columns name, forward, reverse (5'->3', IUPAC allowed)")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--aligner", choices=["mafft", "muscle"], default="mafft")
    p.add_argument("--reference-id", help="reference for coordinates (default: first sequence)")
    p.add_argument("--length-tolerance", type=float, default=0.10,
                   help="allowed deviation from the median length, as a fraction")
    p.add_argument("--outlier-action", choices=["warn", "drop", "abort"], default="drop")
    p.add_argument("--max-mismatches", type=int, default=2, help="max mismatches per primer")
    p.add_argument("--three-prime-window", type=int, default=5,
                   help="3' bases that must match exactly")
    p.add_argument("--min-coverage", type=float, default=95.0, help="minimum coverage (%%)")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging and tracebacks")
    p.add_argument("--progress-file", type=Path, help=argparse.SUPPRESS)
    return p


def check_check_params(a: argparse.Namespace) -> None:
    problems: list[str] = []
    if a.max_mismatches < 0:
        problems.append("--max-mismatches must be >= 0")
    if a.three_prime_window < 1:
        problems.append("--three-prime-window must be >= 1")
    if not 0 <= a.min_coverage <= 100:
        problems.append("--min-coverage must be between 0 and 100")
    if a.length_tolerance <= 0:
        problems.append("--length-tolerance must be > 0")
    if problems:
        raise PipelineError("Invalid parameters: " + "; ".join(problems))


def run_check(a: argparse.Namespace, command: str, log_lines: list[str]) -> int:
    check_check_params(a)
    progress = make_progress(a.progress_file)
    inp = _load_inputs(a)
    assays = load_assays(a.primers, a.three_prime_window)  # before alignment: fail fast
    progress("aligning", 0.05)
    al = _align_inputs(a, inp)
    progress("scoring", 0.5)
    sparams = ScoringParams(a.max_mismatches, a.three_prime_window)
    results = [check_assay(x, al.aln, al.cons, sparams, a.min_coverage) for x in assays]
    for res in results:
        if not res.located:
            log.warning("Pair %s: not located (%s)", res.assay.name, res.note)
        else:
            log.info("Pair %s: %.1f%% coverage", res.assay.name, res.coverage.coverage_pct)
    progress("checking primer dimers", 0.8)
    dimers = cross_dimers(assays)
    for line in structure_warnings(results, dimers):
        log.warning("%s", line)
    write_check_tsv(results, al.cmap, a.outdir / "check.tsv")
    low = [r.assay.name for r in results if r.below_threshold]
    if low:
        log.warning("Pairs below --min-coverage: %s", ", ".join(low))
    params = _report_params(a)
    report = CheckReport(
        command=command, params=params, n_input=len(inp.records),
        median_len=float(statistics.median(inp.lengths)), min_len=min(inp.lengths),
        max_len=max(inp.lengths), length=inp.lf, reference_id=inp.reference_id,
        aln_length=al.aln.length, cons_length=len(al.cons.seq),
        n_conserved=sum(1 for s in al.stats if s.majority_freq == 1.0),
        mean_conservation=statistics.fmean(s.majority_freq for s in al.stats),
        min_coverage=a.min_coverage, results=results, log_lines=log_lines, dimers=dimers,
    )
    write_check_report(report, al.cmap, a.outdir / "report.txt")
    log.info("Wrote results to %s", a.outdir)
    progress("done", 1.0)
    return EXIT_BELOW_THRESHOLD if low else EXIT_OK


TILING_CONFLICTS = (("tile_region", "--tile-region"), ("region", "--region"), ("regions_file", "--regions-file"),
                    ("amplicon_min", "--amplicon-min"), ("amplicon_max", "--amplicon-max"),
                    ("allow_degenerate", "--allow-degenerate"))


def parse_args(argv: list[str]) -> tuple[argparse.Namespace, Callable[..., int], str]:
    """Parse argv for check, tiling or region mode; usage errors exit 1."""
    if argv[:1] == ["check"]:
        args = build_check_parser().parse_args(argv[1:])
        return args, run_check, "vhirprimerlab check " + shlex.join(argv[1:])
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.tile_region or args.amplicons is not None:
        lead, lead_flag = (("tile_region", "--tile-region") if args.tile_region
                           else ("amplicons", "-n"))
        clash = [flag for name, flag in (("amplicons", "-n"), *TILING_CONFLICTS)
                 if name != lead and getattr(args, name)]
        if clash:
            parser.error(f"{lead_flag} cannot be combined with {', '.join(clash)}")
        runner: Callable[..., int] = run_tiling
    elif args.no_degenerate:
        parser.error("--no-degenerate only applies to tiling (-n or --tile-region)")
    elif args.amplicon_min is None or args.amplicon_max is None:
        parser.error("use -n N to tile the whole genome with N amplicons, or --amplicon-min/"
                     "--amplicon-max (and optionally --region) for region design")
    else:
        runner = run
    if args.max_degeneracy is None:
        args.max_degeneracy = 4 if runner is run_tiling else 8
    return args, runner, "vhirprimerlab " + shlex.join(argv)


def validate_args(a: argparse.Namespace) -> None:
    """The parameter checks the chosen mode runs first (raises PipelineError)."""
    if hasattr(a, "primers"):
        check_check_params(a)
    elif a.amplicons is not None or a.tile_region:
        check_tiling_params(a)
    else:
        check_params(a)


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] == ["web"]:
        try:
            from vhirprimerlab.web.server import main as web_main
        except ImportError:
            print("vhirprimerlab web needs: pip install 'vhirprimerlab[web]'", file=sys.stderr)
            return EXIT_ERROR
        return web_main(argv[1:])
    args, runner, command = parse_args(argv)

    fmt = logging.Formatter("%(levelname)s: %(message)s")
    stderr_handler = logging.StreamHandler(sys.stderr)
    memory = _ListHandler()
    for h in (stderr_handler, memory):
        h.setFormatter(fmt)
        log.addHandler(h)
    log.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    log.propagate = False
    try:
        return runner(args, command, memory.lines)
    except PipelineError as exc:
        log.error("%s", exc)
        if args.verbose:
            log.exception("Traceback:")
        return EXIT_ERROR
    finally:
        log.removeHandler(stderr_handler)
        log.removeHandler(memory)


if __name__ == "__main__":
    sys.exit(main())
