from __future__ import annotations

import csv
import dataclasses
from pathlib import Path

import pytest

from tests.synthetic import make_genome_family
from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.io import (
    TILING_TSV_COLUMNS,
    alignment_from_records,
    write_scheme_bed,
    write_tiling_report,
    write_tiling_tsv,
)
from vhirprimerlab.models import LengthFilterResult, PairCoverage, Record, TilingParams, TilingReport
from vhirprimerlab.tiling import tile

PARAMS = TilingParams(
    n_amplicons=5, overlap=50, tm_min=57, tm_opt=60, tm_max=63, max_tm_diff=2,
    len_min=18, len_max=25, gc_min=40, gc_max=60, max_mismatches=2, window=5,
    degenerate=True, max_degeneracy=4, min_coverage=95,
)


@pytest.fixture(scope="module")
def run():
    records = [Record(i, s) for i, s in make_genome_family()]
    aln = alignment_from_records(records)
    stats = score_columns(aln)
    cons = build_consensus(stats)
    cmap = CoordMap(aln, "g01", cons)
    return tile(aln, cons, stats, cmap, PARAMS), cmap, records


def test_tsv(run, tmp_path: Path) -> None:
    result, cmap, _ = run
    write_tiling_tsv(result, cmap, tmp_path / "primers.tsv")
    with (tmp_path / "primers.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        fh.seek(0)
        assert fh.readline().rstrip("\n").split("\t") == TILING_TSV_COLUMNS
    assert len(rows) == 5
    assert rows[0]["forward_name"] == "VHIRPrimerLab_1_LEFT" and rows[0]["reverse_name"] == "VHIRPrimerLab_1_RIGHT"
    assert [r["pool"] for r in rows] == ["1", "2", "1", "2", "1"]
    assert rows[-1]["overlap_next"] == "" and int(rows[0]["overlap_next"]) >= 50
    assert int(rows[0]["amplicon_len"]) == int(rows[0]["amplicon_end"]) - int(rows[0]["amplicon_start"]) + 1
    assert rows[0]["coverage_pct"] == "100.0" and rows[0]["n_total"] == "12"


def test_bed(run, tmp_path: Path) -> None:
    result, cmap, _ = run
    write_scheme_bed(result, cmap, tmp_path / "scheme.bed")
    lines = [ln.split("\t") for ln in (tmp_path / "scheme.bed").read_text().splitlines()]
    assert len(lines) == 10 and all(len(ln) == 7 for ln in lines)
    assert [ln[5] for ln in lines] == ["+", "-"] * 5
    assert [ln[4] for ln in lines] == ["1", "1", "2", "2", "1", "1", "2", "2", "1", "1"]
    for ln, primer in zip(lines, [p for a in result.amplicons for p in (a.fwd, a.rev)]):
        assert ln[0] == "g01"
        start, end = int(ln[1]), int(ln[2])
        assert start < end and end - start == primer.length
        assert ln[6] == primer.seq
    assert lines[0][3] == "VHIRPrimerLab_1_LEFT" and lines[1][3] == "VHIRPrimerLab_1_RIGHT"


def _report(result, records) -> TilingReport:
    lf = LengthFilterResult(kept=records, outliers=[], median=6000.0, action="drop")
    return TilingReport(
        command="vhirprimerlab -n 5", params={"amplicons": 5}, n_input=12, median_len=6000.0,
        min_len=6000, max_len=6000, length=lf, reference_id="g01", aln_length=6000,
        cons_length=6000, n_conserved=5000, mean_conservation=0.99, min_coverage=95.0,
        result=result, log_lines=["INFO: hello"],
    )


def test_report_all_pass(run, tmp_path: Path) -> None:
    result, cmap, records = run
    write_tiling_report(_report(result, records), cmap, tmp_path / "report.txt")
    text = (tmp_path / "report.txt").read_text()
    assert text.startswith("VHIRPrimerLab tiling report")
    for needle in ("Target length:", "Candidate primers:", "Not covered at the ends:",
                   "Primer dimers", "All amplicons meet 95", "INFO: hello"):
        assert needle in text


def test_report_below_threshold(run, tmp_path: Path) -> None:
    result, cmap, records = run
    low = PairCoverage(covered=[f"g{i:02d}" for i in range(1, 7)],
                       not_covered={f"g{i:02d}": "fwd 3' mismatch" for i in range(7, 13)},
                       coverage_pct=50.0)
    amps = list(result.amplicons)
    amps[2] = dataclasses.replace(amps[2], coverage=low, below_threshold=True)
    write_tiling_report(_report(dataclasses.replace(result, amplicons=amps), records), cmap,
                        tmp_path / "report.txt")
    text = (tmp_path / "report.txt").read_text()
    assert "Amplicons below 95 % coverage: 3" in text
    assert "fwd 3' mismatch: 6" in text


def test_tsv_per_primer_degeneracy(run, tmp_path: Path) -> None:
    from vhirprimerlab.degenerate import degenerate_bases, fold, gc_percent
    result, cmap, _ = run
    write_tiling_tsv(result, cmap, tmp_path / "primers.tsv")
    with (tmp_path / "primers.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for side in ("forward", "reverse"):
            seq = r[f"{side}_seq"]
            assert r[f"{side}_degenerate_bases"] == str(degenerate_bases(seq))
            assert r[f"{side}_variants"] == str(fold(seq))
            assert r[f"{side}_gc"] == f"{gc_percent(seq):.1f}"


def test_primer_summary_format() -> None:
    from vhirprimerlab.io import primer_summary
    assert primer_summary("ACGTACGTACGTACGTACGT", 55.04, None) == (
        "ACGTACGTACGTACGTACGT  20 nt, Tm 55.0 °C, GC 50.0 %, no degenerate bases")
    assert primer_summary("ACGTRCGTACGTACGTACGY", 55.0, (54.0, 56.25)) == (
        "ACGTRCGTACGTACGTACGY  20 nt, Tm 55.0 °C (54.0-56.2), GC 55.0 %, "
        "2 degenerate bases (4 variants)")
    assert "Tm not computed" in primer_summary("NNNNNNNNNNNNNNNNNNNN", None, None)


def test_report_lists_every_primer(run, tmp_path: Path) -> None:
    from vhirprimerlab.io import primer_summary
    result, cmap, records = run
    write_tiling_report(_report(result, records), cmap, tmp_path / "report.txt")
    text = (tmp_path / "report.txt").read_text()
    for a in result.amplicons:
        assert f"VHIRPrimerLab_{a.index}_LEFT   {primer_summary(a.fwd.seq, a.fwd.tm, a.fwd.tm_range)}" in text
        assert f"VHIRPrimerLab_{a.index}_RIGHT  {primer_summary(a.rev.seq, a.rev.tm, a.rev.tm_range)}" in text
