from __future__ import annotations

import csv
from pathlib import Path

import pytest

from viraprime import align
from viraprime.cli import main
from tests.conftest import needs_mafft
from tests.synthetic import make_family, write_fasta

AMP = ["--amplicon-min", "100", "--amplicon-max", "300"]
COLUMNS = [
    "region", "rank", "forward_seq", "reverse_seq", "forward_tm", "reverse_tm", "forward_gc",
    "reverse_gc", "forward_len", "reverse_len", "amplicon_start", "amplicon_end", "amplicon_len",
    "coverage_pct", "n_covered", "n_total", "degenerate", "forward_degenerate_bases",
    "reverse_degenerate_bases", "forward_variants", "reverse_variants", "forward_tm_range",
    "reverse_tm_range", "penalty", "below_threshold", "not_covered",
]

pytestmark = needs_mafft


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def test_single_region_end_to_end(family: Path, tmp_path: Path) -> None:
    out = tmp_path / "o"
    assert main(["-i", str(family), "-o", str(out), "--region", "100-500", *AMP]) == 0
    for name in ("alignment.fasta", "primers.tsv", "report.txt"):
        assert (out / name).is_file()
    header, rows = read_tsv(out / "primers.tsv")
    assert header == COLUMNS
    assert 1 <= len(rows) <= 5
    first = rows[0]
    assert first["region"] == "region1" and first["rank"] == "1"
    assert first["coverage_pct"] == "100.0" and first["n_total"] == "10"
    start, end = int(first["amplicon_start"]), int(first["amplicon_end"])
    assert 100 <= start < end <= 500
    assert int(first["amplicon_len"]) == end - start + 1
    assert 100 <= int(first["amplicon_len"]) <= 300


def test_two_regions(family: Path, tmp_path: Path) -> None:
    out = tmp_path / "o"
    rc = main(["-i", str(family), "-o", str(out), "--region", "A:61-300", "--region", "B:301-540", *AMP])
    assert rc == 0
    _, rows = read_tsv(out / "primers.tsv")
    assert {r["region"] for r in rows} == {"A", "B"}
    for name in ("A", "B"):
        assert [r["rank"] for r in rows if r["region"] == name][0] == "1"
    report = (out / "report.txt").read_text()
    assert "Region A (61-300)" in report and "Region B (301-540)" in report


def test_regions_file(family: Path, tmp_path: Path) -> None:
    rf = tmp_path / "regions.tsv"
    rf.write_text("name\tstart\tend\nA\t61\t300\n")
    out = tmp_path / "o"
    assert main(["-i", str(family), "-o", str(out), "--regions-file", str(rf), *AMP]) == 0
    _, rows = read_tsv(out / "primers.tsv")
    assert {r["region"] for r in rows} == {"A"}


def test_mixed_regions_exit_2(family: Path, tmp_path: Path) -> None:
    out = tmp_path / "o"
    rc = main(["-i", str(family), "-o", str(out), "--region", "ok:100-500", "--region", "hv:701-1100", *AMP])
    assert rc == 2
    summary = (out / "report.txt").read_text().split("Summary")[1].split("Log")[0]
    assert "hv" in summary and "ok" not in summary


def test_outlier_dropped_reported(tmp_path: Path) -> None:
    base, recs = make_family()
    fa = write_fasta(tmp_path / "f.fa", recs + [("short", base[:300])])
    out = tmp_path / "o"
    assert main(["-i", str(fa), "-o", str(out), "--region", "100-500", *AMP]) == 0
    assert "short" in (out / "report.txt").read_text()
    assert "short" not in (out / "alignment.fasta").read_text()


def test_aligned_input_skips_aligner(family: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out1 = tmp_path / "o1"
    assert main(["-i", str(family), "-o", str(out1), "--region", "100-500", *AMP]) == 0
    monkeypatch.setattr(align.shutil, "which", lambda _: None)
    out2 = tmp_path / "o2"
    rc = main(["--aligned", str(out1 / "alignment.fasta"), "-o", str(out2), "--region", "100-500", *AMP])
    assert rc == 0
    assert (out1 / "primers.tsv").read_text() == (out2 / "primers.tsv").read_text()


def test_degenerate_rescues_snp(tmp_path: Path) -> None:
    _, recs = make_family()

    # a shared SNP every 15 bases across the conserved block in half the sequences,
    # so with --max-mismatches 0 every primer there needs degeneracy to reach 100 %
    def mutate(s: str) -> str:
        chars = list(s)
        for p in range(100, 500, 15):
            chars[p] = {"A": "G", "G": "A", "C": "T", "T": "C"}[chars[p]]
        return "".join(chars)

    recs = [(rid, mutate(s) if i % 2 else s) for i, (rid, s) in enumerate(recs)]
    fa = write_fasta(tmp_path / "snp.fa", recs)
    common = ["-i", str(fa), "--region", "100-500", *AMP, "--max-mismatches", "0"]
    assert main([*common, "-o", str(tmp_path / "o1")]) == 2
    out2 = tmp_path / "o2"
    assert main([*common, "-o", str(out2), "--allow-degenerate"]) == 0
    _, rows = read_tsv(out2 / "primers.tsv")
    assert rows[0]["degenerate"] == "True" and rows[0]["coverage_pct"] == "100.0"
    assert rows[0]["forward_tm_range"] or rows[0]["reverse_tm_range"]
