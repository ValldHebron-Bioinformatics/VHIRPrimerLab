from __future__ import annotations

from pathlib import Path

import pytest

from vhirprimerlab.web.params import Limits, build_argv, validate_upload


def test_tile_argv() -> None:
    assert build_argv("tile", {"amplicons": 8, "degenerate": False}) == ["-n", "8", "--no-degenerate"]
    assert build_argv("tile", {"amplicons": "8", "degenerate": True, "tm_opt": 61.5}) == [
        "-n", "8", "--tm-opt", "61.5"]
    assert build_argv("tile", {"amplicons": 8, "reference_id": None}) == ["-n", "8"]


def test_regions_and_check_argv() -> None:
    assert build_argv("regions", {"regions": ["cons:61-540", "701-1100"], "amplicon_min": 100,
                                  "amplicon_max": 300, "allow_degenerate": True}) == [
        "--region", "cons:61-540", "--region", "701-1100", "--amplicon-min", "100",
        "--amplicon-max", "300", "--allow-degenerate"]
    assert build_argv("check", {"min_coverage": 90, "aligner": "muscle"}) == [
        "--min-coverage", "90.0", "--aligner", "muscle"]


@pytest.mark.parametrize("mode, params, msg", [
    ("tile", {"outdir": "/etc"}, "unknown setting: outdir"),
    ("check", {"amplicons": 5}, "unknown setting: amplicons"),
    ("tile", {"aligner": "clustal"}, "mafft or muscle"),
    ("tile", {"amplicons": "ten"}, "whole number"),
    ("tile", {"amplicons": True}, "number"),
    ("tile", {"degenerate": "no"}, "true or false"),
    ("nope", {}, "unknown mode"),
])
def test_bad_settings(mode: str, params: dict, msg: str) -> None:
    with pytest.raises(ValueError, match=msg):
        build_argv(mode, params)


def _fasta(path: Path, n: int, length: int) -> Path:
    path.write_text("".join(f">s{i}\n{'ACGT' * (length // 4)}\n" for i in range(n)))
    return path


def test_validate_upload(tmp_path: Path) -> None:
    limits = Limits(max_seqs=5, max_seq_len=1000)
    assert validate_upload(_fasta(tmp_path / "ok.fa", 5, 1000), limits) == (5, 1000)
    with pytest.raises(ValueError, match="Too many sequences: 6"):
        validate_upload(_fasta(tmp_path / "many.fa", 6, 100), limits)
    with pytest.raises(ValueError, match="1,004 bp"):
        validate_upload(_fasta(tmp_path / "long.fa", 2, 1004), limits)
    bad = tmp_path / "bad.fa"
    bad.write_text("hello, this is not a FASTA file\n")
    with pytest.raises(ValueError) as exc:
        validate_upload(bad, limits)
    assert str(tmp_path) not in str(exc.value)


def test_limits_from_env(monkeypatch) -> None:
    monkeypatch.setenv("VHIRPRIMERLAB_MAX_SEQS", "50")
    monkeypatch.setenv("VHIRPRIMERLAB_TIMEOUT_S", "30")
    limits = Limits.from_env()
    assert limits.max_seqs == 50 and limits.timeout_s == 30.0 and limits.workers == 2


@pytest.mark.parametrize("value", ["nan", "inf", float("inf")])
def test_non_finite_numbers_rejected(value) -> None:
    with pytest.raises(ValueError, match="number"):
        build_argv("tile", {"max_tm_diff": value})
