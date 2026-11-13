from __future__ import annotations

from pathlib import Path

import pytest

from vhirprimerlab.web.params import Limits, build_argv, validate_upload
from vhirprimerlab.web.prepare import prepare_job


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


def test_genes_mode_builds_tile_region_flags() -> None:
    argv = build_argv("genes", {"tile_regions": ["F:3012-4631", "L:7074-13090"],
                                "amplicon_len": 400})
    assert argv == ["--tile-region", "F:3012-4631", "--tile-region", "L:7074-13090",
                    "--amplicon-len", "400"]


def test_genes_mode_passes_the_padding_as_the_loss_limits() -> None:
    """The page pads each span, and the engine measures the loss limits from the padded edge."""
    argv = build_argv("genes", {"tile_regions": ["F:2912-4718"], "amplicon_len": 600,
                                "overlap": 50, "max_start_loss": 100, "max_end_loss": 100})
    assert argv == ["--tile-region", "F:2912-4718", "--amplicon-len", "600", "--overlap", "50",
                    "--max-start-loss", "100", "--max-end-loss", "100"]


def test_genes_mode_rejects_tiling_only_settings() -> None:
    with pytest.raises(ValueError, match="unknown setting"):
        build_argv("genes", {"amplicons": 4})


def test_genes_mode_needs_at_least_one_region(tmp_path: Path) -> None:
    fasta = _fasta(tmp_path / "g.fa", 2, 1000).read_bytes()
    with pytest.raises(ValueError, match="Pick at least one region"):
        prepare_job("genes", {"tile_regions": []}, fasta, "", Limits(), ["mafft"])


def test_genes_mode_refuses_a_span_past_the_sequence(tmp_path: Path) -> None:
    fasta = _fasta(tmp_path / "g.fa", 2, 1000).read_bytes()
    with pytest.raises(ValueError, match="ends at 1,200"):
        prepare_job("genes", {"tile_regions": ["F:100-1200"]}, fasta, "", Limits(), ["mafft"])


def test_genes_mode_caps_the_amplicon_count(tmp_path: Path) -> None:
    fasta = _fasta(tmp_path / "g.fa", 2, 1000).read_bytes()
    settings = {"tile_regions": ["F:1-1000"], "amplicon_len": 100, "overlap": 90}
    with pytest.raises(ValueError, match=r"need about 91 amplicons.*up to 50"):
        prepare_job("genes", settings, fasta, "", Limits(max_amplicons=50), ["mafft"])


def test_genes_mode_reports_a_malformed_span(tmp_path: Path) -> None:
    fasta = _fasta(tmp_path / "g.fa", 2, 1000).read_bytes()
    with pytest.raises(ValueError):
        prepare_job("genes", {"tile_regions": ["F:abc"]}, fasta, "", Limits(), ["mafft"])


# --- a background FASTA: sequences the primers must not bind ---------------------------

FASTA = ">a\n" + "ACGTTGCAAGCTTCGGATCC" * 30 + "\n>b\n" + "GGATCCGAAGCTTGCAACGT" * 30 + "\n"
HOST = ">host\n" + "TTGACCGTAGGCATCGATGC" * 20 + "\n"


def test_a_background_is_passed_to_tiling_as_a_file_and_a_flag() -> None:
    argv, files = prepare_job("tile", {"amplicons": 2}, FASTA.encode(), "", Limits(), ["mafft"],
                              background_text=HOST)
    assert argv[argv.index("--background") + 1] == "background.fasta"
    assert files["background.fasta"] == HOST.encode()


def test_a_background_is_passed_to_check_too() -> None:
    argv, files = prepare_job("check", {}, FASTA.encode(), "p\tACGTTGCAAGCTTCGGATCC\tGGATCCGAAGCTTGCAACGT\n",
                              Limits(), ["mafft"], background_text=HOST)
    assert "--background" in argv and files["background.fasta"] == HOST.encode()


def test_no_background_means_no_flag_and_no_file() -> None:
    for text in ("", "  \n"):
        argv, files = prepare_job("tile", {"amplicons": 2}, FASTA.encode(), "", Limits(), ["mafft"],
                                  background_text=text)
        assert "--background" not in argv and "background.fasta" not in files


def test_region_design_has_no_use_for_a_background() -> None:
    with pytest.raises(ValueError, match="tiling and in check mode"):
        prepare_job("regions", {"regions": ["1-300"], "amplicon_min": 100, "amplicon_max": 250},
                    FASTA.encode(), "", Limits(), ["mafft"], background_text=HOST)


@pytest.mark.parametrize("text, message", [
    ("not fasta at all", "FASTA"),
    (">x\nACGT123\n", "invalid character"),
])
def test_a_background_that_is_not_fasta_is_refused_in_words(text: str, message: str) -> None:
    with pytest.raises(ValueError, match=message) as err:
        prepare_job("tile", {"amplicons": 2}, FASTA.encode(), "", Limits(), ["mafft"], background_text=text)
    assert "background" in str(err.value).lower() and "/tmp" not in str(err.value)


def test_a_background_too_large_to_scan_is_refused() -> None:
    big = ">big\n" + "A" * 10_000_001 + "\n"
    with pytest.raises(ValueError, match="10,000,000"):
        prepare_job("tile", {"amplicons": 2}, FASTA.encode(), "", Limits(max_upload_mb=50), ["mafft"],
                    background_text=big)


def test_a_background_over_the_upload_limit_is_refused() -> None:
    with pytest.raises(ValueError, match="larger than"):
        prepare_job("tile", {"amplicons": 2}, FASTA.encode(), "", Limits(max_upload_mb=0.001), ["mafft"],
                    background_text=">h\n" + "A" * 5000 + "\n")
