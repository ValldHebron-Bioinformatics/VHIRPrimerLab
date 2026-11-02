from __future__ import annotations

from pathlib import Path

import pytest

from vhirprimerlab.cli import main
from tests.conftest import needs_mafft
from tests.synthetic import make_family

AMP = ["--amplicon-min", "100", "--amplicon-max", "300"]


def test_amplicon_args_required() -> None:
    with pytest.raises(SystemExit):
        main(["-i", "x.fa", "-o", "out"])


def test_input_and_aligned_exclusive() -> None:
    with pytest.raises(SystemExit):
        main(["-i", "x.fa", "--aligned", "y.fa", "-o", "out", *AMP])


def test_bad_tm_order_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--tm-min", "62", "--tm-opt", "60"])
    assert rc == 1


def test_window_longer_than_primer_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--three-prime-window", "30"])
    assert rc == 1


def test_missing_input_returns_1(tmp_path: Path) -> None:
    assert main(["-i", str(tmp_path / "missing.fa"), "-o", str(tmp_path / "o"), *AMP]) == 1


def test_unknown_reference_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--reference-id", "nope"])
    assert rc == 1


def test_bad_region_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--region", "500-100"])
    assert rc == 1


@needs_mafft
def test_region_outside_reference_returns_1(family: Path, tmp_path: Path) -> None:
    rc = main(["-i", str(family), "-o", str(tmp_path / "o"), *AMP, "--region", "1000-5000"])
    assert rc == 1


@needs_mafft
def test_reference_id_with_description(tmp_path: Path) -> None:
    _, recs = make_family()
    fa = tmp_path / "desc.fa"
    fa.write_text("".join(f">{rid.replace('seq', 'seqX')} some description\n{s}\n" for rid, s in recs))
    rc = main(["-i", str(fa), "-o", str(tmp_path / "o"), *AMP, "--reference-id", "seqX02",
               "--region", "100-500"])
    assert rc == 0


def test_argument_error_exits_1() -> None:
    with pytest.raises(SystemExit) as exc:
        main(["-i", "x.fa", "-o", "out"])  # missing --amplicon-min/max
    assert exc.value.code == 1


@needs_mafft
def test_region_too_short_on_consensus_fails_before_design(
    family: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # seq05 carries a 2 bp insertion after 570: region 560-659 is 100 bp on seq05 but only
    # 98 consensus positions, shorter than --amplicon-min 100
    out = tmp_path / "o"
    rc = main(["-i", str(family), "-o", str(out), "--reference-id", "seq05",
               "--region", "good:100-500", "--region", "short:560-659", *AMP])
    assert rc == 1
    err = capsys.readouterr().err
    assert "short" in err and "amplicon-min" in err
    assert not (out / "primers.tsv").exists()


@needs_mafft
def test_default_reference_is_first_after_filtering(tmp_path: Path) -> None:
    base, recs = make_family()
    fa = tmp_path / "f.fa"
    fa.write_text(f">partial\n{base[:300]}\n" + "".join(f">{r}\n{s}\n" for r, s in recs))
    out = tmp_path / "o"
    assert main(["-i", str(fa), "-o", str(out), "--region", "100-500", *AMP]) == 0
    assert "Reference: seq01" in (out / "report.txt").read_text()
    rc = main(["-i", str(fa), "-o", str(tmp_path / "o2"), "--reference-id", "partial",
               "--region", "100-200", *AMP])
    assert rc == 1
