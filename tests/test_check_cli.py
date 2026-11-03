from __future__ import annotations

import csv
from pathlib import Path

import pytest

from vhirprimerlab.cli import main
from vhirprimerlab.io import CHECK_TSV_COLUMNS
from tests.synthetic import make_aligned_family, make_family, write_fasta

DESIGN = ["--region", "cons:61-540", "--region", "hv:701-1100",
          "--amplicon-min", "100", "--amplicon-max", "300"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture
def aln(tmp_path: Path) -> Path:
    return write_fasta(tmp_path / "aln.fasta", make_aligned_family())


def _rank1(design_out: Path) -> dict[str, dict[str, str]]:
    return {r["region"]: r for r in read_tsv(design_out / "pairs.tsv") if r["rank"] == "1"}


def _primer_file(path: Path, rows: dict[str, dict[str, str]], names: list[str]) -> Path:
    lines = ["name\tforward\treverse"]
    lines += [f"{n}\t{rows[n]['forward_seq']}\t{rows[n]['reverse_seq']}" for n in names]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_make_aligned_family() -> None:
    rows = make_aligned_family()
    assert len({len(s) for _, s in rows}) == 1
    assert [(i, s.replace("-", "")) for i, s in rows] == make_family()[1]


@pytest.mark.parametrize("extra", [[], ["--allow-degenerate"]], ids=["plain", "degenerate"])
def test_check_agrees_with_design(aln: Path, tmp_path: Path, extra: list[str]) -> None:
    d, c = tmp_path / "design", tmp_path / "check"
    main(["--aligned", str(aln), "-o", str(d), *DESIGN, *extra])
    rank1 = _rank1(d)
    p = _primer_file(tmp_path / "p.tsv", rank1, ["cons", "hv"])
    main(["check", "--aligned", str(aln), "--primers", str(p), "-o", str(c)])
    checked = {r["name"]: r for r in read_tsv(c / "check.tsv")}
    for name in ("cons", "hv"):
        for col in ("coverage_pct", "n_covered", "not_covered"):
            assert checked[name][col] == rank1[name][col], (name, col)


def test_check_exit_codes(aln: Path, tmp_path: Path) -> None:
    d = tmp_path / "design"
    main(["--aligned", str(aln), "-o", str(d), *DESIGN])
    rank1 = _rank1(d)
    ok = _primer_file(tmp_path / "ok.tsv", rank1, ["cons"])
    both = _primer_file(tmp_path / "both.tsv", rank1, ["cons", "hv"])
    assert main(["check", "--aligned", str(aln), "--primers", str(ok), "-o", str(tmp_path / "c1")]) == 0
    assert main(["check", "--aligned", str(aln), "--primers", str(both), "-o", str(tmp_path / "c2")]) == 2
    summary = (tmp_path / "c2" / "report.txt").read_text().split("Summary\n-------\n")[1].splitlines()[0]
    assert summary.endswith(": hv")


def test_check_bad_primer_file_exit_1(aln: Path, tmp_path: Path) -> None:
    p = tmp_path / "bad.tsv"
    p.write_text("A1\tACGTXCGTACGTACGTAC\tTTGGCCAATTGGCCAATT\n")
    out = tmp_path / "c"
    assert main(["check", "--aligned", str(aln), "--primers", str(p), "-o", str(out)]) == 1
    assert not (out / "check.tsv").exists() and not (out / "alignment.fasta").exists()


def test_check_rejects_design_option(aln: Path, tmp_path: Path) -> None:
    p = tmp_path / "p.tsv"
    p.write_text("A1\tACGTACGTACGTACGTAC\tTTGGCCAATTGGCCAATT\n")
    with pytest.raises(SystemExit) as exc:
        main(["check", "--aligned", str(aln), "--primers", str(p), "-o", str(tmp_path / "c"),
              "--region", "1-100"])
    assert exc.value.code == 1


@pytest.mark.parametrize("bad", [["--min-coverage", "101"], ["--three-prime-window", "0"]])
def test_check_bad_params_exit_1(aln: Path, tmp_path: Path, bad: list[str]) -> None:
    p = tmp_path / "p.tsv"
    p.write_text("A1\tACGTACGTACGTACGTAC\tTTGGCCAATTGGCCAATT\n")
    assert main(["check", "--aligned", str(aln), "--primers", str(p), "-o", str(tmp_path / "c"), *bad]) == 1


def test_check_outputs(aln: Path, tmp_path: Path) -> None:
    d, c = tmp_path / "design", tmp_path / "check"
    main(["--aligned", str(aln), "-o", str(d), *DESIGN])
    p = _primer_file(tmp_path / "p.tsv", _rank1(d), ["cons"])
    main(["check", "--aligned", str(aln), "--primers", str(p), "-o", str(c)])
    assert (c / "alignment.fasta").is_file() and (c / "report.txt").is_file()
    assert (c / "check.tsv").read_text().splitlines()[0].split("\t") == CHECK_TSV_COLUMNS


def test_design_help_mentions_check(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert "vhirprimerlab check --help" in capsys.readouterr().out


def test_check_reads_the_order_sheet(aln: Path, tmp_path: Path) -> None:
    d, c = tmp_path / "design", tmp_path / "check"
    main(["--aligned", str(aln), "-o", str(d), *DESIGN])
    rank1 = _rank1(d)
    main(["check", "--aligned", str(aln), "--primers", str(d / "primers.tsv"), "-o", str(c)])
    checked = {r["name"]: r for r in read_tsv(c / "check.tsv")}
    assert set(checked) == {"cons", "hv"}
    for name in ("cons", "hv"):
        assert checked[name]["forward_seq"] == rank1[name]["forward_seq"]
        assert checked[name]["coverage_pct"] == rank1[name]["coverage_pct"]
        assert checked[name]["forward_structure"]  # "ok" or what is wrong
        assert float(checked[name]["pair_dimer_dg"]) < 0
    assert "Dimers between pairs" in (c / "report.txt").read_text()
