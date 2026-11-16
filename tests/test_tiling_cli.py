from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from tests.synthetic import make_aligned_family, make_genome_family, write_fasta
from vhirprimerlab.cli import main


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def tiled(tmp_path_factory):
    d = tmp_path_factory.mktemp("tile")
    aln = write_fasta(d / "aln.fasta", make_genome_family())
    out = d / "out"
    rc = main(["--aligned", str(aln), "-n", "5", "-o", str(out), "--progress-file", str(d / "p.json")])
    return rc, aln, out, d


def test_tiling_run_outputs(tiled) -> None:
    rc, _, out, _ = tiled
    assert rc == 0
    assert len(read_tsv(out / "amplicons.tsv")) == 5
    order = read_tsv(out / "primers.tsv")
    assert [list(r) for r in order[:1]] == [["name", "sequence"]]
    assert [r["name"] for r in order[:2]] == ["VHIRPrimerLab_1_LEFT", "VHIRPrimerLab_1_RIGHT"]
    assert len(order) == 10
    assert len((out / "scheme.bed").read_text().splitlines()) == 10
    assert (out / "report.txt").is_file() and (out / "alignment.fasta").is_file()


def test_progress_file(tiled) -> None:
    assert json.loads((tiled[3] / "p.json").read_text()) == {"stage": "done", "progress": 1.0}


def test_tiling_agrees_with_check(tiled, tmp_path: Path) -> None:
    _, aln, out, _ = tiled
    rows = read_tsv(out / "amplicons.tsv")
    p = tmp_path / "p.tsv"
    p.write_text("".join(f"a{r['amplicon']}\t{r['forward_seq']}\t{r['reverse_seq']}\n" for r in rows))
    main(["check", "--aligned", str(aln), "--primers", str(p), "-o", str(tmp_path / "c")])
    checked = {r["name"]: r for r in read_tsv(tmp_path / "c" / "check.tsv")}
    for r in rows:
        c = checked[f"a{r['amplicon']}"]
        for col in ("coverage_pct", "n_covered", "not_covered", "amplicon_start", "amplicon_end"):
            assert c[col] == r[col], (r["amplicon"], col)


def _usage_error(argv: list[str], capsys) -> str:
    with pytest.raises(SystemExit) as exc:
        main(argv)
    assert exc.value.code == 1
    return capsys.readouterr().err


def test_n_with_region_rejected(capsys) -> None:
    err = _usage_error(["-i", "x.fa", "-o", "o", "-n", "5", "--region", "1-100"], capsys)
    assert "-n cannot be combined with --region" in err


def test_no_mode_message(capsys) -> None:
    assert "use -n N" in _usage_error(["-i", "x.fa", "-o", "o"], capsys)


def test_no_degenerate_without_n_rejected(capsys) -> None:
    err = _usage_error(["-i", "x.fa", "-o", "o", "--amplicon-min", "100", "--amplicon-max", "300",
                        "--no-degenerate"], capsys)
    assert "--no-degenerate" in err


def test_region_mode_still_full_reference(tmp_path: Path) -> None:
    aln = write_fasta(tmp_path / "a.fasta", make_aligned_family())
    out = tmp_path / "o"
    main(["--aligned", str(aln), "-o", str(out), "--amplicon-min", "100", "--amplicon-max", "300"])
    assert {r["region"] for r in read_tsv(out / "pairs.tsv")} == {"full"}
    assert [r["name"] for r in read_tsv(out / "primers.tsv")] == ["full_LEFT", "full_RIGHT"]


def test_exit_2_below_threshold(tmp_path: Path) -> None:
    aln = write_fasta(tmp_path / "a.fasta", make_genome_family(hyper=(0, 6000)))
    out = tmp_path / "o"
    rc = main(["--aligned", str(aln), "-n", "5", "-o", str(out), "--no-degenerate",
               "--min-coverage", "100"])
    assert rc == 2 and (out / "primers.tsv").is_file() and (out / "scheme.bed").is_file()


def test_n_too_large(tmp_path: Path, capsys) -> None:
    aln = write_fasta(tmp_path / "a.fasta", make_genome_family())
    assert main(["--aligned", str(aln), "-n", "40", "-o", str(tmp_path / "o")]) == 1
    assert "use at most" in capsys.readouterr().err


def test_bad_tiling_param_exit_1(tmp_path: Path) -> None:
    aln = write_fasta(tmp_path / "a.fasta", make_genome_family())
    assert main(["--aligned", str(aln), "-n", "0", "-o", str(tmp_path / "o")]) == 1
    assert main(["--aligned", str(aln), "-n", "5", "--overlap", "-1", "-o", str(tmp_path / "o")]) == 1


def test_span_tiling_end_to_end(tmp_path) -> None:
    aln = write_fasta(tmp_path / "aln.fasta", make_genome_family())
    out = tmp_path / "out"
    # the synthetic genome has no acceptable primer within 50 nt of position 4001
    code = main(["--aligned", str(aln), "--tile-region", "F:1001-3000",
                 "--tile-region", "SH:4001-4600", "--amplicon-len", "700",
                 "--max-start-loss", "100", "-o", str(out)])
    assert code == 0
    names = [r["forward_name"] for r in read_tsv(out / "amplicons.tsv")]
    assert names[:2] == ["F_1_LEFT", "F_2_LEFT"] and names[-1] == "SH_LEFT"
    assert (out / "scheme.bed").read_text().count("SH_RIGHT") == 1
    assert "Target amplicon length: 700 bp" in (out / "report.txt").read_text()


def test_report_has_a_specificity_section_and_a_clean_run_says_so(tiled) -> None:
    _, _, out, _ = tiled
    text = (out / "report.txt").read_text()
    assert "Specificity\n-----------" in text
    assert "No secondary binding sites found." in text


def test_background_hit_is_reported_and_the_run_still_succeeds(tiled, tmp_path: Path) -> None:
    _, aln, out, _ = tiled
    first = read_tsv(out / "primers.tsv")[0]
    assert first["name"] == "VHIRPrimerLab_1_LEFT"
    host = tmp_path / "host.fasta"
    write_fasta(host, [("host_gene", "ACGT" * 40 + first["sequence"] + "TTGA" * 40)])
    out2 = tmp_path / "out2"
    rc = main(["--aligned", str(aln), "-n", "5", "-o", str(out2), "--background", str(host)])
    assert rc == 0
    text = (out2 / "report.txt").read_text()
    assert "background sequence(s)" in text
    assert "VHIRPrimerLab_1_LEFT: 1 site(s) in the background sequences" in text
    assert "No secondary binding sites found." not in text


def test_an_unreadable_or_oversized_background_is_refused_before_the_run(tiled, tmp_path: Path) -> None:
    _, aln, _, _ = tiled
    assert main(["--aligned", str(aln), "-n", "5", "-o", str(tmp_path / "o"),
                 "--background", str(tmp_path / "missing.fasta")]) == 1
    big = tmp_path / "big.fasta"
    write_fasta(big, [("big", "ACGT" * 2_600_000)])
    assert main(["--aligned", str(aln), "-n", "5", "-o", str(tmp_path / "o2"),
                 "--background", str(big)]) == 1
