"""The edge advice in report.txt, on the command line, and in advice.json."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pytest

from tests.synthetic import make_genome_family, write_fasta
from tests.test_edge_advice import START_POOR, start_poor  # noqa: F401  (fixture)
from tests.test_tiling import _params
from vhirprimerlab.cli import main
from vhirprimerlab.io import write_tiling_report
from vhirprimerlab.models import EdgeAdvice, LengthFilterResult, TilingReport
from vhirprimerlab.tiling import tile

ADVICE = EdgeAdvice(span="F", end="start", limit=176, extra_loss=76,
                    edge_before=43.0, edge_after=94.6)


def _report(result) -> TilingReport:
    lf = LengthFilterResult(kept=[], outliers=[], median=6000.0, action="drop")
    return TilingReport(
        command="x", params={}, n_input=12, median_len=6000.0, min_len=6000, max_len=6000,
        length=lf, reference_id="g01", aln_length=6000, cons_length=6000, n_conserved=1,
        mean_conservation=1.0, min_coverage=95.0, result=result, log_lines=[])


def test_report_section(start_poor, tmp_path: Path) -> None:
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    res = dataclasses.replace(res, edge_advice=[ADVICE])
    write_tiling_report(_report(res), start_poor[3], tmp_path / "r.txt")
    text = (tmp_path / "r.txt").read_text()
    assert "Edge advice" in text
    line = next(ln for ln in text.splitlines() if "--max-start-loss 176" in ln)
    assert line.strip() == ("The first amplicon of F covers 43.0 % of your sequences; giving up 76 "
                            "more bases at its start would reach 94.6 %: add --max-start-loss 176")
    assert text.index("Edge advice") < text.index("Summary")


def test_report_without_advice_has_no_section(start_poor, tmp_path: Path) -> None:
    res = tile(*start_poor, _params(n_amplicons=5, max_start_loss=300, max_end_loss=60))
    write_tiling_report(_report(res), start_poor[3], tmp_path / "r.txt")
    assert "Edge advice" not in (tmp_path / "r.txt").read_text()


@pytest.fixture(scope="module")
def poor_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("adv")
    aln = write_fasta(d / "aln.fasta", make_genome_family(**START_POOR))
    return aln, d


def test_cli_warns_and_writes_advice_json(poor_run, capsys) -> None:
    aln, d = poor_run
    out = d / "out"
    main(["--aligned", str(aln), "-n", "5", "--max-end-loss", "60", "-o", str(out)])
    err = capsys.readouterr().err
    data = json.loads((out / "advice.json").read_text())
    [item] = data
    assert set(item) == {f.name for f in dataclasses.fields(EdgeAdvice)}
    assert item["span"] is None and item["end"] == "start"
    assert item["limit"] == 50 + item["extra_loss"]
    msgs = [ln for ln in err.splitlines() if "--max-start-loss" in ln]
    assert len(msgs) == 1 and msgs[0].startswith("WARNING: ")
    assert msgs[0].endswith(f"--max-start-loss {item['limit']}")
    assert "Edge advice" in (out / "report.txt").read_text()


def test_cli_without_advice_writes_no_advice_json(tmp_path: Path) -> None:
    aln = write_fasta(tmp_path / "aln.fasta", make_genome_family())
    main(["--aligned", str(aln), "-n", "5", "-o", str(tmp_path / "out")])
    assert not (tmp_path / "out" / "advice.json").exists()


def test_cli_removes_a_stale_advice_json(tmp_path: Path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    (out / "advice.json").write_text('[{"stale": true}]')
    aln = write_fasta(tmp_path / "aln.fasta", make_genome_family())
    main(["--aligned", str(aln), "-n", "5", "-o", str(out)])
    assert not (out / "advice.json").exists()
