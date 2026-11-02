from __future__ import annotations

import csv
from pathlib import Path

from vhirprimerlab.cli import main
from tests.conftest import needs_mafft
from tests.synthetic import make_family, write_fasta

AMP = ["--amplicon-min", "100", "--amplicon-max", "300"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@needs_mafft
def test_hypervariable_region_exit_2(family: Path, tmp_path: Path) -> None:
    out = tmp_path / "o"
    rc = main(["-i", str(family), "-o", str(out), "--region", "hv:701-1100", *AMP])
    assert rc == 2
    rows = read_tsv(out / "primers.tsv")
    assert rows, "best available pairs must still be reported"
    assert all(r["below_threshold"] == "True" for r in rows)
    assert all(r["region"] == "hv" for r in rows)
    report = (out / "report.txt").read_text()
    assert "hv" in report and "best coverage" in report.lower()


@needs_mafft
def test_zero_candidates_region_explained(tmp_path: Path) -> None:
    _, recs = make_family()
    recs = [(rid, "A" * 400 + s[400:]) for rid, s in recs]
    fa = write_fasta(tmp_path / "polya.fa", recs)
    out = tmp_path / "o"
    rc = main(["-i", str(fa), "-o", str(out), "--region", "polyA:1-400", "--region", "ok:401-600",
               *AMP])
    assert rc == 2
    rows = read_tsv(out / "primers.tsv")
    assert rows and {r["region"] for r in rows} == {"ok"}
    report = (out / "report.txt").read_text()
    section = report.split("Region polyA")[1].split("Region ok")[0]
    assert "no candidate" in section and "considered" in section
