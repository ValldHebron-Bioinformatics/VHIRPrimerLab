from __future__ import annotations

import json
import shutil
from pathlib import Path

from tests.synthetic import make_genome_family, write_fasta
from vhirprimerlab import align
from vhirprimerlab.align import register_aligner
from vhirprimerlab.cli import main, set_progress_reporter
from vhirprimerlab.io import alignment_from_records


def test_registered_aligner_is_used(tmp_path: Path, monkeypatch) -> None:
    fasta = write_fasta(tmp_path / "in.fasta", make_genome_family(length=3000))
    calls = []

    def fake(records):
        calls.append(len(records))
        return alignment_from_records(records)

    monkeypatch.setattr(shutil, "which", lambda name: None)  # no MAFFT program
    register_aligner("mafft", fake)
    try:
        assert main(["-i", str(fasta), "-n", "2", "-o", str(tmp_path / "out")]) == 0
    finally:
        register_aligner("mafft", None)
    assert calls == [12]
    assert (tmp_path / "out" / "primers.tsv").is_file()


def test_registered_aligner_missing_ids_error(tmp_path: Path, capsys) -> None:
    fasta = write_fasta(tmp_path / "in.fasta", make_genome_family(length=3000))
    register_aligner("mafft", lambda records: alignment_from_records(records[:-1]))
    try:
        assert main(["-i", str(fasta), "-n", "2", "-o", str(tmp_path / "out")]) == 1
    finally:
        register_aligner("mafft", None)
    assert "g12" in capsys.readouterr().err


def test_unregistered_falls_back_to_program() -> None:
    assert "mafft" not in align._REGISTERED


def test_progress_reporter_receives_stages(tmp_path: Path) -> None:
    aln = write_fasta(tmp_path / "aln.fasta", make_genome_family(length=3000))
    seen: list[str] = []
    set_progress_reporter(lambda stage, fraction: seen.append(stage))
    try:
        rc = main(["--aligned", str(aln), "-n", "2", "-o", str(tmp_path / "out"),
                   "--progress-file", str(tmp_path / "p.json")])
    finally:
        set_progress_reporter(None)
    assert rc == 0
    order = [s for i, s in enumerate(seen) if i == 0 or seen[i - 1] != s]
    for stage in ("aligning", "finding primers", "scoring", "choosing amplicons", "done"):
        assert stage in order
    assert order.index("aligning") < order.index("scoring") < order.index("done")
    assert json.loads((tmp_path / "p.json").read_text())["stage"] == "done"
