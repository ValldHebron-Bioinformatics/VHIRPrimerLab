from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from vhirprimerlab import align
from vhirprimerlab.align import check_binary, run_aligner
from vhirprimerlab.io import alignment_from_records, load_alignment, write_alignment
from vhirprimerlab.models import Alignment, PipelineError, Record
from tests.conftest import needs_mafft
from tests.synthetic import write_fasta


def test_missing_binary_hint(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(align.shutil, "which", lambda _: None)
    with pytest.raises(PipelineError, match="conda install -c bioconda mafft"):
        check_binary("mafft")


@needs_mafft
def test_mafft_aligns_indel() -> None:
    recs = [Record("a", "ACGTACGTTTGCAACGTGGA"), Record("b", "ACGTACTTGCAACGTGGA")]
    aln = run_aligner(recs, "mafft")
    assert aln.ids == ["a", "b"]
    assert len(aln.rows[0]) == len(aln.rows[1])
    assert "-" in aln.row("b")
    assert aln.row("a").replace("-", "") == "ACGTACGTTTGCAACGTGGA"


@pytest.mark.skipif(shutil.which("muscle") is None, reason="muscle not installed")
def test_muscle_aligns() -> None:
    recs = [Record("a", "ACGTACGTTTGCAACGTGGA"), Record("b", "ACGTACTTGCAACGTGGA")]
    aln = run_aligner(recs, "muscle")
    assert aln.ids == ["a", "b"]
    assert len(aln.rows[0]) == len(aln.rows[1])


def test_load_alignment_unequal_lengths(tmp_path: Path) -> None:
    p = write_fasta(tmp_path / "a.fa", [("a", "AC-GT"), ("b", "ACGT")])
    with pytest.raises(PipelineError, match="b"):
        load_alignment(p)


def test_alignment_from_records_drops_all_gap_columns() -> None:
    aln = alignment_from_records([Record("x", "A-C"), Record("y", "A-G")])
    assert aln.rows == ["AC", "AG"]


def test_write_alignment_roundtrip(tmp_path: Path) -> None:
    aln = Alignment(["x", "y"], ["AC-T", "ACGT"])
    write_alignment(aln, tmp_path / "out.fa")
    assert load_alignment(tmp_path / "out.fa") == aln
