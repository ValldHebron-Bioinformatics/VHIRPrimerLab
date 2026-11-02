from __future__ import annotations

from pathlib import Path

import pytest

from vhirprimerlab.io import load_fasta
from vhirprimerlab.models import PipelineError, Record
from tests.synthetic import write_fasta


def test_empty_file_raises(tmp_path: Path) -> None:
    (tmp_path / "e.fa").write_text("")
    with pytest.raises(PipelineError, match="empty"):
        load_fasta(tmp_path / "e.fa")


def test_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="not found"):
        load_fasta(tmp_path / "nope.fa")


def test_duplicate_ids_listed(tmp_path: Path) -> None:
    p = write_fasta(tmp_path / "d.fa", [("a", "ACGT"), ("b", "ACGT"), ("a", "ACGT"), ("b", "ACGT")])
    with pytest.raises(PipelineError, match="Duplicate IDs: a, b"):
        load_fasta(p)


def test_invalid_char_reports_id_and_position(tmp_path: Path) -> None:
    p = write_fasta(tmp_path / "x.fa", [("s1", "ACGTXA")])
    with pytest.raises(PipelineError, match=r"s1.*position 5.*'X'"):
        load_fasta(p)


def test_gap_rejected_unless_allowed(tmp_path: Path) -> None:
    p = write_fasta(tmp_path / "g.fa", [("s1", "AC-GT")])
    with pytest.raises(PipelineError):
        load_fasta(p)
    assert load_fasta(p, allow_gaps=True) == [Record("s1", "AC-GT")]


def test_empty_record_raises(tmp_path: Path) -> None:
    p = tmp_path / "r.fa"
    p.write_text(">s1\nACGT\n>s2\n\n")
    with pytest.raises(PipelineError, match="s2.*empty"):
        load_fasta(p)


def test_normalises_case_u_and_crlf(tmp_path: Path) -> None:
    p = tmp_path / "c.fa"
    p.write_bytes(b">s1\r\nacgu\r\n")
    assert load_fasta(p) == [Record("s1", "ACGT")]


def test_family_fixture_loads(family: Path) -> None:
    assert len(load_fasta(family)) == 10
