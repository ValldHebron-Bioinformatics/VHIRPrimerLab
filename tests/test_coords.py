from __future__ import annotations

import pytest

from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.models import Alignment, PipelineError, Region


def _cmap(aln: Alignment, ref: str = "ref") -> CoordMap:
    return CoordMap(aln, ref, build_consensus(score_columns(aln)))


ALN = Alignment(["ref", "s2", "s3"], ["AC--GTA", "ACTTGTA", "AC--GAA"])


def test_ref_length() -> None:
    assert _cmap(ALN).ref_length == 5


def test_ref_to_col_skips_gaps() -> None:
    cm = _cmap(ALN)
    assert cm.ref_to_col(1) == 0
    assert cm.ref_to_col(3) == 4  # 1-based ref pos 3 = 'G' at column 4


def test_col_to_ref_snaps_on_gap() -> None:
    cm = _cmap(ALN)
    assert cm.col_to_ref(2, "start") == 3
    assert cm.col_to_ref(3, "end") == 2


def test_col_to_cons_snaps_on_dropped_column() -> None:
    cm = _cmap(ALN)
    assert cm.col_to_cons(2, "start") == 2
    assert cm.col_to_cons(3, "end") == 1


def test_roundtrip_ref_col_ref() -> None:
    cm = _cmap(ALN)
    for pos in range(1, 6):
        assert cm.col_to_ref(cm.ref_to_col(pos), "start") == pos


def test_region_to_included() -> None:
    assert _cmap(ALN).region_to_included(Region("r", 2, 4)) == (1, 3)


def test_region_outside_reference_raises() -> None:
    with pytest.raises(PipelineError, match="outside"):
        _cmap(ALN).validate_region(Region("r", 1, 99), amplicon_min=2)


def test_region_shorter_than_amplicon_min_raises() -> None:
    with pytest.raises(PipelineError, match="amplicon-min"):
        _cmap(ALN).validate_region(Region("r", 1, 4), amplicon_min=10)


def test_unknown_reference_raises() -> None:
    with pytest.raises(PipelineError, match="nope"):
        _cmap(ALN, ref="nope")


def test_ref_gap_at_primer_edge() -> None:
    # reference lacks the base in column 2, which the consensus keeps;
    # a primer starting in column 2 must snap to the next reference base
    aln = Alignment(["ref", "s2", "s3"], ["AC-GTAC", "ACTGTAC", "ACTGTAC"])
    cm = _cmap(aln)
    assert build_consensus(score_columns(aln)).cons_to_col[2] == 2
    start = cm.col_to_ref(2, "start")
    assert start == 3
    assert 1 <= start <= cm.ref_length
