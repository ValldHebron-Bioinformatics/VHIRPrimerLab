from __future__ import annotations

import pytest

from vhirprimerlab.lengths import filter_by_length
from vhirprimerlab.models import PipelineError, Record


def _recs(extra: list[Record] | None = None) -> list[Record]:
    base = [Record(f"s{i}", "A" * 100) for i in range(1, 5)]
    return base + (extra if extra is not None else [Record("long", "A" * 150)])


def test_drop_removes_outlier() -> None:
    r = filter_by_length(_recs(), 0.10, "drop", "s1")
    assert [x.id for x in r.kept] == ["s1", "s2", "s3", "s4"]
    assert r.outliers == [("long", 150)]
    assert r.median == 100


def test_warn_keeps_outlier() -> None:
    r = filter_by_length(_recs(), 0.10, "warn", "s1")
    assert len(r.kept) == 5
    assert r.outliers == [("long", 150)]


def test_abort_raises() -> None:
    with pytest.raises(PipelineError, match="long"):
        filter_by_length(_recs(), 0.10, "abort", "s1")


def test_boundary_is_inclusive() -> None:
    r = filter_by_length(_recs([Record("edge", "A" * 110)]), 0.10, "drop", "s1")
    assert r.outliers == []
    assert len(r.kept) == 5


def test_gaps_not_counted() -> None:
    r = filter_by_length(_recs([Record("gappy", "A" * 100 + "-" * 50)]), 0.10, "drop", "s1")
    assert r.outliers == []


def test_reference_outlier_dropped_raises() -> None:
    with pytest.raises(PipelineError, match="reference"):
        filter_by_length(_recs(), 0.10, "drop", "long")


def test_fewer_than_two_left_raises() -> None:
    recs = [Record("a", "A" * 100), Record("b", "A" * 300), Record("c", "A" * 10)]
    with pytest.raises(PipelineError, match="2"):
        filter_by_length(recs, 0.10, "drop", "a")
