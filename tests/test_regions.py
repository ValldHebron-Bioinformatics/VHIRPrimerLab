from __future__ import annotations

from pathlib import Path

import pytest

from vhirprimerlab.models import PipelineError, Region
from vhirprimerlab.regions import load_regions_file, parse_region_spec, resolve_regions


def test_plain_spec() -> None:
    assert parse_region_spec("120-480", 1) == Region("region1", 120, 480)


def test_named_spec() -> None:
    assert parse_region_spec("E1:900-1350", 2) == Region("E1", 900, 1350)


@pytest.mark.parametrize("bad", ["abc", "10-", "5-5", "9-3", "x y:1-9", "0-10", ":1-9"])
def test_bad_spec(bad: str) -> None:
    with pytest.raises(PipelineError):
        parse_region_spec(bad, 1)


def test_file_with_header_and_comments(tmp_path: Path) -> None:
    p = tmp_path / "r.tsv"
    p.write_text("name\tstart\tend\n# c\n\nA\t1\t50\nB\t40\t90\n")
    assert load_regions_file(p) == [Region("A", 1, 50), Region("B", 40, 90)]


def test_file_without_header(tmp_path: Path) -> None:
    p = tmp_path / "r.tsv"
    p.write_text("A\t1\t50\n")
    assert load_regions_file(p) == [Region("A", 1, 50)]


def test_file_bad_line_reports_line_number(tmp_path: Path) -> None:
    p = tmp_path / "r.tsv"
    p.write_text("A\t1\t50\nB\t1\n")
    with pytest.raises(PipelineError, match="line 2"):
        load_regions_file(p)


def test_duplicate_names_raise() -> None:
    with pytest.raises(PipelineError, match="A"):
        resolve_regions(["A:1-10", "A:5-20"], None, 100)


def test_both_sources_raise(tmp_path: Path) -> None:
    p = tmp_path / "r.tsv"
    p.write_text("A\t1\t50\n")
    with pytest.raises(PipelineError, match="--regions-file"):
        resolve_regions(["1-10"], p, 100)


def test_default_full_region() -> None:
    assert resolve_regions(None, None, 600) == [Region("full", 1, 600)]


def test_auto_names_follow_order() -> None:
    assert [r.name for r in resolve_regions(["1-10", "X:5-20", "30-40"], None, 100)] == [
        "region1", "X", "region3",
    ]
