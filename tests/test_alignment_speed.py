"""alignment_from_records drops all-gap columns; the fast version must agree with the plain one."""

from __future__ import annotations

import random

import pytest

from vhirprimerlab.io import alignment_from_records
from vhirprimerlab.models import Alignment, PipelineError, Record


def reference(records: list[Record]) -> Alignment:
    width = len(records[0].seq)
    keep = [c for c in range(width) if any(r.seq[c] != "-" for r in records)]
    return Alignment(ids=[r.id for r in records], rows=["".join(r.seq[c] for c in keep) for r in records])


@pytest.mark.parametrize("seed", range(60))
def test_same_columns_as_the_plain_version(seed: int) -> None:
    rng = random.Random(seed)
    width = rng.randint(1, 120)
    gappy = rng.choice([0.0, 0.3, 0.9, 1.0])      # how many columns are gap-only in every row
    dead = {c for c in range(width) if rng.random() < gappy}
    records = []
    for i in range(rng.randint(1, 12)):
        seq = "".join("-" if c in dead or rng.random() < 0.1 else rng.choice("ACGTNRYK") for c in range(width))
        records.append(Record(f"s{i}", seq))
    got = alignment_from_records(records)
    want = reference(records)
    assert got.ids == want.ids and got.rows == want.rows


def test_rows_of_different_length_are_still_refused() -> None:
    with pytest.raises(PipelineError, match="differ in length"):
        alignment_from_records([Record("a", "ACGT"), Record("b", "ACG")])


def test_lower_case_and_unusual_characters_survive() -> None:
    rows = alignment_from_records([Record("a", "AC-t"), Record("b", "A--t")]).rows
    assert rows == ["ACt", "A-t"]          # the column of two gaps goes, the lower-case t stays


def test_a_column_that_is_all_gaps_in_every_row_is_dropped() -> None:
    rows = alignment_from_records([Record("a", "A-C"), Record("b", "G-T")]).rows
    assert rows == ["AC", "GT"]
