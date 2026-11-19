"""The README's column lists must be the ones the writers emit.

A reader building a parser from the README breaks when a column is added to the code and not to the
docs. This reads each list out of the README and compares it with the constant the writer uses.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from vhirprimerlab import io

README = (Path(__file__).resolve().parent.parent / "README.md").read_text()


def _listed(after: str) -> list[str]:
    """The comma-separated names in the first backtick span after `after`, newlines folded."""
    match = re.search(re.escape(after) + r"[^`]*`([^`]+)`", README)
    assert match, f"README has no column list after {after!r}"
    return [name.strip() for name in " ".join(match.group(1).split()).split(",") if name.strip()]


@pytest.mark.parametrize("marker, columns", [
    ("`pairs.tsv` — one row per reported pair", io.TSV_COLUMNS),
])
def test_the_readme_lists_the_columns_the_writer_emits(marker: str, columns: list[str]) -> None:
    assert _listed(marker) == list(columns)
