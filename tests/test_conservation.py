from __future__ import annotations

import pytest

from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.models import Alignment

ALN = Alignment(["ref", "s2", "s3"], ["AC--GTA", "ACTTGTA", "AC--GAA"])


def test_gap_majority_columns_dropped() -> None:
    cons = build_consensus(score_columns(ALN))
    assert cons.seq == "ACGTA"
    assert cons.cons_to_col == [0, 1, 4, 5, 6]


def test_column_stats_values() -> None:
    stats = score_columns(ALN)
    assert stats[5].majority_base == "T"
    assert stats[5].majority_freq == pytest.approx(2 / 3)
    assert stats[2].gap_frac == pytest.approx(2 / 3)


def test_majority_tie_order() -> None:
    stats = score_columns(Alignment(["a", "b"], ["C", "A"]))
    assert stats[0].majority_base == "A"


def test_ambiguous_counted() -> None:
    stats = score_columns(Alignment(["a", "b", "c"], ["A", "N", "A"]))
    assert stats[0].n_ambiguous == 1
    assert stats[0].majority_freq == pytest.approx(2 / 3)
