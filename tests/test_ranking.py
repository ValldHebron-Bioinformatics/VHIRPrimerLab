from __future__ import annotations

from vhirprimerlab.models import CandidatePair, PairCoverage, Primer, ScoredPair
from vhirprimerlab.ranking import failure_summary, rank


def sp(cov: float, pen: float, not_covered: dict[str, str] | None = None) -> ScoredPair:
    p = Primer("ACGT", 0, 4, 60.0, 50.0)
    return ScoredPair(
        pair=CandidatePair(p, p, pen), fwd_seq="ACGT", rev_seq="ACGT", fwd_cols=[0], rev_cols=[0],
        coverage=PairCoverage([], not_covered or {}, cov), degenerate=False,
        fwd_tm_range=None, rev_tm_range=None,
    )


def test_coverage_then_penalty() -> None:
    a, b, c = sp(100, 2), sp(100, 1), sp(96, 0)
    assert rank([a, b, c], 95, 5).pairs == [b, a, c]


def test_filters_below_min() -> None:
    good, bad = sp(100, 1), sp(90, 0)
    r = rank([bad, good], 95, 5)
    assert r.pairs == [good] and not r.below_threshold


def test_none_pass_returns_best_flagged() -> None:
    pairs = [sp(50, 0), sp(80, 3), sp(70, 1)]
    r = rank(pairs, 95, 2)
    assert r.below_threshold
    assert [p.coverage.coverage_pct for p in r.pairs] == [80, 70]


def test_truncates_to_num_pairs() -> None:
    assert len(rank([sp(100, i) for i in range(10)], 95, 5).pairs) == 5


def test_empty_input() -> None:
    r = rank([], 95, 5)
    assert r.pairs == [] and r.below_threshold


def test_failure_summary_counts() -> None:
    best = sp(70, 0, {"a": "fwd 3' mismatch", "b": "fwd 3' mismatch", "c": "rev too many mismatches"})
    assert failure_summary([sp(10, 0, {"x": "rev 3' mismatch"}), best]) == [
        ("fwd 3' mismatch", 2), ("rev too many mismatches", 1),
    ]


def test_failure_summary_splits_combined_reasons() -> None:
    best = sp(50, 0, {"a": "fwd 3' mismatch; rev 3' mismatch"})
    assert failure_summary([best]) == [("fwd 3' mismatch", 1), ("rev 3' mismatch", 1)]
