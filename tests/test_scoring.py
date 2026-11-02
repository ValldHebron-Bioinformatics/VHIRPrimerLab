from __future__ import annotations

import pytest

from vhirprimerlab.models import Alignment, ScoringParams
from vhirprimerlab.scoring import match_primer, score_pair

P = ScoringParams(max_mismatches=2, window=5)
SITE = "ACGTTGCAAGCTTCGGATCC"  # 20-mer
COLS = list(range(20))


def mut(seq: str, *positions: int, to: str | None = None) -> str:
    s = list(seq)
    for p in positions:
        s[p] = to if to is not None else {"A": "C", "C": "G", "G": "T", "T": "A"}[s[p]]
    return "".join(s)


def test_exact_match() -> None:
    r = match_primer(SITE, COLS, SITE, "right", P)
    assert r.ok and r.mismatches == 0 and r.reason is None


def test_two_mismatches_outside_window_ok() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 0, 3), "right", P)
    assert r.ok and r.mismatches == 2


def test_three_mismatches_rejected() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 0, 3, 7), "right", P)
    assert not r.ok and r.reason == "too many mismatches"


def test_single_3prime_mismatch_forward_rejected() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 19), "right", P)
    assert not r.ok and r.reason == "3' mismatch" and r.three_prime_mismatches == 1


def test_window_edge_forward() -> None:
    assert not match_primer(SITE, COLS, mut(SITE, 15), "right", P).ok
    assert match_primer(SITE, COLS, mut(SITE, 14), "right", P).ok


def test_single_3prime_mismatch_reverse_rejected() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 0), "left", P)
    assert not r.ok and r.reason == "3' mismatch"


def test_same_mismatch_outside_reverse_window_ok() -> None:
    assert match_primer(SITE, COLS, mut(SITE, 19), "left", P).ok


def test_gap_and_n_are_mismatches() -> None:
    row = mut(mut(SITE, 2, to="-"), 6, to="N")
    r = match_primer(SITE, COLS, row, "right", P)
    assert r.mismatches == 2


def test_insertion_counted() -> None:
    row = SITE[:10] + "G" + SITE[10:]
    cols = list(range(10)) + list(range(11, 21))
    r = match_primer(SITE, cols, row, "right", P)
    assert r.mismatches == 1 and r.ok


def test_gap_in_skipped_column_not_counted() -> None:
    row = SITE[:10] + "-" + SITE[10:]
    cols = list(range(10)) + list(range(11, 21))
    assert match_primer(SITE, cols, row, "right", P).mismatches == 0


def test_insertion_inside_window_rejected() -> None:
    row = SITE[:17] + "G" + SITE[17:]
    cols = list(range(17)) + list(range(18, 21))
    r = match_primer(SITE, cols, row, "right", P)
    assert not r.ok and r.reason == "3' mismatch"


def test_degenerate_site_base_matches() -> None:
    site = SITE[:2] + "R" + SITE[3:]
    assert match_primer(site, COLS, SITE[:2] + "A" + SITE[3:], "right", P).mismatches == 0
    assert match_primer(site, COLS, SITE[:2] + "G" + SITE[3:], "right", P).mismatches == 0
    assert match_primer(site, COLS, SITE[:2] + "C" + SITE[3:], "right", P).mismatches == 1


def test_score_pair_reasons() -> None:
    rev_site = "TTGACCGTAGGCATCGATGC"
    good = SITE + "A" * 20 + rev_site
    bad = mut(SITE, 19) + "A" * 20 + rev_site
    aln = Alignment(["s1", "s2", "s3"], [good, good, bad])
    cov = score_pair(SITE, COLS, rev_site, list(range(40, 60)), aln, P)
    assert cov.covered == ["s1", "s2"]
    assert cov.not_covered == {"s3": "fwd 3' mismatch"}
    assert cov.coverage_pct == pytest.approx(66.67, abs=0.01)


def test_score_pair_both_fail() -> None:
    rev_site = "TTGACCGTAGGCATCGATGC"
    good = SITE + "A" * 20 + rev_site
    bad = mut(SITE, 19) + "A" * 20 + mut(rev_site, 10, 12, 14)
    aln = Alignment(["s1", "s2"], [good, bad])
    cov = score_pair(SITE, COLS, rev_site, list(range(40, 60)), aln, P)
    assert cov.not_covered == {"s2": "fwd 3' mismatch; rev too many mismatches"}
