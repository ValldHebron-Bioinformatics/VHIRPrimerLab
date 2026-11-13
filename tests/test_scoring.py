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


def test_gap_is_a_mismatch_and_n_is_unknown() -> None:
    row = mut(mut(SITE, 2, to="-"), 6, to="N")
    r = match_primer(SITE, COLS, row, "right", P)
    assert r.mismatches == 1 and r.ambiguous == 1 and r.ok


def test_n_outside_the_window_is_not_a_mismatch() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 3, to="N"), "right", P)
    assert r.ok and r.mismatches == 0 and r.ambiguous == 1


def test_a_code_that_includes_the_primer_base_is_unknown_one_that_excludes_it_a_mismatch() -> None:
    assert SITE[0] == "A"
    compatible = match_primer(SITE, COLS, mut(SITE, 0, to="R"), "right", P)   # R = A or G
    assert compatible.mismatches == 0 and compatible.ambiguous == 1
    excluded = match_primer(SITE, COLS, mut(SITE, 0, to="Y"), "right", P)    # Y = C or T
    assert excluded.mismatches == 1 and excluded.ambiguous == 0


def test_a_code_inside_a_degenerate_primer_base_is_a_plain_match() -> None:
    site = "R" + SITE[1:]                                                    # primer base A or G
    row = mut(SITE, 0, to="R")
    r = match_primer(site, COLS, row, "right", P)
    assert r.mismatches == 0 and r.ambiguous == 0


def test_unknown_bases_past_the_allowance_count_as_mismatches() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 0, 2, 4, to="N"), "right", P)      # allowance is 2
    assert r.ambiguous == 2 and r.mismatches == 1 and r.ok


def test_a_run_of_n_still_fails_the_primer() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 0, 1, 2, 3, 4, 5, to="N"), "right", P)
    assert not r.ok and r.reason == "too many mismatches"


def test_an_unknown_base_at_the_3_prime_end_fails_and_says_why() -> None:
    r = match_primer(SITE, COLS, mut(SITE, 19, to="N"), "right", P)
    assert not r.ok and r.reason == "ambiguous 3' base"
    left = match_primer(SITE, COLS, mut(SITE, 0, to="N"), "left", P)
    assert not left.ok and left.reason == "ambiguous 3' base"


def test_a_real_3_prime_mismatch_keeps_its_own_reason() -> None:
    row = mut(mut(SITE, 19), 18, to="N")
    assert match_primer(SITE, COLS, row, "right", P).reason == "3' mismatch"


def test_max_ambiguous_zero_counts_every_unknown_base_as_a_mismatch() -> None:
    strict = ScoringParams(max_mismatches=2, window=5, max_ambiguous=0)
    r = match_primer(SITE, COLS, mut(SITE, 3, to="N"), "right", strict)
    assert r.mismatches == 1 and r.ambiguous == 0


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


def test_a_sequence_with_an_n_inside_a_primer_site_is_still_covered() -> None:
    rev_site = "TTGACCGTAGGCATCGATGC"
    good = SITE + "A" * 20 + rev_site
    draft = mut(SITE, 7, to="N") + "A" * 20 + mut(rev_site, 9, to="N")   # a draft assembly
    ends_in_n = SITE[:19] + "N" + "A" * 20 + rev_site
    aln = Alignment(["s1", "s2", "s3"], [good, draft, ends_in_n])
    cov = score_pair(SITE, COLS, rev_site, list(range(40, 60)), aln, P)
    assert cov.covered == ["s1", "s2"]
    assert cov.not_covered == {"s3": "fwd ambiguous 3' base"}
