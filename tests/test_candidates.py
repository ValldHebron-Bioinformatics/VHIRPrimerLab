from __future__ import annotations

from collections import Counter

import pytest
from Bio.Seq import reverse_complement

from tests.synthetic import make_genome_family
from viraprime.candidates import enumerate_primers, pick_conserved, score_primers
from viraprime.conservation import build_consensus, score_columns
from viraprime.degenerate import degenerate_site, fold
from viraprime.io import alignment_from_records
from viraprime.models import Record, ScoringParams, TilingParams
from viraprime.scoring import match_primer

PARAMS = TilingParams(
    n_amplicons=5, overlap=50, tm_min=57, tm_opt=60, tm_max=63, max_tm_diff=2,
    len_min=18, len_max=25, gc_min=40, gc_max=60, max_mismatches=2, window=5,
    degenerate=True, max_degeneracy=4, min_coverage=95,
)


@pytest.fixture(scope="module")
def family():
    aln = alignment_from_records([Record(i, s) for i, s in make_genome_family()])
    stats = score_columns(aln)
    return aln, stats, build_consensus(stats)


@pytest.fixture(scope="module")
def raw(family):
    return enumerate_primers(family[2].seq, PARAMS)


def test_enumerate_respects_limits(raw) -> None:
    fwds, revs = raw
    assert len(fwds) > 100 and len(revs) > 100
    for c in fwds + revs:
        assert 18 <= c.length <= 25
        assert 57 <= c.tm <= 63
        assert 40 <= c.gc <= 60
    assert all(c.strand == "+" and c.seq == c.site == c.base_seq for c in fwds)
    assert all(c.strand == "-" and c.seq == reverse_complement(c.site) for c in revs)


def test_enumerate_skips_non_acgt(family) -> None:
    seq = family[2].seq[:1000]
    seq = seq[:100] + "N" + seq[101:]
    fwds, revs = enumerate_primers(seq, PARAMS)
    assert not [c for c in fwds + revs if c.start <= 100 <= c.end]


def test_pick_conserved_keeps_the_best_of_every_window(family, raw) -> None:
    """Variable regions keep their least-variable primers; the tiling search decides."""
    _, stats, cons = family
    for cands in raw:
        picked = pick_conserved(cands, stats, cons, PARAMS)
        anchor = (lambda c: c.end) if cands[0].strand == "+" else (lambda c: c.start)
        assert max(Counter(anchor(c) // 40 for c in picked).values()) <= 3
        # every 40-bp window that had candidates still has some, inside the variable block too
        assert {anchor(c) // 40 for c in picked} == {anchor(c) // 40 for c in cands}
        assert [c for c in picked if c.start >= 3000 and c.end < 3300]


def test_score_matches_match_primer(family, raw) -> None:
    aln, stats, cons = family
    plain = TilingParams(**{**PARAMS.__dict__, "degenerate": False})
    picked = pick_conserved(raw[0], stats, cons, plain)[:10] + pick_conserved(raw[1], stats, cons, plain)[:10]
    for c in score_primers(picked, aln, cons, plain):
        assert c.cols == cons.cons_to_col[c.start:c.start + c.length]
        side = "right" if c.strand == "+" else "left"
        for i, row in enumerate(aln.rows):
            ok = match_primer(c.site, c.cols, row, side, plain.scoring).ok
            assert bool(c.mask >> i & 1) == ok


def _aln(rows: list[str]):
    return alignment_from_records([Record(f"s{i}", r) for i, r in enumerate(rows)])


def test_degenerate_site_rescues_and_keeps_3prime() -> None:
    site = "ACGTACGTACATGCATGCAT"  # position 10 is A
    variant = site[:10] + "G" + site[11:]
    aln = _aln([site] * 5 + [variant] * 5)
    cols = list(range(20))
    strict = ScoringParams(0, 5)
    out = degenerate_site(site, cols, aln, strict, 4, "right")
    assert out == site[:10] + "R" + site[11:] and fold(out) == 2
    assert all(match_primer(out, cols, r, "right", strict).ok for r in aln.rows)

    last = site[:19] + "G"  # variable 3'-terminal base: never made degenerate
    aln2 = _aln([site] * 5 + [last] * 5)
    assert degenerate_site(site, cols, aln2, strict, 4, "right") == site


def test_degenerate_site_respects_max_fold() -> None:
    site = "ACGTACGTACATGCATGCAT"
    aln = _aln([site] * 5 + [site[:10] + "G" + site[11:]] * 5)
    assert degenerate_site(site, list(range(20)), aln, ScoringParams(0, 5), 1, "right") == site


def test_score_primers_degenerate_rescues(family, raw) -> None:
    aln, stats, cons = family
    picked = pick_conserved(raw[0], stats, cons, PARAMS)
    plain = score_primers(picked, aln, cons, TilingParams(**{**PARAMS.__dict__, "degenerate": False}))
    degen = score_primers(picked, aln, cons, PARAMS)
    full = (1 << len(aln.ids)) - 1
    assert sum(c.mask == full for c in degen) >= sum(c.mask == full for c in plain)
    for c in degen:
        assert fold(c.site) <= 4
        if c.site != c.base_seq:
            assert c.tm_range is not None


def test_degenerate_site_keeps_gc_within_limits() -> None:
    # 8 G/C of 20 = 40 % GC; making position 10 (G) into R (A/G) drops the average to 37.5 %
    site = "AATTGCAGCTGCATGCATAT"
    variant = site[:10] + "A" + site[11:]
    aln = _aln([site] * 5 + [variant] * 5)
    strict = ScoringParams(0, 5)
    assert degenerate_site(site, list(range(20)), aln, strict, 4, "right") != site
    assert degenerate_site(site, list(range(20)), aln, strict, 4, "right", gc_limits=(40, 60)) == site
