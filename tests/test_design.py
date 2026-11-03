from __future__ import annotations

from Bio.Seq import reverse_complement

from vhirprimerlab.design import run_primer3
from vhirprimerlab.models import Consensus, DesignParams


def params(amin: int, amax: int, candidates: int = 50) -> DesignParams:
    return DesignParams(
        tm_min=57, tm_opt=60, tm_max=63, max_tm_diff=2,
        len_min=18, len_max=25, gc_min=40, gc_max=60,
        amplicon_min=amin, amplicon_max=amax, candidates=candidates,
    )


def test_pairs_inside_included_region(family_base: str) -> None:
    cons = Consensus(family_base, list(range(len(family_base))))
    r = run_primer3(cons, (100, 400), params(100, 300))
    assert r.pairs
    for p in r.pairs:
        assert p.forward.start >= 100
        assert p.reverse.start + p.reverse.length <= 500
        amplicon = p.reverse.start + p.reverse.length - p.forward.start
        assert 100 <= amplicon <= 300


def test_primer_sequences_match_consensus(family_base: str) -> None:
    cons = Consensus(family_base, list(range(len(family_base))))
    p = run_primer3(cons, (100, 400), params(100, 300)).pairs[0]
    assert family_base[p.forward.start: p.forward.start + p.forward.length] == p.forward.seq
    assert family_base[p.reverse.start: p.reverse.start + p.reverse.length] == reverse_complement(
        p.reverse.seq
    )
    assert 57 <= p.forward.tm <= 63 and 40 <= p.forward.gc <= 60


def test_penalties_ascending(family_base: str) -> None:
    cons = Consensus(family_base, list(range(len(family_base))))
    pens = [p.penalty for p in run_primer3(cons, (100, 400), params(100, 300)).pairs]
    assert pens == sorted(pens)


def test_no_candidates_returns_explain() -> None:
    cons = Consensus("A" * 300, list(range(300)))
    r = run_primer3(cons, (0, 300), params(100, 200))
    assert r.pairs == []
    assert "PRIMER_LEFT_EXPLAIN" in r.explain


def test_primers_keep_the_3_prime_gc_limit(family_base: str) -> None:
    from vhirprimerlab.structure import MAX_END_GC, end_gc

    cons = Consensus(family_base, list(range(len(family_base))))
    pairs = run_primer3(cons, (0, len(family_base)), params(100, 300)).pairs
    assert pairs
    for p in pairs:
        assert end_gc(p.forward.seq) <= MAX_END_GC and end_gc(p.reverse.seq) <= MAX_END_GC
