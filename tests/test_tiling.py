from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from tests.synthetic import make_genome_family
from vhirprimerlab import structure, tiling
from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.io import alignment_from_records
from vhirprimerlab.models import PipelineError, PrimerCandidate, Record, TilingParams, primer_names
from vhirprimerlab.tiling import choose, geometry, tile

PARAMS = TilingParams(
    n_amplicons=5, overlap=50, tm_min=57, tm_opt=60, tm_max=63, max_tm_diff=2,
    len_min=18, len_max=25, gc_min=40, gc_max=60, max_mismatches=2, window=5,
    degenerate=True, max_degeneracy=4, min_coverage=95,
)


def _prepare(rows):
    aln = alignment_from_records([Record(i, s) for i, s in rows])
    stats = score_columns(aln)
    cons = build_consensus(stats)
    return aln, cons, stats, CoordMap(aln, rows[0][0], cons)


def _params(**kw) -> TilingParams:
    return TilingParams(**{**PARAMS.__dict__, **kw})


@pytest.fixture(scope="module")
def family():
    return _prepare(make_genome_family())


@pytest.fixture(scope="module")
def result(family):
    return tile(*family, PARAMS)


# --- choose() on hand-built candidates -------------------------------------------

def _p(strand: str, start: int, bits: int, n: int = 10) -> PrimerCandidate:
    mask = (1 << bits) - 1 if bits <= n else (1 << n) - 1
    return PrimerCandidate(strand, start, 20, "A" * 20, "A" * 20, "A" * 20, 60.0, None, 50.0, 0.0, mask, [])


def test_choose_maximises_worst_amplicon() -> None:
    # genome 1000, n=2: L=525, lengths 420-630, end slack 131, max outer overlap 209
    geom = geometry(1000, 2, 50, 0.20, 25)
    fwds = [_p("+", 0, 10), _p("+", 350, 7), _p("+", 450, 8)]
    revs = [_p("-", 480, 10), _p("-", 600, 8), _p("-", 960, 10)]
    # path A: f0-r480 (10) then f350-r960 (7): sum 17, worst 7
    # path B: f0-r600 (8)  then f450-r960 (8): sum 16, worst 8  <- expected
    path = choose(fwds, revs, 2, geom, _params(n_amplicons=2), (0, 999))
    assert [(f.start, r.start) for f, r in path] == [(0, 600), (450, 960)]


def test_choose_returns_none_when_no_tiling() -> None:
    geom = geometry(1000, 2, 50, 0.20, 25)
    fwds = [_p("+", 0, 10)]
    revs = [_p("-", 480, 10), _p("-", 960, 10)]
    assert choose(fwds, revs, 2, geom, _params(n_amplicons=2), (0, 999)) is None


# --- tile() on the genome family ---------------------------------------------------

def test_exact_count_pools_and_order(result) -> None:
    amps = result.amplicons
    assert [a.index for a in amps] == [1, 2, 3, 4, 5]
    assert [a.pool for a in amps] == [1, 2, 1, 2, 1]
    starts = [a.fwd.start for a in amps]
    assert starts == sorted(set(starts))


def test_overlap_and_same_pool_spacing(result) -> None:
    amps = result.amplicons
    for a, b in zip(amps, amps[1:]):
        overlap = a.rev.start - b.fwd.end - 1
        assert overlap >= 50 and a.overlap_next == overlap
    assert amps[-1].overlap_next is None
    for a, c in zip(amps, amps[2:]):
        assert c.fwd.start > a.rev.end


def test_reaches_both_ends(result) -> None:
    assert result.amplicons[0].fwd.end + 1 <= 50
    assert 6000 - result.amplicons[-1].rev.start <= 50
    assert all(n <= 50 for n in result.uncovered_ends)


def test_choose_needs_primers_within_the_end_limits() -> None:
    geom = geometry(1000, 2, 50, 0.20, 25)
    fwds = [_p("+", 0, 10), _p("+", 450, 10)]  # first insert starts 20 nt in
    revs = [_p("-", 600, 10), _p("-", 960, 10)]  # last insert ends 40 nt before the end
    assert choose(fwds, revs, 2, geom, _params(n_amplicons=2), (0, 999))
    assert choose(fwds, revs, 2, geom, _params(n_amplicons=2, max_end_loss=39), (0, 999)) is None
    assert choose(fwds, revs, 2, geom, _params(n_amplicons=2, max_start_loss=19), (0, 999)) is None


def test_the_two_end_limits_are_independent() -> None:
    geom = geometry(1000, 2, 50, 0.20, 25)
    fwds = [_p("+", 100, 10), _p("+", 450, 10)]  # first insert starts 120 nt in
    revs = [_p("-", 600, 10), _p("-", 960, 10)]  # last insert ends 40 nt before the end
    assert choose(fwds, revs, 2, geom, _params(n_amplicons=2), (0, 999)) is None
    loose = _params(n_amplicons=2, max_start_loss=120)
    assert choose(fwds, revs, 2, geom, loose, (0, 999))


def test_choose_prefers_primers_that_match_on_their_own() -> None:
    # both last amplicons cover 5 sequences (the reverse misses the rest), but one
    # forward primer matches all 10 sequences on its own
    geom = geometry(1000, 2, 50, 0.20, 25)
    weak, strong = _p("+", 440, 5), _p("+", 450, 10)
    fwds = [_p("+", 0, 10), weak, strong]
    revs = [_p("-", 600, 10), _p("-", 960, 5)]
    path = choose(fwds, revs, 2, geom, _params(n_amplicons=2), (0, 999))
    assert path[1][0] is strong


def test_no_fit_names_the_unreachable_end(family) -> None:
    aln, cons, stats, cmap = family
    with pytest.raises(PipelineError, match="within 0 nt of the genome start"):
        tile(aln, cons, stats, cmap, _params(max_start_loss=0))
    with pytest.raises(PipelineError, match="within 0 nt of the genome end"):
        tile(aln, cons, stats, cmap, _params(max_end_loss=0))


def test_avoids_hypervariable_block(result) -> None:
    for a in result.amplicons:
        for p in (a.fwd, a.rev):
            assert p.end < 3000 or p.start >= 3300
        assert a.coverage.coverage_pct == 100.0 and not a.below_threshold


def test_mask_agrees_with_score_pair(result) -> None:
    for a in result.amplicons:
        assert bin(a.fwd.mask & a.rev.mask).count("1") == len(a.coverage.covered)


def test_n1_single_amplicon() -> None:
    res = tile(*_prepare(make_genome_family(length=1500, n_seqs=6, hyper=(700, 800))),
               _params(n_amplicons=1))
    assert len(res.amplicons) == 1
    assert res.amplicons[0].pool == 1 and res.amplicons[0].overlap_next is None


def test_widening_reported(family, monkeypatch, caplog) -> None:
    real = tiling.choose

    def strict_fails(fwds, revs, n, geom, params, span):
        return None if geom.tolerance == 0.20 else real(fwds, revs, n, geom, params, span)

    monkeypatch.setattr(tiling, "choose", strict_fails)
    with caplog.at_level(logging.INFO, logger="vhirprimerlab"):
        res = tile(*family, PARAMS)
    assert res.geometry.tolerance == 0.35
    assert "±35 %" in caplog.text


def test_no_fit_message() -> None:
    rows = [(i, s[:1400] + "A" * 3200 + s[4600:]) for i, s in make_genome_family()]
    with pytest.raises(PipelineError) as exc:
        tile(*_prepare(rows), _params(n_amplicons=3))
    msg = str(exc.value)
    assert "no acceptable primers between" in msg and "try -n" in msg


def test_n_too_large(family) -> None:
    with pytest.raises(PipelineError, match="use at most"):
        tile(*family, _params(n_amplicons=40))


def test_dimer_warning_format(family, monkeypatch) -> None:
    monkeypatch.setattr(structure.primer3, "calc_heterodimer",
                        lambda a, b: SimpleNamespace(dg=-12000.0, tm=0.0))
    res = tile(*family, PARAMS)
    assert res.dimer_warnings
    assert res.dimer_warnings[0].startswith("pool 1: VHIRPrimerLab_1_LEFT × ")
    assert "-12.0 kcal/mol" in res.dimer_warnings[0]


def test_dimer_warning_quotes_the_amplicon_names(family, monkeypatch) -> None:
    monkeypatch.setattr(structure.primer3, "calc_heterodimer",
                        lambda a, b: SimpleNamespace(dg=-12000.0, tm=0.0))
    monkeypatch.setattr(tiling, "_pair_name", lambda k: f"GENE_{k}")
    res = tile(*family, PARAMS)
    assert res.amplicons[0].name == "GENE_1"
    assert res.dimer_warnings[0].startswith("pool 1: GENE_1_LEFT × GENE_1_RIGHT ")
    known = {p for a in res.amplicons for p in primer_names(a.name)}
    for w in res.dimer_warnings:
        assert set(w.split(": ", 1)[1].split(" ΔG")[0].split(" × ")) <= known


def test_choose_prefers_reaching_the_ends() -> None:
    geom = geometry(1000, 1, 50, 0.20, 25)  # one amplicon, end slack 250
    import dataclasses
    # primers at the very ends have a worse penalty; reaching the ends must still win
    fwds = [dataclasses.replace(_p("+", 0, 10), penalty=1.0), _p("+", 150, 10)]
    revs = [_p("-", 830, 10), dataclasses.replace(_p("-", 980, 10), penalty=1.0)]
    path = choose(fwds, revs, 1, geom, _params(n_amplicons=1), (0, 999))
    assert [(f.start, r.start) for f, r in path] == [(0, 980)]


def test_no_fit_suggests_a_smaller_n_that_works() -> None:
    rows = [(i, s[:2000] + "N" * 1500 + s[3500:]) for i, s in make_genome_family()]
    prepared = _prepare(rows)
    with pytest.raises(PipelineError) as exc:
        tile(*prepared, _params(n_amplicons=10))
    hint = str(exc.value).split("try -n ", 1)[1]
    suggested = [int(x) for x in hint.split(", ")]
    assert suggested and all(m < 10 for m in suggested)
    assert len(tile(*prepared, _params(n_amplicons=suggested[0])).amplicons) == suggested[0]


def test_wide_variable_block_gives_a_flagged_tiling_not_an_error() -> None:
    # 1,200 bp where every base varies (p=0.25): too wide for 10 amplicons to step over
    res = tile(*_prepare(make_genome_family(hyper=(2500, 3700))), _params(n_amplicons=10))
    assert len(res.amplicons) == 10
    low = [a for a in res.amplicons if a.below_threshold]
    assert low
    near = lambda p: p.end >= 2400 and p.start <= 3800
    for a in low:  # only amplicons with a primer in or near the block fall short
        assert near(a.fwd) or near(a.rev)
    for a in res.amplicons:
        if a.rev.end < 2400 or a.fwd.start > 3800:
            assert not a.below_threshold


# --- dimer avoidance on hand-built candidates ---------------------------------------

def _dimer_setup(alt_bits: int):
    import dataclasses
    geom = geometry(1000, 1, 50, 0.20, 25)  # one amplicon, lengths 800-1200, end slack 250
    bad = dataclasses.replace(_p("+", 0, 10), seq="C" * 20, penalty=1.0)
    alt = dataclasses.replace(_p("+", 100, alt_bits), seq="G" * 20)
    rev = dataclasses.replace(_p("-", 980, 10), seq="T" * 20)
    return geom, [bad, alt], [rev]


NAMES = ["VHIRPrimerLab_1"]


def _c_sticks(a: str, b: str) -> float:
    return -12000.0 if "C" in a + b else -3000.0


def test_avoid_dimers_swaps_out_a_dimer_primer(monkeypatch) -> None:
    monkeypatch.setattr(tiling, "dimer_dg", _c_sticks)
    geom, fwds, revs = _dimer_setup(alt_bits=10)
    params = _params(n_amplicons=1, max_start_loss=150)  # the alternative starts at 100
    path, g = tiling._choose_any(fwds, revs, 1, [geom], params, (0, 999), set())
    assert path[0][0].start == 0 and tiling._dimers(path, NAMES)
    path, g = tiling._avoid_dimers(path, NAMES, g, fwds, revs, 1, [geom], params, (0, 999))
    assert path[0][0].start == 100
    assert tiling._dimers(path, NAMES) == []


def test_avoid_dimers_keeps_pair_when_swap_loses_coverage(monkeypatch) -> None:
    monkeypatch.setattr(tiling, "dimer_dg", _c_sticks)
    geom, fwds, revs = _dimer_setup(alt_bits=6)  # the alternative covers fewer sequences
    params = _params(n_amplicons=1)
    path, g = tiling._choose_any(fwds, revs, 1, [geom], params, (0, 999), set())
    path, g = tiling._avoid_dimers(path, NAMES, g, fwds, revs, 1, [geom], params, (0, 999))
    assert path[0][0].start == 0
    assert tiling._dimers(path, NAMES) == [
        "pool 1: VHIRPrimerLab_1_LEFT × VHIRPrimerLab_1_RIGHT ΔG -12.0 kcal/mol"]


def test_tile_within_a_span_stays_inside_it(family) -> None:
    aln, cons, stats, cmap = family
    result = tile(aln, cons, stats, cmap, _params(n_amplicons=2), span=(1000, 2999))
    for a in result.amplicons:
        assert a.fwd.start >= 1000 and a.rev.end <= 2999
    assert len(result.amplicons) == 2


def test_the_span_edges_carry_the_end_rules(family) -> None:
    aln, cons, stats, cmap = family
    result = tile(aln, cons, stats, cmap, _params(n_amplicons=2), span=(1000, 2999))
    assert result.amplicons[0].fwd.end + 1 - 1000 <= 50
    assert 2999 + 1 - result.amplicons[-1].rev.start <= 50


def test_a_span_with_no_primers_names_the_span_edge(family) -> None:
    aln, cons, stats, cmap = family
    with pytest.raises(PipelineError, match="of the span start"):
        tile(aln, cons, stats, cmap, _params(n_amplicons=1, max_start_loss=0), span=(1000, 2999))
