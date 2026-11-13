from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from tests.synthetic import make_genome_family
from vhirprimerlab import structure, tiling
from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.coords import CoordMap
from vhirprimerlab.io import alignment_from_records
from vhirprimerlab.models import PipelineError, PrimerCandidate, Record, Region, TilingParams, primer_names
from vhirprimerlab.tiling import choose, geometry, tile, tile_spans

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


# --- tile_spans(): several spans, one scheme ---------------------------------------

def _span_params() -> TilingParams:
    # the fixture's sequence leaves 75 bases before any primer site at the SH start (only the start limit needs loosening)
    return _params(amplicon_len=700, max_start_loss=100)


def test_tile_spans_numbers_and_names_across_spans(family) -> None:
    aln, cons, stats, cmap = family
    spans = [Region("F", 1001, 3000), Region("SH", 4001, 4600)]
    res = tile_spans(aln, cons, stats, cmap, _span_params(), spans)
    assert [a.name for a in res.amplicons] == ["F_1", "F_2", "F_3", "SH"]
    assert [a.index for a in res.amplicons] == [1, 2, 3, 4]
    assert [a.pool for a in res.amplicons] == [1, 2, 1, 2]
    assert res.spanned is True


def test_tile_spans_overlaps_inside_a_span_only(family) -> None:
    aln, cons, stats, cmap = family
    spans = [Region("F", 1001, 3000), Region("SH", 4001, 4600)]
    res = tile_spans(aln, cons, stats, cmap, _span_params(), spans)
    assert res.amplicons[0].overlap_next >= 50          # inside F
    assert res.amplicons[2].overlap_next is None        # last of F: no neighbour
    assert res.amplicons[-1].overlap_next is None


def test_tile_spans_checks_dimers_across_spans(family, monkeypatch) -> None:
    aln, cons, stats, cmap = family
    monkeypatch.setattr(tiling, "dimer_dg", lambda a, b: -20_000)
    spans = [Region("F", 1001, 3000), Region("SH", 4001, 4600)]
    res = tile_spans(aln, cons, stats, cmap, _span_params(), spans)
    # F_2 and SH share pool 2, so only a scheme-wide check can pair them
    assert any("F_2" in w and "SH" in w for w in res.dimer_warnings)


def test_tile_spans_sums_candidates_and_rejects_no_spans(family) -> None:
    aln, cons, stats, cmap = family
    f_only = tile_spans(aln, cons, stats, cmap, _span_params(), [Region("F", 1001, 3000)])
    both = tile_spans(aln, cons, stats, cmap, _span_params(),
                      [Region("F", 1001, 3000), Region("SH", 4001, 4600)])
    assert both.n_fwd > f_only.n_fwd and both.n_rev > f_only.n_rev
    with pytest.raises(PipelineError, match="No regions"):
        tile_spans(aln, cons, stats, cmap, _span_params(), [])


def test_tile_spans_names_the_failing_span(family) -> None:
    aln, cons, stats, cmap = family
    # default end/start loss of 50 cannot be met at the SH start in this fixture
    with pytest.raises(PipelineError, match="Region 'SH'"):
        tile_spans(aln, cons, stats, cmap, _params(amplicon_len=700), [Region("SH", 4001, 4600)])


# --- span mode talks about --amplicon-len, the flag it derives N from -------------------

def test_min_amplicon_len_is_the_shortest_that_still_tiles() -> None:
    # One shorter always needs more amplicons than the geometry allows.
    for span_len, overlap in ((1807, 50), (1807, 120), (2000, 50), (13232, 50)):
        least = tiling.min_amplicon_len(span_len, overlap, 25)
        n = tiling.amplicons_for(span_len, least, overlap)
        assert tiling.geometry(span_len, n, overlap, tiling.TOLERANCES[0], 25) is not None
        worse = tiling.amplicons_for(span_len, least - 1, overlap)
        assert tiling.geometry(span_len, worse, overlap, tiling.TOLERANCES[0], 25) is None
    assert tiling.min_amplicon_len(1807, 50, 25) == 270
    assert tiling.min_amplicon_len(1807, 120, 25) == 458
    assert tiling.min_amplicon_len(40, 0, 25) == 0      # too short for even one amplicon


def test_span_geometry_error_names_the_amplicon_length(family) -> None:
    aln, cons, stats, cmap = family
    with pytest.raises(PipelineError) as exc:
        tile_spans(aln, cons, stats, cmap, _params(amplicon_len=200), [Region("F", 1001, 3000)])
    msg = str(exc.value)
    assert "--amplicon-len 267 or more" in msg
    assert "-n" not in msg and "use at most" not in msg


def test_whole_genome_geometry_error_still_names_n(family) -> None:
    with pytest.raises(PipelineError, match="N too large: amplicons would be"):
        tile(*family, _params(n_amplicons=40))


def test_span_no_fit_suggests_an_amplicon_length(family) -> None:
    rows = [(i, s[:1400] + "A" * 2000 + s[3400:]) for i, s in make_genome_family()]
    aln, cons, stats, cmap = _prepare(rows)
    with pytest.raises(PipelineError) as exc:
        tile_spans(aln, cons, stats, cmap, _params(amplicon_len=700), [Region("F", 1001, 4000)])
    msg = str(exc.value)
    assert "no acceptable primers between" in msg
    assert "try --amplicon-len 3000" in msg and "try -n" not in msg


def test_tile_spans_refuses_two_amplicons_with_the_same_name(family) -> None:
    aln, cons, stats, cmap = family
    # F takes three amplicons, F_1 to F_3. A sibling feature with no /gene whose product
    # sanitises to F_1 takes one, also named F_1: primers.tsv and scheme.bed would then carry
    # two different oligos called F_1_LEFT.
    spans = [Region("F", 1001, 3000), Region("F_1", 4001, 4600)]
    with pytest.raises(PipelineError, match="F_1"):
        tile_spans(aln, cons, stats, cmap, _span_params(), spans)


def test_tile_spans_keeps_distinct_names(family) -> None:
    aln, cons, stats, cmap = family
    res = tile_spans(aln, cons, stats, cmap, _span_params(),
                     [Region("F", 1001, 3000), Region("SH", 4001, 4600)])
    names = [a.name for a in res.amplicons]
    assert len(set(names)) == len(names)
    assert len({n for a in names for n in primer_names(a)}) == 2 * len(names)


# --- pools of overlapping spans ----------------------------------------------------

def _same_pool_overlaps(amplicons) -> list[tuple[str, str]]:
    return [(a.name, b.name) for i, a in enumerate(amplicons) for b in amplicons[i + 1:]
            if a.pool == b.pool and a.fwd.start <= b.rev.end and b.fwd.start <= a.rev.end]


def test_assign_pools_flips_a_span_to_avoid_a_clash() -> None:
    # span A alternates pools 1, 2; span B's first amplicon overlaps A's first, so with the
    # numbering's own pools (B_1 = pool 1) they would share a tube: B starts in pool 2 instead
    pools, warnings = tiling.assign_pools([
        [("A_1", 0, 100), ("A_2", 200, 300)],
        [("B_1", 50, 150), ("B_2", 400, 500)],
    ])
    assert pools == [1, 2, 2, 1] and warnings == []


def test_assign_pools_keeps_the_numbering_when_nothing_overlaps() -> None:
    pools, warnings = tiling.assign_pools([[("A_1", 0, 100), ("A_2", 150, 250)], [("B", 400, 500)]])
    assert pools == [1, 2, 1] and warnings == []


def test_assign_pools_warns_about_an_unavoidable_clash() -> None:
    # three amplicons that overlap one another cannot be split over two pools
    pools, warnings = tiling.assign_pools([[("A", 0, 100)], [("B", 50, 150)], [("C", 60, 160)]])
    assert len(warnings) == 1 and "C" in warnings[0] and "same pool" in warnings[0]


def test_overlapping_regions_never_share_a_pool_silently(family) -> None:
    aln, cons, stats, cmap = family
    spans = [Region("F", 1001, 3000), Region("G", 2001, 3000)]
    res = tile_spans(aln, cons, stats, cmap, _span_params(), spans)
    clashes = _same_pool_overlaps(res.amplicons)
    assert clashes, "this scheme cannot be two-coloured, so the test needs a clash to report"
    assert len(res.pool_warnings) == len(clashes)
    assert all(a in w and b in w for (a, b), w in zip(clashes, res.pool_warnings))


def test_dimers_are_checked_in_the_pools_the_amplicons_are_in(family, monkeypatch) -> None:
    aln, cons, stats, cmap = family
    monkeypatch.setattr(tiling, "dimer_dg", lambda a, b: -20_000)
    res = tile_spans(aln, cons, stats, cmap, _span_params(),
                     [Region("F", 1001, 3000), Region("G", 2001, 3000)])
    pool_of = {a.name: a.pool for a in res.amplicons}
    for warning in res.dimer_warnings:
        pool = int(warning.split(":")[0].removeprefix("pool "))
        names = [n for n in pool_of if n + "_LEFT" in warning or n + "_RIGHT" in warning]
        assert names and all(pool_of[n] == pool for n in names), warning
