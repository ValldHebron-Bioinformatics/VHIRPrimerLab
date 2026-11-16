"""Advice on loosening --max-start-loss / --max-end-loss when an edge amplicon is poor."""

from __future__ import annotations

import dataclasses
import random

import pytest

from tests.synthetic import _substitute, make_genome_family
from tests.test_tiling import PARAMS, _params, _prepare
from vhirprimerlab import tiling
from vhirprimerlab.models import Region
from vhirprimerlab.tiling import tile, tile_spans

# 20 sequences whose first 70 bases are hypervariable: a primer inside the default 50 base
# window matches few of them, one a little further in matches all.
START_POOR = dict(hyper=(0, 70), seed=7, n_seqs=20)
END_POOR = dict(hyper=(5930, 6000), seed=7, n_seqs=20)


@pytest.fixture(scope="module")
def start_poor():
    return _prepare(make_genome_family(**START_POOR))


@pytest.fixture(scope="module")
def end_poor():
    return _prepare(make_genome_family(**END_POOR))


def _worst(res) -> float:
    return min(a.coverage.coverage_pct for a in res.amplicons)


def _edge(res, end: str) -> float:
    return res.amplicons[0 if end == "start" else -1].coverage.coverage_pct


def _interior_poor():
    """START_POOR with a wide variable block inside: interior amplicons end up weaker than the edge."""
    rng = random.Random(1)
    rows = make_genome_family(**START_POOR)
    out = [rows[0]]
    for rid, seq in rows[1:]:
        seq = list(seq)
        for i in range(1500, 4500):
            if rng.random() < 0.4:
                seq[i] = _substitute(rng, seq[i])
        out.append((rid, "".join(seq)))
    return _prepare(out)


@pytest.fixture(scope="module")
def interior_poor():
    return _interior_poor()


def _coords(res):
    return [(a.fwd.start, a.fwd.seq, a.rev.start, a.rev.seq) for a in res.amplicons]


def test_advice_when_loosening_the_start_helps(start_poor) -> None:
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert res.amplicons[0].below_threshold
    [adv] = res.edge_advice
    assert (adv.span, adv.end) == (None, "start")
    assert adv.edge_before == pytest.approx(_edge(res, "start")) and adv.edge_before < 95
    assert adv.edge_after > adv.edge_before
    assert adv.extra_loss == adv.limit - 50 > 0


def test_advice_is_reproduced_by_rerunning_with_its_flag(start_poor) -> None:
    base = _params(n_amplicons=5, max_end_loss=60)
    adv = tile(*start_poor, base).edge_advice[0]
    again = tile(*start_poor, dataclasses.replace(base, max_start_loss=adv.limit))
    assert _edge(again, "start") == pytest.approx(adv.edge_after)
    assert again.edge_advice == []  # nothing further to offer once it is taken
    assert _edge(again, "start") > adv.edge_before


def test_advice_for_the_end(end_poor) -> None:
    res = tile(*end_poor, _params(n_amplicons=5, max_start_loss=60))
    [adv] = res.edge_advice
    assert adv.end == "end" and res.amplicons[-1].below_threshold
    again = tile(*end_poor, _params(n_amplicons=5, max_start_loss=60, max_end_loss=adv.limit))
    assert _edge(again, "end") == pytest.approx(adv.edge_after)


def test_the_result_is_what_was_asked_for(start_poor, monkeypatch) -> None:
    asked = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert asked.edge_advice
    # the same input with the advice switched off: the primers must be identical
    monkeypatch.setattr(tiling, "_edge_advice", lambda *a, **k: [])
    plain = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert plain.edge_advice == []
    assert _coords(asked) == _coords(plain)
    assert asked.amplicons[0].fwd.end + 1 <= 50  # still inside the asked-for window


def test_the_first_answer_is_reached_by_its_own_limit(start_poor, monkeypatch) -> None:
    # The limit named by the screen on the primers as first scored must be one that really
    # delivers: allowing it must not swap the primer it was named for. It used to (the end rules
    # replaced the ordinary primer once the limit reached it), so 104 reached only 40 % and a
    # second pass had to find 108. Now one screen names 104 and 104 reaches 100 %.
    limits = []
    original = tiling._candidates
    def spy(aln, cons, stats, params, span, progress=None, memo=None):
        limits.append(params.max_start_loss)
        return original(aln, cons, stats, params, span, progress, memo)
    monkeypatch.setattr(tiling, "_candidates", spy)
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert len(limits) == 2 and limits[0] == 50  # the run, then the screen under the named limit
    [adv] = res.edge_advice
    assert adv.limit == limits[1] and adv.edge_after == 100.0


def test_advice_is_judged_on_the_edge_not_the_weakest(interior_poor) -> None:
    res = tile(*interior_poor, _params(n_amplicons=5, max_end_loss=60))
    edge = _edge(res, "start")
    assert _worst(res) < edge < 95  # an interior amplicon is the weakest
    [adv] = res.edge_advice
    assert adv.edge_before == pytest.approx(edge) and adv.edge_after > edge
    again = tile(*interior_poor, _params(n_amplicons=5, max_end_loss=60, max_start_loss=adv.limit))
    assert _edge(again, "start") == pytest.approx(adv.edge_after)
    assert _worst(again) >= _worst(res)  # the weakest is no worse


def test_advice_promises_what_dimer_avoidance_leaves(interior_poor, monkeypatch) -> None:
    base = _params(n_amplicons=5, max_end_loss=60)
    first = tile(*interior_poor, base).edge_advice[0]
    clean = tile(*interior_poor, dataclasses.replace(base, max_start_loss=first.limit))
    bad = clean.amplicons[0].fwd.seq
    # that edge primer now forms a dimer with everything, so a re-run swaps it out (the
    # interior amplicon is the weakest, so the swap costs no coverage and is kept)
    monkeypatch.setattr(tiling, "dimer_dg", lambda a, b, t=None: -12000.0 if bad in (a, b) else -3000.0)
    [adv] = tile(*interior_poor, base).edge_advice
    assert adv.edge_after < first.edge_after  # no longer the 100 % of the dimer-prone primer
    rerun = tile(*interior_poor, dataclasses.replace(base, max_start_loss=adv.limit))
    assert rerun.amplicons[0].fwd.seq != bad
    assert _edge(rerun, "start") == pytest.approx(adv.edge_after)


def test_only_the_poor_end_is_loosened_and_to_one_amplicon(start_poor, monkeypatch) -> None:
    seen = []
    original = tiling._choose_any
    def spy(fwds, revs, n, geoms, params, span, banned):
        seen.append(params)
        return original(fwds, revs, n, geoms, params, span, banned)
    monkeypatch.setattr(tiling, "_choose_any", spy)
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    asked, explored, *checks = seen
    assert asked.max_start_loss == 50 and explored.max_start_loss == res.geometry.target_len
    assert all(p.max_end_loss == 60 for p in seen)  # the end that is fine is never loosened
    assert res.edge_advice[0].limit in [p.max_start_loss for p in checks]  # the answer was checked


def test_a_clean_run_costs_no_extra_search(family_clean, monkeypatch) -> None:
    calls = []
    original = tiling._choose_any
    monkeypatch.setattr(tiling, "_choose_any", lambda *a: calls.append(a) or original(*a))
    res = tile(*family_clean, PARAMS)
    assert not any(a.below_threshold for a in res.amplicons)
    assert res.edge_advice == [] and len(calls) == 1


@pytest.fixture(scope="module")
def family_clean():
    return _prepare(make_genome_family())


def test_no_advice_when_loosening_does_not_improve(start_poor, monkeypatch) -> None:
    original = tiling._choose_any
    first = []
    def same(fwds, revs, n, geoms, params, span, banned):
        out = original(fwds, revs, n, geoms, params, span, banned)
        first.append(out)
        return first[0]  # the loosened search finds nothing better than the first
    monkeypatch.setattr(tiling, "_choose_any", same)
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert res.amplicons[0].below_threshold and res.edge_advice == []


def test_no_advice_when_loosening_finds_no_tiling(start_poor, monkeypatch) -> None:
    original = tiling._choose_any
    calls = []
    def none_second(fwds, revs, n, geoms, params, span, banned):
        calls.append(1)
        out = original(fwds, revs, n, geoms, params, span, banned)
        return out if len(calls) == 1 else (None, geoms[-1])
    monkeypatch.setattr(tiling, "_choose_any", none_second)
    assert tile(*start_poor, _params(n_amplicons=5, max_end_loss=60)).edge_advice == []


def test_no_advice_when_the_check_finds_it_worse(start_poor, monkeypatch) -> None:
    original = tiling._choose_any
    calls: list = []
    def worse_check(fwds, revs, n, geoms, params, span, banned):
        out = original(fwds, revs, n, geoms, params, span, banned)
        calls.append(out)
        return calls[0] if len(calls) >= 3 else out  # the check turns out no better than the run
    monkeypatch.setattr(tiling, "_choose_any", worse_check)
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert len(calls) == 3 and res.edge_advice == []  # screened, checked, nothing to advise


def test_no_advice_when_the_check_finds_no_tiling(start_poor, monkeypatch) -> None:
    original = tiling._choose_any
    calls = []
    def nothing_on_check(fwds, revs, n, geoms, params, span, banned):
        calls.append(1)
        if len(calls) >= 3:
            return None, geoms[-1]
        return original(fwds, revs, n, geoms, params, span, banned)
    monkeypatch.setattr(tiling, "_choose_any", nothing_on_check)
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert len(calls) == 3 and res.edge_advice == []


def test_no_advice_when_it_would_make_the_weakest_worse(start_poor, monkeypatch) -> None:
    calls = []
    original = tiling._choose_any
    monkeypatch.setattr(tiling, "_choose_any", lambda *a: calls.append(a) or original(*a))
    # the explored tiling lifts the edge but drops an interior amplicon below today's weakest
    monkeypatch.setattr(tiling, "_coverages", lambda path, aln, params: [100.0, 1.0] + [100.0] * (len(path) - 2))
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert res.edge_advice == [] and len(calls) == 2  # run, screen: no costly check


def test_no_advice_when_the_limit_is_already_one_amplicon(start_poor, monkeypatch) -> None:
    calls = []
    original = tiling._choose_any
    monkeypatch.setattr(tiling, "_choose_any", lambda *a: calls.append(a) or original(*a))
    res = tile(*start_poor, _params(n_amplicons=5, max_start_loss=2000, max_end_loss=60))
    assert res.edge_advice == [] and len(calls) == 1


def test_tile_spans_names_the_span(start_poor) -> None:
    res = tile_spans(*start_poor, _params(amplicon_len=1000, max_end_loss=60),
                     [Region("F", 1, 3000), Region("G", 3201, 5800)])
    assert [(a.span, a.end) for a in res.edge_advice] == [("F", "start")]
    assert res.edge_advice[0].edge_after > res.edge_advice[0].edge_before


# --- the candidate memo: same primers, found once -------------------------------------

def test_a_memo_gives_the_primers_a_fresh_run_gives(start_poor) -> None:
    from vhirprimerlab.candidates import CandidateMemo

    aln, cons, stats, _ = start_poor
    base = _params(n_amplicons=5, max_end_loss=60)
    span = (0, len(cons.seq) - 1)
    memo = CandidateMemo(base, span, aln)
    tiling._candidates(aln, cons, stats, base, span, memo=memo)           # fills the memo
    for start_loss in (50, 80, 104, 200):                                 # the allowances advice explores
        params = dataclasses.replace(base, max_start_loss=start_loss)
        assert (tiling._candidates(aln, cons, stats, params, span, memo=memo)
                == tiling._candidates(aln, cons, stats, params, span))


def test_a_memo_is_only_used_for_the_parameters_it_was_made_for(start_poor) -> None:
    from vhirprimerlab.candidates import CandidateMemo

    aln, cons, stats, _ = start_poor
    base = _params(n_amplicons=5)
    span = (0, len(cons.seq) - 1)
    memo = CandidateMemo(base, span, aln)
    assert memo.fits(dataclasses.replace(base, max_start_loss=999), span, aln)      # allowances may move
    assert not memo.fits(dataclasses.replace(base, tm_max=base.tm_max + 1), span, aln)
    assert not memo.fits(base, (5, span[1]), aln)
    other = tiling._candidates(aln, cons, stats, dataclasses.replace(base, tm_max=base.tm_max + 1), span, memo=memo)
    assert other == tiling._candidates(aln, cons, stats, dataclasses.replace(base, tm_max=base.tm_max + 1), span)
