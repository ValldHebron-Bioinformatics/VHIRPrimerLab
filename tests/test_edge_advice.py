"""Advice on loosening --max-start-loss / --max-end-loss when an edge amplicon is poor."""

from __future__ import annotations

import dataclasses

import pytest

from tests.synthetic import make_genome_family
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


def test_advice_when_loosening_the_start_helps(start_poor) -> None:
    res = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert res.amplicons[0].below_threshold
    [adv] = res.edge_advice
    assert (adv.span, adv.end) == (None, "start")
    assert adv.worst_before == pytest.approx(_worst(res)) and adv.worst_before < 95
    assert adv.worst_after > adv.worst_before
    assert adv.extra_loss == adv.limit - 50 > 0


def test_advice_is_reproduced_by_rerunning_with_its_flag(start_poor) -> None:
    base = _params(n_amplicons=5, max_end_loss=60)
    adv = tile(*start_poor, base).edge_advice[0]
    again = tile(*start_poor, dataclasses.replace(base, max_start_loss=adv.limit))
    assert _worst(again) == pytest.approx(adv.worst_after)
    assert again.edge_advice == []  # nothing further to offer once it is taken
    assert _worst(again) > adv.worst_before


def test_advice_for_the_end(end_poor) -> None:
    res = tile(*end_poor, _params(n_amplicons=5, max_start_loss=60))
    [adv] = res.edge_advice
    assert adv.end == "end" and res.amplicons[-1].below_threshold
    again = tile(*end_poor, _params(n_amplicons=5, max_start_loss=60, max_end_loss=adv.limit))
    assert _worst(again) == pytest.approx(adv.worst_after)


def test_the_result_is_what_was_asked_for(start_poor) -> None:
    asked = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    assert asked.edge_advice
    # same primers as a run in which the advice cannot be explored (limit at the cap)
    spy_calls = []
    original = tiling._choose_any
    def spy(*args):
        spy_calls.append(args)
        return original(*args)
    tiling._choose_any = spy
    try:
        again = tile(*start_poor, _params(n_amplicons=5, max_end_loss=60))
    finally:
        tiling._choose_any = original
    assert [(a.fwd.start, a.rev.start) for a in again.amplicons] == \
           [(a.fwd.start, a.rev.start) for a in asked.amplicons]
    assert again.amplicons[0].fwd.start <= 50  # still inside the asked-for window


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
    assert res.edge_advice[0].worst_after > res.edge_advice[0].worst_before
