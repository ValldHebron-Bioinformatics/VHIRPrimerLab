"""Checking a primer's structure only when it is about to be kept must keep the same primers."""

from __future__ import annotations

import dataclasses

import pytest

from tests.synthetic import make_genome_family
from tests.test_tiling import PARAMS, _params, _prepare
from vhirprimerlab import candidates, tiling
from vhirprimerlab.candidates import enumerate_primers, pick_conserved
from vhirprimerlab.structure import structure_ok


@pytest.fixture(scope="module")
def family():
    return _prepare(make_genome_family())


def eager_picks(cons, stats, params, span):
    """Filter every anchor's primer by structure first, then rank: the behaviour to preserve."""
    raw_f, raw_r = enumerate_primers(cons.seq, params, span)
    return (pick_conserved(raw_f, stats, cons, params), pick_conserved(raw_r, stats, cons, params))


@pytest.mark.parametrize("start_loss, end_loss", [(50, 50), (0, 0), (200, 120), (1, 400)])
def test_same_primers_when_the_structure_check_is_deferred(family, start_loss, end_loss) -> None:
    aln, cons, stats, _ = family
    params = dataclasses.replace(_params(n_amplicons=5), max_start_loss=start_loss, max_end_loss=end_loss)
    span = (0, len(cons.seq) - 1)
    want = eager_picks(cons, stats, params, span)
    raw_f, raw_r = enumerate_primers(cons.seq, params, span, check_structure=False)
    valid = lambda c: structure_ok(c.seq, params.thermo)  # noqa: E731
    got = (pick_conserved(raw_f, stats, cons, params, valid=valid),
           pick_conserved(raw_r, stats, cons, params, valid=valid))
    assert got == want


def test_deferring_the_check_asks_about_far_fewer_primers(family, monkeypatch) -> None:
    aln, cons, stats, _ = family
    params = _params(n_amplicons=5)
    span = (0, len(cons.seq) - 1)
    calls = {"eager": 0, "lazy": 0}
    real = candidates.structure_ok

    def counting(kind):
        def fn(seq, thermo=None):
            calls[kind] += 1
            return real(seq) if thermo is None else real(seq, thermo)
        return fn

    monkeypatch.setattr(candidates, "structure_ok", counting("eager"))
    eager_picks(cons, stats, params, span)
    monkeypatch.setattr(candidates, "structure_ok", counting("lazy"))
    raw_f, raw_r = enumerate_primers(cons.seq, params, span, check_structure=False)
    valid = lambda c: candidates.structure_ok(c.seq, params.thermo)  # noqa: E731
    pick_conserved(raw_f, stats, cons, params, valid=valid)
    pick_conserved(raw_r, stats, cons, params, valid=valid)
    assert calls["lazy"] < calls["eager"] / 3


def test_the_tiling_picks_the_primers_the_eager_check_would_have(family) -> None:
    """tile() defers the check; the primers it works with must be those the eager one picks."""
    aln, cons, stats, cmap = family
    params, span = PARAMS, (0, len(cons.seq) - 1)
    wanted_f, wanted_r = eager_picks(cons, stats, params, span)
    scored_f, scored_r = tiling._candidates(aln, cons, stats, params, span)
    assert [(c.start, c.length) for c in scored_f] == [(c.start, c.length) for c in wanted_f]
    assert [(c.start, c.length) for c in scored_r] == [(c.start, c.length) for c in wanted_r]
