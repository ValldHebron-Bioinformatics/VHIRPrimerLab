"""score_pair and the mismatch burden must equal a plain row-by-row evaluation.

Both are the hot path of degenerate primer design: they score a primer pair against every sequence,
thousands of times. They may reuse the result for rows that are identical under the primer, because
match_primer reads nothing else; this test holds them to the plain loop on random alignments with
gaps, ambiguity codes and degenerate primers, so a shortcut that changed an answer would be caught.
"""

from __future__ import annotations

import random

import pytest

from vhirprimerlab.degenerate import _mismatch_burden
from vhirprimerlab.models import Alignment, PairCoverage, ScoringParams
from vhirprimerlab.scoring import match_primer, score_pair

WIDTH = 90


def _alignment(rng: random.Random, n_rows: int) -> Alignment:
    base = [rng.choice("ACGT") for _ in range(WIDTH)]
    rows = []
    for _ in range(n_rows):
        row = list(base)
        for i in range(WIDTH):
            r = rng.random()
            if r < 0.06:
                row[i] = rng.choice("ACGT")
            elif r < 0.08:
                row[i] = "-"
            elif r < 0.09:
                row[i] = rng.choice("NRY")
        rows.append("".join(row))
    # many rows identical to each other under the primers, as in a real set of related genomes
    rows += rows[: n_rows // 2]
    return Alignment([f"s{i}" for i in range(len(rows))], rows)


def _site(rng: random.Random, aln: Alignment, cols: list[int]) -> str:
    return "".join(rng.choice([aln.rows[0][c] if aln.rows[0][c] in "ACGT" else "A", "R", "Y", "N"])
                   for c in cols)


def _plain_score_pair(fsite, fcols, rsite, rcols, aln, params) -> PairCoverage:
    covered, not_covered = [], {}
    for sid, row in zip(aln.ids, aln.rows):
        f = match_primer(fsite, fcols, row, "right", params)
        r = match_primer(rsite, rcols, row, "left", params)
        if f.ok and r.ok:
            covered.append(sid)
        else:
            reasons = [f"fwd {f.reason}"] if f.reason else []
            reasons += [f"rev {r.reason}"] if r.reason else []
            not_covered[sid] = "; ".join(reasons)
    return PairCoverage(covered, not_covered, 100.0 * len(covered) / len(aln.ids) if aln.ids else 0.0)


@pytest.mark.parametrize("seed", range(12))
def test_score_pair_and_burden_equal_the_plain_loop(seed: int) -> None:
    rng = random.Random(seed)
    params = ScoringParams(max_mismatches=rng.choice([1, 2, 3]), window=rng.choice([3, 5]),
                           max_ambiguous=rng.choice([0, 2]))
    aln = _alignment(rng, rng.choice([6, 20, 41]))
    fcols = list(range(8, 20)) + list(range(21, 29))        # column 20 is an insertion for the primer
    rcols = list(range(55, 62)) + list(range(63, 76))
    fsite, rsite = _site(rng, aln, fcols), _site(rng, aln, rcols)
    got = score_pair(fsite, fcols, rsite, rcols, aln, params)
    want = _plain_score_pair(fsite, fcols, rsite, rcols, aln, params)
    assert (got.covered, got.not_covered, got.coverage_pct) == (want.covered, want.not_covered, want.coverage_pct)
    burden = _mismatch_burden({"fwd": fsite, "rev": rsite}, {"fwd": fcols, "rev": rcols}, aln.rows, params)
    assert burden == sum(match_primer(fsite, fcols, row, "right", params).mismatches
                         + match_primer(rsite, rcols, row, "left", params).mismatches for row in aln.rows)
