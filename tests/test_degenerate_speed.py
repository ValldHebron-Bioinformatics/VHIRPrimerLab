"""degenerate_site keeps its answers however it is computed: compare it with the plain version.

``reference_degenerate_site`` is the straightforward greedy rule (every trial re-scored from
scratch with match_primer). The shipped version updates each sequence's mismatch count from the
one base a trial changes, so it must give the same site on any alignment, gaps and ambiguity
codes included.
"""

from __future__ import annotations

import random
from collections import Counter

import pytest

from vhirprimerlab.degenerate import BASES, _gc_ok, degenerate_site, fold
from vhirprimerlab.models import Alignment, ScoringParams
from vhirprimerlab.scoring import CODE_FOR, IUPAC, match_primer


def reference_degenerate_site(site, cols, aln, params, max_fold, three_prime, gc_limits=None, accept=None):
    n = len(site)
    terminal = n - 1 if three_prime == "right" else 0
    lo, hi = cols[0], cols[-1] + 1
    weight = Counter(row[lo:hi] for row in aln.rows)
    offset = [c - lo for c in cols]

    def misses(s, slices):
        return [x for x in slices if not match_primer(s, offset, x, three_prime, params).ok]

    def burden(s, slices):
        return sum(weight[x] * match_primer(s, offset, x, three_prime, params).mismatches for x in slices)

    unmatched = misses(site, list(weight))
    current = site
    best_site, best_missing = site, sum(weight[x] for x in unmatched)
    while unmatched:
        current_burden = burden(current, unmatched)
        n_missing = sum(weight[x] for x in unmatched)
        trials = []
        for i in range(n):
            if i == terminal:
                continue
            allowed = IUPAC[current[i]]
            seen = {x[offset[i]] for x in unmatched}
            for b in sorted((seen & set(BASES)) - allowed, key=BASES.index):
                trial = current[:i] + CODE_FOR[allowed | {b}] + current[i + 1:]
                if fold(trial) > max_fold or not _gc_ok(trial, gc_limits):
                    continue
                still = misses(trial, unmatched)
                gain = n_missing - sum(weight[x] for x in still)
                removed = current_burden - burden(trial, unmatched)
                dist = n - 1 - i if three_prime == "right" else i
                if gain > 0 or removed > 0:
                    trials.append(((gain, removed, dist, -BASES.index(b)), trial, still))
        trials.sort(key=lambda t: t[0], reverse=True)
        best = next((t for t in trials if accept is None or accept(t[1])), None)
        if best is None:
            break
        _, current, unmatched = best
        missing = sum(weight[x] for x in unmatched)
        if missing < best_missing:
            best_site, best_missing = current, missing
    return best_site


def random_alignment(rng: random.Random, n_rows: int, width: int, noise: float) -> Alignment:
    """Rows that differ from one base sequence by substitutions, gaps, N and other codes."""
    base = [rng.choice(BASES) for _ in range(width)]
    rows = []
    for _ in range(n_rows):
        row = base[:]
        for p in range(width):
            r = rng.random()
            if r < noise:
                row[p] = rng.choice(BASES)
            elif r < noise * 1.25:
                row[p] = rng.choice("RYKMSWN")
            elif r < noise * 1.4:
                row[p] = "-"
        rows.append("".join(row))
    return Alignment([f"s{i}" for i in range(n_rows)], rows)


@pytest.mark.parametrize("seed", range(120))
def test_same_site_as_the_plain_version(seed: int) -> None:
    rng = random.Random(seed)
    width = 70
    aln = random_alignment(rng, rng.randint(3, 40), width, rng.choice([0.02, 0.06, 0.15]))
    start = rng.randint(5, 30)
    length = rng.randint(18, 25)
    cols = list(range(start, start + length))

    def consensus_base(col: int) -> str:
        top = Counter(r[col] for r in aln.rows).most_common(1)[0][0]
        return top if top in BASES else "A"

    site = "".join(consensus_base(c) for c in cols)  # the primer is written in consensus bases
    params = ScoringParams(rng.randint(0, 3), rng.randint(2, 5), rng.choice([0, 2]))
    side = rng.choice(["right", "left"])
    max_fold = rng.choice([2, 4, 8])
    gc = rng.choice([None, (30.0, 70.0), (0.0, 100.0)])
    accept = rng.choice([None, lambda s: s.count("R") + s.count("Y") <= 1])
    want = reference_degenerate_site(site, cols, aln, params, max_fold, side, gc, accept)
    assert degenerate_site(site, cols, aln, params, max_fold, side, gc, accept) == want
    # the caller may hand over the starting scores of every distinct slice
    lo, hi = cols[0], cols[-1] + 1
    offset = [c - lo for c in cols]
    initial = {}
    for row in aln.rows:
        x = row[lo:hi]
        if x not in initial:
            r = match_primer(site, offset, x, side, params)
            initial[x] = (r.mismatches, r.three_prime_mismatches)
    assert degenerate_site(site, cols, aln, params, max_fold, side, gc, accept, initial=initial) == want
