"""Greedy introduction of IUPAC degenerate bases to rescue uncovered sequences."""

from __future__ import annotations

import itertools
import math
from collections import Counter

import primer3

from vhirprimerlab.models import Alignment, PairCoverage, ScoringParams
from vhirprimerlab.scoring import CODE_FOR, IUPAC, ThreePrime, match_primer, score_pair

BASES = "ACGT"


def fold(site: str) -> int:
    """Number of distinct oligos a degenerate primer stands for."""
    return math.prod(len(IUPAC[b]) for b in site)


def degenerate_bases(primer: str) -> int:
    """Number of IUPAC degenerate (non-ACGT) positions in a primer."""
    return sum(b not in BASES for b in primer)


def gc_percent(primer: str) -> float:
    """GC % averaged over every variant of a (possibly degenerate) primer."""
    return 100 * sum(len(IUPAC[b] & {"G", "C"}) / len(IUPAC[b]) for b in primer) / len(primer)


def tm_range(primer: str) -> tuple[float, float]:
    """Min and max primer3 Tm over every expansion of a (possibly degenerate) primer."""
    options = [sorted(IUPAC[b]) for b in primer]
    tms = [primer3.calc_tm("".join(v)) for v in itertools.product(*options)]
    return min(tms), max(tms)


def _gc_ok(site: str, gc_limits: tuple[float, float] | None) -> bool:
    return gc_limits is None or gc_limits[0] <= gc_percent(site) <= gc_limits[1]


def _mismatch_burden(
    sites: dict[str, str], cols: dict[str, list[int]], rows: list[str], params: ScoringParams
) -> int:
    """Total primer mismatches over all sequences (lower is better)."""
    return sum(
        match_primer(sites["fwd"], cols["fwd"], row, "right", params).mismatches
        + match_primer(sites["rev"], cols["rev"], row, "left", params).mismatches
        for row in rows
    )


def degenerate_pair(
    fwd_site: str,
    fwd_cols: list[int],
    rev_site: str,
    rev_cols: list[int],
    aln: Alignment,
    params: ScoringParams,
    max_fold: int,
    gc_limits: tuple[float, float] | None = None,
) -> tuple[str, str, PairCoverage]:
    """Add one IUPAC base at a time, choosing the change that rescues most sequences.

    Changes are ranked by sequences rescued, then by mismatches removed (so a sequence
    needing several changes can still be rescued), then by distance from the 3' end.
    The result is cut back to the last change that rescued a sequence, so no
    degeneracy is added without benefit. Sites are in forward-strand orientation; the
    3'-terminal base of each primer (last base of the forward site, first base of the
    reverse site) is never changed. With gc_limits, changes that take a primer's GC %
    (averaged over its variants) outside the limits are skipped.
    """
    sites = {"fwd": fwd_site, "rev": rev_site}
    cols = {"fwd": fwd_cols, "rev": rev_cols}
    rows = dict(zip(aln.ids, aln.rows))
    cov = score_pair(sites["fwd"], fwd_cols, sites["rev"], rev_cols, aln, params)
    burden = _mismatch_burden(sites, cols, aln.rows, params)
    best_state = (dict(sites), cov)

    while cov.not_covered:
        best: tuple[tuple[int, int, int, bool, int], dict[str, str], PairCoverage, int] | None = None
        for which in ("fwd", "rev"):
            site = sites[which]
            n = len(site)
            terminal = n - 1 if which == "fwd" else 0
            for i in range(n):
                if i == terminal:
                    continue
                current = IUPAC[site[i]]
                seen = {rows[sid][cols[which][i]] for sid in cov.not_covered}
                for b in sorted((seen & set(BASES)) - current, key=BASES.index):
                    trial_site = site[:i] + CODE_FOR[current | {b}] + site[i + 1:]
                    if fold(trial_site) > max_fold or not _gc_ok(trial_site, gc_limits):
                        continue
                    trial = dict(sites, **{which: trial_site})
                    trial_cov = score_pair(trial["fwd"], fwd_cols, trial["rev"], rev_cols, aln, params)
                    trial_burden = _mismatch_burden(trial, cols, aln.rows, params)
                    gain = len(trial_cov.covered) - len(cov.covered)
                    dist_from_3p = n - 1 - i if which == "fwd" else i
                    key = (gain, burden - trial_burden, dist_from_3p, which == "fwd", -BASES.index(b))
                    if best is None or key > best[0]:
                        best = (key, trial, trial_cov, trial_burden)
        if best is None or (best[0][0] <= 0 and best[0][1] <= 0):
            break
        _, sites, cov, burden = best
        if len(cov.covered) > len(best_state[1].covered):
            best_state = (dict(sites), cov)

    final_sites, final_cov = best_state
    return final_sites["fwd"], final_sites["rev"], final_cov


def degenerate_site(
    site: str,
    cols: list[int],
    aln: Alignment,
    params: ScoringParams,
    max_fold: int,
    three_prime: ThreePrime,
    gc_limits: tuple[float, float] | None = None,
) -> str:
    """degenerate_pair's greedy rule for a single primer site (forward-strand orientation).

    Adding a degenerate base never breaks a match, so only rows not yet matched are
    re-scored. The 3'-terminal base is never changed, and with gc_limits no change may
    take the GC % outside them; the state matching the most
    rows (earliest on ties) is returned.
    """
    n = len(site)
    terminal = n - 1 if three_prime == "right" else 0
    # identical row slices give identical results: score each distinct slice once
    lo, hi = cols[0], cols[-1] + 1
    weight = Counter(row[lo:hi] for row in aln.rows)
    offset = [c - lo for c in cols]

    def misses(s: str, slices: list[str]) -> list[str]:
        return [x for x in slices if not match_primer(s, offset, x, three_prime, params).ok]

    def burden(s: str, slices: list[str]) -> int:
        return sum(weight[x] * match_primer(s, offset, x, three_prime, params).mismatches for x in slices)

    unmatched = misses(site, list(weight))
    current = site
    best_site, best_missing = site, sum(weight[x] for x in unmatched)
    while unmatched:
        current_burden = burden(current, unmatched)
        n_missing = sum(weight[x] for x in unmatched)
        best: tuple[tuple[int, int, int, int], str, list[str]] | None = None
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
                key = (gain, removed, dist, -BASES.index(b))
                if best is None or key > best[0]:
                    best = (key, trial, still)
        if best is None or (best[0][0] <= 0 and best[0][1] <= 0):
            break
        _, current, unmatched = best
        missing = sum(weight[x] for x in unmatched)
        if missing < best_missing:
            best_site, best_missing = current, missing
    return best_site
