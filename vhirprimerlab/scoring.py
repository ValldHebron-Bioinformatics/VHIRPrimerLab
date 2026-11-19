"""Re-scoring of primers against every aligned sequence.

A primer "site" is the oligo as it reads on the forward strand: the forward primer as-is,
the reverse primer reverse-complemented. The forward primer's 3' end is the right end of
its site; the reverse primer's 3' end is the left end.
"""

from __future__ import annotations

from typing import Literal

from vhirprimerlab.models import Alignment, MatchResult, PairCoverage, ScoringParams

IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("CG"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}
CODE_FOR: dict[frozenset[str], str] = {v: k for k, v in IUPAC.items()}

ThreePrime = Literal["right", "left"]


def match_primer(
    site: str,
    cols: list[int],
    row: str,
    three_prime: ThreePrime,
    params: ScoringParams,
) -> MatchResult:
    """Count mismatches of a primer site against one alignment row."""
    n = len(site)
    w = min(params.window, n)
    if three_prime == "right":
        window_idx = range(n - w, n)
        win_lo, win_hi = cols[n - w], cols[-1]
    else:
        window_idx = range(0, w)
        win_lo, win_hi = cols[0], cols[w - 1]

    mismatches = 0
    in_window = 0
    unknown_in_window = 0
    forgiven = 0
    for i, col in enumerate(cols):
        base = row[col]
        wanted = IUPAC[site[i]]
        have = IUPAC.get(base) if base != "-" else None
        if have is not None and have <= wanted:
            continue  # the base, or a code that only names bases the primer allows here
        if have is not None and have & wanted:
            # N, R, ...: the sequence may or may not match here; it is unknown, not different
            if i in window_idx:
                unknown_in_window += 1  # at the 3' end an unknown base cannot be trusted
            elif forgiven < params.max_ambiguous:
                forgiven += 1
            else:
                mismatches += 1
            continue
        mismatches += 1
        if i in window_idx:
            in_window += 1

    primer_cols = set(cols)
    for col in range(cols[0], cols[-1] + 1):
        if col not in primer_cols and row[col] != "-":
            mismatches += 1  # base inserted in this sequence relative to the primer
            if win_lo < col < win_hi:
                in_window += 1

    if in_window:
        reason: str | None = "3' mismatch"
    elif unknown_in_window:
        reason = "ambiguous 3' base"
    elif mismatches > params.max_mismatches:
        reason = "too many mismatches"
    else:
        reason = None
    return MatchResult(reason is None, mismatches, in_window + unknown_in_window, reason, forgiven)


class RowMatcher:
    """match_primer, worked out once per distinct primer site and stretch of sequence under it.

    match_primer reads nothing of a row but the columns from the primer's first to its last, so two
    rows that agree there get the same answer. Related genomes agree there far more often than not,
    and degenerate design asks for the same site against every row thousands of times, so a matcher
    kept for the length of one design saves most of those calls without changing a result.
    """

    def __init__(self, params: ScoringParams) -> None:
        self.params = params
        self._memo: dict[tuple, MatchResult] = {}
        self._cols: dict[int, tuple[list[int], tuple[int, ...]]] = {}

    def __call__(self, site: str, cols: list[int], row: str, three_prime: ThreePrime) -> MatchResult:
        entry = self._cols.get(id(cols))
        if entry is None or entry[0] is not cols:   # keeps the list alive, so its id stays its own
            entry = self._cols[id(cols)] = (cols, tuple(cols))
        key = (site, three_prime, entry[1], row[cols[0]:cols[-1] + 1])
        found = self._memo.get(key)
        if found is None:
            found = self._memo[key] = match_primer(site, cols, row, three_prime, self.params)
        return found


def score_pair(
    fwd_site: str,
    fwd_cols: list[int],
    rev_site: str,
    rev_cols: list[int],
    aln: Alignment,
    params: ScoringParams,
    matcher: RowMatcher | None = None,
) -> PairCoverage:
    """Coverage of a primer pair: sequences where both primers match.

    A matcher kept by the caller across calls (see RowMatcher) makes repeated scoring of the same
    primers against the same rows cheap; without one, rows that agree under a primer share a result
    within this call."""
    covered: list[str] = []
    not_covered: dict[str, str] = {}
    match = matcher or RowMatcher(params)
    for sid, row in zip(aln.ids, aln.rows):
        fwd = match(fwd_site, fwd_cols, row, "right")
        rev = match(rev_site, rev_cols, row, "left")
        if fwd.ok and rev.ok:
            covered.append(sid)
        else:
            reasons = [f"fwd {fwd.reason}"] if fwd.reason else []
            reasons += [f"rev {rev.reason}"] if rev.reason else []
            not_covered[sid] = "; ".join(reasons)
    pct = 100.0 * len(covered) / len(aln.ids) if aln.ids else 0.0
    return PairCoverage(covered, not_covered, pct)
