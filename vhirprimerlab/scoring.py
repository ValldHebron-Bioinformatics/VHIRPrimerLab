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
    for i, col in enumerate(cols):
        s = row[col]
        if s not in IUPAC[site[i]] or s not in "ACGT":
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
    elif mismatches > params.max_mismatches:
        reason = "too many mismatches"
    else:
        reason = None
    return MatchResult(reason is None, mismatches, in_window, reason)


def score_pair(
    fwd_site: str,
    fwd_cols: list[int],
    rev_site: str,
    rev_cols: list[int],
    aln: Alignment,
    params: ScoringParams,
) -> PairCoverage:
    """Coverage of a primer pair: sequences where both primers match."""
    covered: list[str] = []
    not_covered: dict[str, str] = {}
    for sid, row in zip(aln.ids, aln.rows):
        fwd = match_primer(fwd_site, fwd_cols, row, "right", params)
        rev = match_primer(rev_site, rev_cols, row, "left", params)
        if fwd.ok and rev.ok:
            covered.append(sid)
        else:
            reasons = [f"fwd {fwd.reason}"] if fwd.reason else []
            reasons += [f"rev {rev.reason}"] if rev.reason else []
            not_covered[sid] = "; ".join(reasons)
    pct = 100.0 * len(covered) / len(aln.ids) if aln.ids else 0.0
    return PairCoverage(covered, not_covered, pct)
