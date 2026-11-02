"""Per-column conservation scoring and consensus construction."""

from __future__ import annotations

from collections import Counter

from vhirprimerlab.models import Alignment, ColumnStats, Consensus

BASE_ORDER = "ACGT"


def score_columns(aln: Alignment) -> list[ColumnStats]:
    """Majority base, its frequency, gap fraction and ambiguity count per column."""
    n = len(aln.rows)
    stats: list[ColumnStats] = []
    for col in zip(*aln.rows):
        counts = Counter(col)
        gaps = counts.get("-", 0)
        ambiguous = n - gaps - sum(counts.get(b, 0) for b in BASE_ORDER)
        best = max(BASE_ORDER, key=lambda b: (counts.get(b, 0), -BASE_ORDER.index(b)))
        best_count = counts.get(best, 0)
        stats.append(
            ColumnStats(
                majority_base=best if best_count else "N",
                majority_freq=best_count / n,
                gap_frac=gaps / n,
                n_ambiguous=ambiguous,
            )
        )
    return stats


def build_consensus(stats: list[ColumnStats]) -> Consensus:
    """Majority-base consensus of every column whose gap fraction is below 0.5."""
    cols = [i for i, s in enumerate(stats) if s.gap_frac < 0.5]
    return Consensus(seq="".join(stats[i].majority_base for i in cols), cons_to_col=cols)
