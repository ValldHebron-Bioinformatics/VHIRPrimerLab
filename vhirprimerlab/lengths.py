"""Length outlier detection relative to the median sequence length."""

from __future__ import annotations

import logging
import statistics
from typing import Literal

from vhirprimerlab.models import LengthFilterResult, PipelineError, Record

log = logging.getLogger(__name__)

OutlierAction = Literal["warn", "drop", "abort"]


def ungapped_length(seq: str) -> int:
    return len(seq) - seq.count("-")


def filter_by_length(
    records: list[Record],
    tolerance: float,
    action: OutlierAction,
    reference_id: str | None,
) -> LengthFilterResult:
    """Flag sequences whose length differs from the median by more than tolerance × median.

    reference_id is the user-chosen reference (None when defaulted); dropping it is an error.
    """
    lengths = {r.id: ungapped_length(r.seq) for r in records}
    median = float(statistics.median(lengths.values()))
    limit = tolerance * median
    outliers = [(r.id, lengths[r.id]) for r in records if abs(lengths[r.id] - median) > limit]
    outlier_ids = {rid for rid, _ in outliers}

    for rid, length in outliers:
        log.warning(
            "Length outlier: %s (%d bp; median %.0f bp, tolerance ±%.0f bp)",
            rid, length, median, limit,
        )
    if outliers and action == "abort":
        listed = ", ".join(f"{rid} ({n} bp)" for rid, n in outliers)
        raise PipelineError(f"Length outliers found (--outlier-action abort): {listed}")

    if action == "drop":
        if reference_id is not None and reference_id in outlier_ids:
            raise PipelineError(
                f"The reference sequence '{reference_id}' is a length outlier and would be "
                "dropped; choose another --reference-id or adjust --length-tolerance"
            )
        kept = [r for r in records if r.id not in outlier_ids]
        if outliers:
            log.info("Dropped %d length outlier(s): %s", len(outliers), ", ".join(sorted(outlier_ids)))
    else:
        kept = list(records)

    if len(kept) < 2:
        raise PipelineError(f"Only {len(kept)} sequence(s) left after length filtering; at least 2 are needed")
    return LengthFilterResult(kept=kept, outliers=outliers, median=median, action=action)
