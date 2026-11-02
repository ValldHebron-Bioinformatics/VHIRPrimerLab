"""Conversion between reference positions, alignment columns and consensus positions.

Reference positions are 1-based; alignment columns and consensus positions are 0-based.
"""

from __future__ import annotations

import bisect
import logging
from typing import Literal

from vhirprimerlab.models import Alignment, Consensus, PipelineError, Region

log = logging.getLogger(__name__)

Side = Literal["start", "end"]


def _snap(sorted_cols: list[int], col: int, side: Side, what: str) -> int:
    """Index in sorted_cols of col, or of the nearest column right (start) / left (end)."""
    if side == "start":
        i = bisect.bisect_left(sorted_cols, col)
        if i == len(sorted_cols):
            raise PipelineError(f"No {what} position at or after alignment column {col + 1}")
    else:
        i = bisect.bisect_right(sorted_cols, col) - 1
        if i < 0:
            raise PipelineError(f"No {what} position at or before alignment column {col + 1}")
    return i


class CoordMap:
    def __init__(self, aln: Alignment, reference_id: str, consensus: Consensus) -> None:
        if reference_id not in aln.ids:
            raise PipelineError(
                f"Reference sequence '{reference_id}' not found among the (filtered) sequences"
            )
        self.reference_id = reference_id
        ref_row = aln.row(reference_id)
        self._ref_cols = [c for c, ch in enumerate(ref_row) if ch != "-"]
        self._cons_cols = consensus.cons_to_col
        self.ref_length = len(self._ref_cols)

    def ref_to_col(self, pos: int) -> int:
        if not 1 <= pos <= self.ref_length:
            raise PipelineError(
                f"Position {pos} is outside the reference '{self.reference_id}' (1-{self.ref_length})"
            )
        return self._ref_cols[pos - 1]

    def col_to_ref(self, col: int, side: Side) -> int:
        return _snap(self._ref_cols, col, side, "reference") + 1

    def col_to_cons(self, col: int, side: Side) -> int:
        return _snap(self._cons_cols, col, side, "consensus")

    def validate_region(self, region: Region, amplicon_min: int) -> None:
        if region.start < 1 or region.end > self.ref_length:
            raise PipelineError(
                f"Region '{region.name}' ({region.start}-{region.end}) is outside the reference "
                f"'{self.reference_id}' (1-{self.ref_length})"
            )
        if region.start >= region.end:
            raise PipelineError(f"Region '{region.name}': start must be smaller than end")
        if region.end - region.start + 1 < amplicon_min:
            raise PipelineError(
                f"Region '{region.name}' is {region.end - region.start + 1} bp, shorter than "
                f"--amplicon-min ({amplicon_min})"
            )

    def region_to_included(self, region: Region) -> tuple[int, int]:
        """primer3 SEQUENCE_INCLUDED_REGION (0-based start, length) on the consensus."""
        start_col = self.ref_to_col(region.start)
        end_col = self.ref_to_col(region.end)
        cons_start = self.col_to_cons(start_col, "start")
        cons_end = self.col_to_cons(end_col, "end")
        if cons_end < cons_start:
            raise PipelineError(f"Region '{region.name}' contains no consensus positions")
        if self._cons_cols[cons_start] != start_col or self._cons_cols[cons_end] != end_col:
            log.info(
                "Region '%s': edges fall on gap-majority columns; snapped inward to consensus "
                "positions %d-%d",
                region.name, cons_start + 1, cons_end + 1,
            )
        return cons_start, cons_end - cons_start + 1
