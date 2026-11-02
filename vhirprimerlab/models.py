"""Shared data types for every pipeline stage."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


class PipelineError(Exception):
    """A user-facing error: bad input, bad parameters or a failed external tool."""


@dataclass(frozen=True)
class Record:
    id: str
    seq: str


@dataclass
class Alignment:
    ids: list[str]
    rows: list[str]

    @property
    def length(self) -> int:
        return len(self.rows[0]) if self.rows else 0

    def row(self, id: str) -> str:
        try:
            return self.rows[self.ids.index(id)]
        except ValueError:
            raise PipelineError(f"Sequence '{id}' not found in alignment") from None


@dataclass(frozen=True)
class Region:
    name: str
    start: int
    end: int


@dataclass
class LengthFilterResult:
    kept: list[Record]
    outliers: list[tuple[str, int]]
    median: float
    action: str


@dataclass(frozen=True)
class ColumnStats:
    majority_base: str
    majority_freq: float
    gap_frac: float
    n_ambiguous: int


@dataclass
class Consensus:
    seq: str
    cons_to_col: list[int]


@dataclass(frozen=True)
class ScoringParams:
    max_mismatches: int
    window: int


@dataclass(frozen=True)
class MatchResult:
    ok: bool
    mismatches: int
    three_prime_mismatches: int
    reason: str | None


@dataclass
class PairCoverage:
    covered: list[str]
    not_covered: dict[str, str]
    coverage_pct: float


@dataclass(frozen=True)
class Primer:
    seq: str
    start: int  # 0-based leftmost consensus position (for both primers)
    length: int
    tm: float
    gc: float


@dataclass(frozen=True)
class CandidatePair:
    forward: Primer
    reverse: Primer
    penalty: float


@dataclass(frozen=True)
class DesignParams:
    tm_min: float
    tm_opt: float
    tm_max: float
    max_tm_diff: float
    len_min: int
    len_max: int
    gc_min: float
    gc_max: float
    amplicon_min: int
    amplicon_max: int
    candidates: int


@dataclass
class DesignResult:
    pairs: list[CandidatePair]
    explain: dict[str, str]


@dataclass
class ScoredPair:
    pair: CandidatePair
    fwd_seq: str  # 5'->3', possibly degenerate
    rev_seq: str  # 5'->3', possibly degenerate
    fwd_cols: list[int]
    rev_cols: list[int]
    coverage: PairCoverage
    degenerate: bool
    fwd_tm_range: tuple[float, float] | None
    rev_tm_range: tuple[float, float] | None


@dataclass
class RankResult:
    pairs: list[ScoredPair]
    below_threshold: bool


@dataclass
class RegionResult:
    region: Region
    included: tuple[int, int]
    n_candidates: int
    n_passing: int
    rank: RankResult
    explain: dict[str, str]
    failures: list[tuple[str, int]]


@dataclass
class RunReport:
    command: str
    params: dict[str, object]
    n_input: int
    median_len: float
    min_len: int
    max_len: int
    length: LengthFilterResult
    reference_id: str
    aln_length: int
    cons_length: int
    n_conserved: int
    mean_conservation: float
    min_coverage: float
    regions: list[RegionResult]
    log_lines: list[str]


@dataclass(frozen=True)
class Assay:
    """A primer pair to check, both oligos written 5'->3' (upper-case, U as T)."""

    name: str
    forward: str
    reverse: str


@dataclass(frozen=True)
class Location:
    """Best binding site of a primer on the consensus (0-based start)."""

    start: int
    mismatches: int


@dataclass
class CheckResult:
    assay: Assay
    located: bool
    note: str
    fwd: Location | None
    rev: Location | None
    fwd_cols: list[int]
    rev_cols: list[int]
    coverage: PairCoverage
    fwd_tm: float | None
    rev_tm: float | None
    fwd_tm_range: tuple[float, float] | None
    rev_tm_range: tuple[float, float] | None
    fwd_gc: float
    rev_gc: float
    degenerate: bool
    below_threshold: bool


@dataclass
class CheckReport:
    command: str
    params: dict[str, object]
    n_input: int
    median_len: float
    min_len: int
    max_len: int
    length: LengthFilterResult
    reference_id: str
    aln_length: int
    cons_length: int
    n_conserved: int
    mean_conservation: float
    min_coverage: float
    results: list[CheckResult]
    log_lines: list[str]


@dataclass(frozen=True)
class TilingParams:
    """Settings for tiling the whole genome with N overlapping amplicons."""

    n_amplicons: int
    overlap: int
    tm_min: float
    tm_opt: float
    tm_max: float
    max_tm_diff: float
    len_min: int
    len_max: int
    gc_min: float
    gc_max: float
    max_mismatches: int
    window: int
    degenerate: bool
    max_degeneracy: int
    min_coverage: float

    @property
    def scoring(self) -> ScoringParams:
        return ScoringParams(self.max_mismatches, self.window)


@dataclass(frozen=True)
class Geometry:
    """Amplicon size limits for one tiling attempt (consensus positions)."""

    target_len: int
    len_min: int
    len_max: int
    tolerance: float
    end_slack: int
    max_outer_overlap: int


@dataclass
class PrimerCandidate:
    """One primer considered for tiling, located on the consensus.

    site is the oligo as it reads on the forward strand (the reverse primer
    reverse-complemented); seq is the primer 5'->3' as ordered. base_seq is seq
    before any degenerate bases were added. Bit i of mask is set when the primer
    matches aln.ids[i].
    """

    strand: Literal["+", "-"]
    start: int
    length: int
    site: str
    base_seq: str
    seq: str
    tm: float
    tm_range: tuple[float, float] | None
    gc: float
    penalty: float
    mask: int = 0
    cols: list[int] = field(default_factory=list)

    @property
    def end(self) -> int:
        return self.start + self.length - 1


@dataclass
class Amplicon:
    index: int  # 1-based
    pool: int
    fwd: PrimerCandidate
    rev: PrimerCandidate
    coverage: PairCoverage
    below_threshold: bool
    overlap_next: int | None  # insert overlap with the next amplicon (consensus bp)


@dataclass
class TilingResult:
    amplicons: list[Amplicon]
    geometry: Geometry
    n_fwd: int
    n_rev: int
    dimer_warnings: list[str]
    uncovered_ends: tuple[int, int]  # reference bases outside the first / last insert


@dataclass
class TilingReport:
    command: str
    params: dict[str, object]
    n_input: int
    median_len: float
    min_len: int
    max_len: int
    length: LengthFilterResult
    reference_id: str
    aln_length: int
    cons_length: int
    n_conserved: int
    mean_conservation: float
    min_coverage: float
    result: TilingResult
    log_lines: list[str]
