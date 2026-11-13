"""Where else a primer scheme can bind: secondary sites and the products they could make.

A primer binds a sequence at a site when it follows the same rule as coverage scoring: at most
``max_mismatches`` mismatches and none in the last ``window`` bases of its 3' end. Every such site
other than the primer's designed one is a secondary site. Two primers of one pool whose sites face
each other within ``max_product`` bases could amplify an unintended product when at least one of
the two sites is secondary. Sequences are scanned on both strands; ``N`` and other ambiguity codes
in a sequence never form a site.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass
from typing import Iterable, Literal

from vhirprimerlab.models import Alignment, ScoringParams
from vhirprimerlab.scoring import IUPAC

MIN_OLIGO = 15  # shorter than this and the 3' window is most of the primer
MAX_FIND_VARIANTS = 16  # up to this many seeds are looked up with str.find, more with one regex
MAX_SEED_VARIANTS = 256  # a degenerate 3' end with more variants than this is scanned base by base
MAX_REPORTED_PRODUCTS = 50

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class Oligo:
    """A primer as ordered (5'->3') and where it was designed to bind.

    ``strand`` is "+" for a forward primer and "-" for a reverse one; ``first_col`` and
    ``last_col`` are the alignment columns its designed site spans.
    """

    name: str
    pool: int
    seq: str
    strand: Literal["+", "-"]
    first_col: int
    last_col: int


@dataclass(frozen=True)
class PrimerSpecificity:
    name: str
    pool: int
    off_target_sites: int  # secondary sites in the input sequences
    sequences: int  # input sequences that have at least one
    fewest_mismatches: int | None  # of the best secondary site
    background_sites: int = 0  # sites in the background sequences


@dataclass(frozen=True)
class OffTargetProduct:
    pool: int
    fwd: str  # the primer that extends to the right
    rev: str  # the primer that extends to the left
    fwd_secondary: bool  # the fwd primer binds a secondary site
    rev_secondary: bool
    sequences: int
    length_min: int
    length_max: int
    example: str  # a sequence that has it
    background: bool  # found in a background sequence


@dataclass
class SpecificityReport:
    primers: list[PrimerSpecificity]
    products: list[OffTargetProduct]
    n_products: int  # before the list was cut to MAX_REPORTED_PRODUCTS
    n_background: int = 0  # background sequences scanned

    @property
    def clean(self) -> bool:
        return not self.products and all(
            p.off_target_sites == 0 and p.background_sites == 0 for p in self.primers)


@dataclass(frozen=True)
class _Hit:
    oligo: int  # index into the oligos
    strand: str  # "+": the primer extends to the right on the sequence as written
    start: int  # 0-based, on the sequence as written
    end: int  # exclusive
    mismatches: int
    secondary: bool


class _Scanner:
    """One oligo's pattern: per-position allowed bases and the regex of its 3' seed."""

    def __init__(self, seq: str, window: int, max_mismatches: int) -> None:
        self.allowed = [IUPAC[c] for c in seq]
        self.length = len(seq)
        self.max_mismatches = max_mismatches
        self.seed_len = min(window, len(seq))
        seeds = ["".join(v) for v in itertools.islice(
            itertools.product(*(sorted(b) for b in self.allowed[-self.seed_len:])),
            MAX_SEED_VARIANTS + 1)]
        self.seeds = seeds if len(seeds) <= MAX_FIND_VARIANTS else []
        self.pattern = (re.compile("(?=(?:" + "|".join(seeds) + "))")
                        if MAX_FIND_VARIANTS < len(seeds) <= MAX_SEED_VARIANTS else None)

    def _seed_positions(self, text: str) -> Iterable[int]:
        """Where the 3' seed occurs in ``text``, as the start of the oligo it would belong to."""
        offset = self.length - self.seed_len
        if self.seeds:
            for seed in self.seeds:
                at = text.find(seed, offset)
                while at != -1:
                    yield at - offset
                    at = text.find(seed, at + 1)
        elif self.pattern is not None:
            for m in self.pattern.finditer(text, offset):
                yield m.start() - offset
        else:
            yield from range(0, len(text) - self.length + 1)

    def scan(self, text: str) -> list[tuple[int, int]]:
        """(start, mismatches) of every site on ``text``, as written."""
        length, seed, found = self.length, self.seed_len, []
        n = len(text)
        for i in self._seed_positions(text):
            if i < 0 or i + length > n:
                continue
            mismatches = 0
            for allowed, base in zip(self.allowed, text[i:i + length]):
                if base not in allowed or base not in "ACGT":
                    mismatches += 1
                    if mismatches > self.max_mismatches:
                        break
            else:
                # the seed is the 3' window; a base outside ACGT there is no site either
                if all(b in "ACGT" and b in a for a, b in
                       zip(self.allowed[-seed:], text[i + length - seed:i + length])):
                    found.append((i, mismatches))
        return found


def _ungapped_interval(row: str, first_col: int, last_col: int) -> tuple[int, int]:
    """The half-open interval an alignment column range covers once the gaps are removed."""
    return (first_col - row.count("-", 0, first_col),
            last_col + 1 - row.count("-", 0, last_col + 1))


def check_specificity(
    oligos: list[Oligo],
    aln: Alignment,
    params: ScoringParams,
    max_product: int,
    background: list[tuple[str, str]] | None = None,
) -> SpecificityReport:
    """Scan every oligo against every sequence of ``aln`` (and the background) on both strands.

    ``max_product`` is the longest product worth reporting; use the longest amplicon allowed.
    """
    for o in oligos:
        if len(o.seq) < MIN_OLIGO or any(c not in IUPAC for c in o.seq):
            raise ValueError(f"primer {o.name!r} is not a usable oligo of at least {MIN_OLIGO} bases")
    scanners = [_Scanner(o.seq, params.window, params.max_mismatches) for o in oligos]

    # (id, text, is_background, designed intervals per oligo or None)
    targets: list[tuple[str, str, bool, list[tuple[int, int]] | None]] = []
    for sid, row in zip(aln.ids, aln.rows):
        targets.append((sid, row.replace("-", "").upper(), False,
                        [_ungapped_interval(row, o.first_col, o.last_col) for o in oligos]))
    for sid, text in background or []:
        targets.append((sid, text.replace("-", "").upper(), True, None))

    cache: dict[str, list[tuple[int, int, int, str]]] = {}  # identical texts are scanned once
    hits_by_target: list[list[_Hit]] = []
    for sid, text, _, designed in targets:
        if text not in cache:
            rc, n, found = reverse_complement(text), len(text), []
            for k, scanner in enumerate(scanners):
                found += [(k, i, mm, "+") for i, mm in scanner.scan(text)]
                found += [(k, n - j - scanner.length, mm, "-") for j, mm in scanner.scan(rc)]
            cache[text] = found
        hits: list[_Hit] = []
        for k, start, mm, strand in cache[text]:
            end = start + scanners[k].length
            own = False
            if designed is not None and strand == oligos[k].strand:
                d_start, d_end = designed[k]
                own = start < d_end and d_start < end
            hits.append(_Hit(k, strand, start, end, mm, not own))
        hits_by_target.append(hits)

    counts: list[list] = [[0, set(), None, 0] for _ in oligos]  # sites, sequences, fewest, background
    products: dict[tuple, dict] = {}
    for (sid, _, is_background, _), hits in zip(targets, hits_by_target):
        for h in hits:
            if not h.secondary:
                continue
            c = counts[h.oligo]
            if is_background:
                c[3] += 1
            else:
                c[0] += 1
                c[1].add(sid)
                c[2] = h.mismatches if c[2] is None else min(c[2], h.mismatches)
        rights = [h for h in hits if h.strand == "+"]
        lefts = [h for h in hits if h.strand == "-"]
        for a, b in itertools.product(rights, lefts):
            if not (a.secondary or b.secondary) or oligos[a.oligo].pool != oligos[b.oligo].pool:
                continue
            length = b.end - a.start
            if b.end <= a.end or not MIN_OLIGO <= length <= max_product:
                continue
            key = (oligos[a.oligo].pool, a.oligo, b.oligo, a.secondary, b.secondary, is_background)
            entry = products.setdefault(key, {"seqs": set(), "min": length, "max": length, "ex": sid})
            entry["seqs"].add(sid)
            entry["min"], entry["max"] = min(entry["min"], length), max(entry["max"], length)

    primers = [PrimerSpecificity(o.name, o.pool, c[0], len(c[1]), c[2], c[3])
               for o, c in zip(oligos, counts)]
    ordered = sorted(products.items(), key=lambda kv: (-len(kv[1]["seqs"]), kv[0]))
    reported = [OffTargetProduct(pool, oligos[a].name, oligos[b].name, fa, fb, len(e["seqs"]),
                                 e["min"], e["max"], e["ex"], bg)
                for (pool, a, b, fa, fb, bg), e in ordered[:MAX_REPORTED_PRODUCTS]]
    return SpecificityReport(primers, reported, len(ordered), len(background or []))
