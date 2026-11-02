"""Hairpin and primer-dimer checks that cover every variant of a degenerate primer.

A degenerate primer is ordered as a mix of plain oligos (CACR... is CACA... plus
CACG...), so every check runs on each variant, and on each pair of variants when
two oligos of the mix could stick to each other.
"""

from __future__ import annotations

import itertools

import primer3

from vhirprimerlab.scoring import IUPAC

MAX_STRUCTURE_TM = 47.0  # primer3's default limit for hairpin and self-dimer Tm
MAX_END_GC = 3  # at most 3 G/C among the last 5 bases (3' end)
DIMER_DG = -9000.0  # cal/mol; two primers in one pool below this form a dimer


def variants(primer: str) -> list[str]:
    """Every plain (ACGT) oligo a possibly degenerate primer stands for."""
    return ["".join(v) for v in itertools.product(*(sorted(IUPAC[b]) for b in primer))]


def end_gc(primer: str) -> int:
    """G/C count in the last 5 bases of the worst variant."""
    return sum(bool(IUPAC[b] & {"G", "C"}) for b in primer[-5:])


def structure_ok(primer: str) -> bool:
    """True when no variant folds into a hairpin or sticks to itself or to another variant,
    and no variant has more than MAX_END_GC G/C in its last 5 bases. primer is 5'->3'."""
    if end_gc(primer) > MAX_END_GC:
        return False
    vs = variants(primer)
    if any(primer3.calc_hairpin(v).tm > MAX_STRUCTURE_TM for v in vs):
        return False
    for i, a in enumerate(vs):
        if primer3.calc_homodimer(a).tm > MAX_STRUCTURE_TM:
            return False
        if any(primer3.calc_heterodimer(a, b).tm > MAX_STRUCTURE_TM for b in vs[i + 1:]):
            return False
    return True


def dimer_dg(a: str, b: str) -> float:
    """Lowest heterodimer ΔG (cal/mol) between any variant of a and any variant of b."""
    return min(primer3.calc_heterodimer(x, y).dg for x in variants(a) for y in variants(b))
