"""Hairpin and primer-dimer checks that cover every variant of a degenerate primer.

A degenerate primer is ordered as a mix of plain oligos (CACR... is CACA... plus
CACG...), so every check runs on each variant, and on each pair of variants when
two oligos of the mix could stick to each other.
"""

from __future__ import annotations

import itertools

from vhirprimerlab.scoring import IUPAC
from vhirprimerlab.thermo import DEFAULT_THERMO, DIMER_DG, Thermo

MAX_STRUCTURE_TM = 47.0  # primer3's default limit for hairpin and self-dimer Tm
MAX_END_GC = 3  # at most 3 G/C among the last 5 bases (3' end)


def variants(primer: str) -> list[str]:
    """Every plain (ACGT) oligo a possibly degenerate primer stands for."""
    return ["".join(v) for v in itertools.product(*(sorted(IUPAC[b]) for b in primer))]


def end_gc(primer: str) -> int:
    """G/C count in the last 5 bases of the worst variant."""
    return sum(bool(IUPAC[b] & {"G", "C"}) for b in primer[-5:])


def structure_ok(primer: str, thermo: Thermo = DEFAULT_THERMO) -> bool:
    """True when no variant folds into a hairpin or sticks to itself or to another variant,
    and no variant has more than MAX_END_GC G/C in its last 5 bases. primer is 5'->3'."""
    if end_gc(primer) > MAX_END_GC:
        return False
    vs = variants(primer)
    if any(thermo.hairpin_tm(v) > MAX_STRUCTURE_TM for v in vs):
        return False
    for i, a in enumerate(vs):
        if thermo.homodimer_tm(a) > MAX_STRUCTURE_TM:
            return False
        if any(thermo.heterodimer(a, b)[0] > MAX_STRUCTURE_TM for b in vs[i + 1:]):
            return False
    return True


def structure_report(primer: str, thermo: Thermo = DEFAULT_THERMO) -> list[str]:
    """What structure_ok rejects, with the worst variant's value; empty when it passes."""
    out: list[str] = []
    n = end_gc(primer)
    if n > MAX_END_GC:
        out.append(f"{n} G/C in the last 5 bases (max {MAX_END_GC})")
    vs = variants(primer)
    worst = (
        ("hairpin", max(thermo.hairpin_tm(v) for v in vs)),
        ("self-dimer", max(thermo.homodimer_tm(v) for v in vs)),
        ("dimer between its variants", max(
            (thermo.heterodimer(a, b)[0] for i, a in enumerate(vs) for b in vs[i + 1:]),
            default=float("-inf"),
        )),
    )
    out += [f"{what} Tm {tm:.1f} °C (max {MAX_STRUCTURE_TM:g})" for what, tm in worst
            if tm > MAX_STRUCTURE_TM]
    return out


def dimer_dg(a: str, b: str, thermo: Thermo = DEFAULT_THERMO) -> float:
    """Lowest heterodimer ΔG (cal/mol) between any variant of a and any variant of b."""
    return min(thermo.heterodimer(x, y)[1] for x in variants(a) for y in variants(b))
