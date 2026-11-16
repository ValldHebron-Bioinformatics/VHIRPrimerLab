"""The reaction conditions every melting-temperature and structure calculation is made under.

One ``Thermo`` carries them, so a primer's Tm, its hairpins and its dimers are all computed for
the same buffer. The defaults are primer3-py's own (50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM
oligo, SantaLucia), written out here so that a different primer3-py version cannot change them
silently. Calls to primer3 are cached: the same primer is asked about many times while tiling.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import primer3

TM_METHODS = ("santalucia", "breslauer")
SALT_METHODS = ("santalucia", "owczarzy", "schildkraut")
# primer3's own numbering of the same choices (PRIMER_TM_FORMULA, PRIMER_SALT_CORRECTIONS)
TM_FORMULA = {"breslauer": 0, "santalucia": 1}
SALT_CORRECTIONS = {"schildkraut": 0, "santalucia": 1, "owczarzy": 2}
DIMER_DG = -9000.0  # cal/mol; two primers in one pool below this form a dimer
CACHE_SIZE = 200_000


@dataclass(frozen=True)
class Thermo:
    mv_conc: float = 50.0  # mM monovalent cations (Na⁺, K⁺)
    dv_conc: float = 1.5  # mM divalent cations (Mg²⁺)
    dntp_conc: float = 0.6  # mM of each dNTP, which bind Mg²⁺
    dna_conc: float = 50.0  # nM oligo
    tm_method: str = "santalucia"
    salt_method: str = "santalucia"
    structure_temp: float = 37.0  # °C at which hairpin and dimer ΔG are evaluated
    max_dimer_dg: float = DIMER_DG  # cal/mol; two primers of a pool below this form a dimer

    def problems(self) -> list[str]:
        """What is unusable in these conditions, in words for the user; empty when fine."""
        out = []
        for value, label, ok in (
            (self.mv_conc, "monovalent salt concentration", lambda v: v >= 0),
            (self.dv_conc, "divalent salt concentration", lambda v: v >= 0),
            (self.dntp_conc, "dNTP concentration", lambda v: v >= 0),
            (self.dna_conc, "oligo concentration", lambda v: v > 0),
            (self.structure_temp, "structure temperature", lambda v: 0 <= v <= 100),
            (self.max_dimer_dg, "dimer ΔG limit", lambda v: v < 0),
        ):
            if not math.isfinite(value) or not ok(value):
                out.append(f"the {label} is out of range ({value:g})")
        if self.tm_method not in TM_METHODS:
            out.append(f"the Tm formula must be one of {', '.join(TM_METHODS)}")
        if self.salt_method not in SALT_METHODS:
            out.append(f"the salt correction must be one of {', '.join(SALT_METHODS)}")
        return out

    def tm(self, seq: str) -> float:
        return _tm(self, seq)

    def hairpin_tm(self, seq: str) -> float:
        return _hairpin(self, seq)

    def homodimer_tm(self, seq: str) -> float:
        return _homodimer(self, seq)

    def heterodimer(self, a: str, b: str) -> tuple[float, float]:
        """(Tm in °C, ΔG in cal/mol) of the strongest dimer of two oligos."""
        return _heterodimer(self, *sorted((a, b)))

    def _salts(self) -> dict[str, float]:
        return {"mv_conc": self.mv_conc, "dv_conc": self.dv_conc, "dntp_conc": self.dntp_conc,
                "dna_conc": self.dna_conc}


DEFAULT_THERMO = Thermo()


@lru_cache(maxsize=CACHE_SIZE)
def _tm(t: Thermo, seq: str) -> float:
    return primer3.calc_tm(seq, **t._salts(), tm_method=t.tm_method,
                           salt_corrections_method=t.salt_method)


@lru_cache(maxsize=CACHE_SIZE)
def _hairpin(t: Thermo, seq: str) -> float:
    return primer3.calc_hairpin(seq, **t._salts(), temp_c=t.structure_temp).tm


@lru_cache(maxsize=CACHE_SIZE)
def _homodimer(t: Thermo, seq: str) -> float:
    return primer3.calc_homodimer(seq, **t._salts(), temp_c=t.structure_temp).tm


@lru_cache(maxsize=CACHE_SIZE)
def _heterodimer(t: Thermo, a: str, b: str) -> tuple[float, float]:
    res = primer3.calc_heterodimer(a, b, **t._salts(), temp_c=t.structure_temp)
    return res.tm, res.dg


def clear_caches() -> None:
    """Forget every cached answer (for tests that replace primer3 with a stand-in)."""
    for cached in (_tm, _hairpin, _homodimer, _heterodimer):
        cached.cache_clear()
