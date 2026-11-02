from __future__ import annotations

from types import SimpleNamespace

from vhirprimerlab import structure
from vhirprimerlab.structure import dimer_dg, end_gc, structure_ok, variants


def test_variants_expand_iupac() -> None:
    assert variants("ACGT") == ["ACGT"]
    assert sorted(variants("ARY")) == ["AAC", "AAT", "AGC", "AGT"]


def test_end_gc_counts_worst_variant() -> None:
    assert end_gc("AAAAAATTGCA") == 2
    assert end_gc("AAAAASATTA") == 1  # S = C or G
    assert end_gc("AAAAATCRSC") == 4  # R can be G


def test_structure_ok_plain_primer() -> None:
    assert structure_ok("CACAGCTCTAGGATCAGAAAGAGT")
    assert not structure_ok("GAATTCGCGAATTCGCGAATTC")  # palindrome: strong self-dimer


def test_structure_ok_checks_every_variant(monkeypatch) -> None:
    seen: list[str] = []

    def hairpin(seq):
        seen.append(seq)
        return SimpleNamespace(tm=60.0 if "GGG" in seq else 0.0)

    monkeypatch.setattr(structure.primer3, "calc_hairpin", hairpin)
    # only the G variant of R "folds": the primer must be rejected
    assert not structure_ok("CACAGCTCTAGGRTCAGAAAGAGT")
    assert {"CACAGCTCTAGGATCAGAAAGAGT", "CACAGCTCTAGGGTCAGAAAGAGT"} <= set(seen)


def test_structure_ok_checks_variant_pairs(monkeypatch) -> None:
    monkeypatch.setattr(structure.primer3, "calc_hairpin", lambda s: SimpleNamespace(tm=0.0))
    monkeypatch.setattr(structure.primer3, "calc_homodimer", lambda s: SimpleNamespace(tm=0.0))
    calls: list[tuple[str, str]] = []

    def hetero(a, b):
        calls.append((a, b))
        return SimpleNamespace(tm=50.0)

    monkeypatch.setattr(structure.primer3, "calc_heterodimer", hetero)
    assert structure_ok("ACGTACGTACGTACGTAT")  # one variant: no pairs to check
    assert not structure_ok("ACGTACGTRCGTACGTAT")
    assert calls == [("ACGTACGTACGTACGTAT", "ACGTACGTGCGTACGTAT")]


def test_dimer_dg_takes_worst_variant_pair(monkeypatch) -> None:
    def hetero(a, b):
        return SimpleNamespace(dg=-12000.0 if a.startswith("G") and b.startswith("T") else -3000.0)

    monkeypatch.setattr(structure.primer3, "calc_heterodimer", hetero)
    assert dimer_dg("RAAA", "YCCC") == -12000.0
    assert dimer_dg("AAAA", "YCCC") == -3000.0
