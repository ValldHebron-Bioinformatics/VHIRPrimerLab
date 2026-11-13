"""Reaction conditions (salt, dNTP, oligo concentration, Tm formula) and the thermodynamics cache."""

from __future__ import annotations

import csv
from pathlib import Path

import primer3
import pytest

from tests.synthetic import make_genome_family, write_fasta
from vhirprimerlab import structure, thermo
from vhirprimerlab.cli import main
from vhirprimerlab.degenerate import tm_range
from vhirprimerlab.design import _global_args
from vhirprimerlab.models import DesignParams
from vhirprimerlab.thermo import DEFAULT_THERMO, Thermo

SEQ = "ACGTTGCAAGCTTCGGATCC"
OTHER = "GGATCCGAAGCTTGCAACGT"          # SEQ reverse-complemented: the strongest possible partner


def test_defaults_are_what_primer3_uses_without_arguments() -> None:
    """A change in primer3's own defaults would silently move every Tm, so it is pinned here."""
    assert DEFAULT_THERMO.tm(SEQ) == primer3.calc_tm(SEQ)
    assert DEFAULT_THERMO.hairpin_tm(SEQ) == primer3.calc_hairpin(SEQ).tm
    assert DEFAULT_THERMO.homodimer_tm(SEQ) == primer3.calc_homodimer(SEQ).tm
    het = primer3.calc_heterodimer(SEQ, OTHER)
    assert DEFAULT_THERMO.heterodimer(SEQ, OTHER) == (het.tm, het.dg)


def test_salt_and_oligo_concentration_move_the_tm_in_the_expected_direction() -> None:
    base = DEFAULT_THERMO.tm(SEQ)
    assert Thermo(dv_conc=4.0).tm(SEQ) > base                 # more Mg2+ stabilises the duplex
    assert Thermo(mv_conc=20.0).tm(SEQ) < base                # less Na+ destabilises it
    assert Thermo(dna_conc=500.0).tm(SEQ) > base              # more oligo, higher Tm
    assert Thermo(tm_method="breslauer").tm(SEQ) != base
    assert Thermo(salt_method="owczarzy").tm(SEQ) != base


def test_the_conditions_reach_primer3_unchanged() -> None:
    t = Thermo(mv_conc=30.0, dv_conc=2.5, dntp_conc=0.4, dna_conc=200.0, salt_method="owczarzy")
    assert t.tm(SEQ) == primer3.calc_tm(SEQ, mv_conc=30.0, dv_conc=2.5, dntp_conc=0.4, dna_conc=200.0,
                                        tm_method="santalucia", salt_corrections_method="owczarzy")
    assert t.hairpin_tm(SEQ) == primer3.calc_hairpin(SEQ, mv_conc=30.0, dv_conc=2.5, dntp_conc=0.4,
                                                     dna_conc=200.0).tm


def test_the_structure_temperature_reaches_the_dimer_calculation() -> None:
    cold, warm = Thermo(structure_temp=37.0), Thermo(structure_temp=60.0)
    assert warm.heterodimer(SEQ, OTHER)[1] > cold.heterodimer(SEQ, OTHER)[1]   # weaker when warmer
    assert structure.dimer_dg(SEQ, OTHER, warm) == warm.heterodimer(SEQ, OTHER)[1]


def test_repeated_questions_are_answered_from_the_cache_in_either_order() -> None:
    t = Thermo(dv_conc=1.75)                       # a setting nothing else uses, so the cache is fresh
    t.tm(SEQ)
    before = thermo._tm.cache_info().hits
    t.tm(SEQ)
    assert thermo._tm.cache_info().hits == before + 1
    t.heterodimer(SEQ, OTHER)
    before = thermo._heterodimer.cache_info().hits
    assert t.heterodimer(OTHER, SEQ) == t.heterodimer(SEQ, OTHER)
    assert thermo._heterodimer.cache_info().hits == before + 2


@pytest.mark.parametrize("kwargs, text", [
    ({"mv_conc": -1.0}, "monovalent"),
    ({"dv_conc": -0.1}, "divalent"),
    ({"dntp_conc": -1.0}, "dNTP"),
    ({"dna_conc": 0.0}, "oligo"),
    ({"tm_method": "wallace"}, "Tm formula"),
    ({"salt_method": "magic"}, "salt correction"),
    ({"structure_temp": 200.0}, "temperature"),
    ({"max_dimer_dg": 1000.0}, "dimer"),
    ({"mv_conc": float("nan")}, "monovalent"),
])
def test_unusable_conditions_are_named(kwargs: dict, text: str) -> None:
    problems = Thermo(**kwargs).problems()
    assert problems and text in " ".join(problems)


def test_default_conditions_have_no_problems() -> None:
    assert DEFAULT_THERMO.problems() == []


def test_a_degenerate_primers_tm_range_follows_the_conditions() -> None:
    lo, hi = tm_range("ACGTTGCAAGCTTCGGATCR")
    lo2, hi2 = tm_range("ACGTTGCAAGCTTCGGATCR", Thermo(dv_conc=4.0))
    assert lo2 > lo and hi2 > hi


def test_primer3_design_gets_the_same_conditions() -> None:
    params = DesignParams(57, 60, 63, 2, 18, 25, 40, 60, 100, 300, 10,
                          Thermo(mv_conc=30.0, dv_conc=2.5, dntp_conc=0.4, dna_conc=200.0,
                                 tm_method="breslauer", salt_method="owczarzy"))
    args = _global_args(params)
    assert (args["PRIMER_SALT_MONOVALENT"], args["PRIMER_SALT_DIVALENT"], args["PRIMER_DNTP_CONC"],
            args["PRIMER_DNA_CONC"]) == (30.0, 2.5, 0.4, 200.0)
    assert args["PRIMER_TM_FORMULA"] == 0 and args["PRIMER_SALT_CORRECTIONS"] == 2
    default = _global_args(DesignParams(57, 60, 63, 2, 18, 25, 40, 60, 100, 300, 10))
    assert (default["PRIMER_TM_FORMULA"], default["PRIMER_SALT_CORRECTIONS"]) == (1, 1)


# --- end to end ---------------------------------------------------------------------

def _rows(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _tiles(tmp_path: Path, name: str, *flags: str) -> Path:
    aln = write_fasta(tmp_path / "aln.fasta", make_genome_family())
    out = tmp_path / name
    assert main(["--aligned", str(aln), "-n", "5", "-o", str(out), *flags]) in (0, 3)
    return out


def test_stating_the_default_conditions_changes_nothing(tmp_path: Path) -> None:
    plain = _tiles(tmp_path, "plain")
    stated = _tiles(tmp_path, "stated", "--mv-conc", "50", "--dv-conc", "1.5", "--dntp-conc", "0.6",
                    "--dna-conc", "50", "--tm-method", "santalucia", "--salt-correction", "santalucia",
                    "--structure-temp", "37", "--max-dimer-dg", "-9")
    assert (plain / "primers.tsv").read_text() == (stated / "primers.tsv").read_text()
    assert (plain / "amplicons.tsv").read_text() == (stated / "amplicons.tsv").read_text()


def test_the_conditions_decide_the_reported_tm(tmp_path: Path) -> None:
    plain = _tiles(tmp_path, "plain")
    salty = _tiles(tmp_path, "salty", "--dv-conc", "4", "--tm-max", "70")
    t = Thermo(dv_conc=4.0)
    checked = 0
    for row in _rows(salty / "amplicons.tsv"):
        seq = row["forward_seq"]
        if not set(seq) - set("ACGT"):
            assert float(row["forward_tm"]) == pytest.approx(t.tm(seq), abs=0.05)
            checked += 1
    assert checked
    assert ([r["forward_tm"] for r in _rows(plain / "amplicons.tsv")]
            != [r["forward_tm"] for r in _rows(salty / "amplicons.tsv")])
    assert "dv_conc: 4.0" in (salty / "report.txt").read_text()


@pytest.mark.parametrize("flag, value", [
    ("--dv-conc", "-1"), ("--mv-conc", "-5"), ("--dna-conc", "0"), ("--structure-temp", "300"),
    ("--max-dimer-dg", "5"),
])
def test_a_bad_condition_is_refused(tmp_path: Path, flag: str, value: str) -> None:
    aln = write_fasta(tmp_path / "aln.fasta", make_genome_family())
    assert main(["--aligned", str(aln), "-n", "5", "-o", str(tmp_path / "o"), flag, value]) == 1


def test_check_mode_reports_tm_and_judges_dimers_under_the_given_conditions(tmp_path: Path) -> None:
    designed = _tiles(tmp_path, "design")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("".join(f"a{r['amplicon']}\t{r['forward_seq']}\t{r['reverse_seq']}\n"
                             for r in _rows(designed / "amplicons.tsv")))

    def check(name: str, *flags: str) -> list[dict[str, str]]:
        out = tmp_path / name
        assert main(["check", "--aligned", str(tmp_path / "aln.fasta"), "--primers", str(pairs),
                     "-o", str(out), *flags]) in (0, 3)
        return _rows(out / "check.tsv")

    plain, salty = check("c1"), check("c2", "--dv-conc", "4")
    assert [r["forward_tm"] for r in plain] != [r["forward_tm"] for r in salty]
    # a limit above every pair's ΔG calls every pair a dimer; the default calls none
    assert "forms a dimer" not in (tmp_path / "c1" / "report.txt").read_text()
    check("c3", "--max-dimer-dg", "-0.1")
    assert "forms a dimer" in (tmp_path / "c3" / "report.txt").read_text()
