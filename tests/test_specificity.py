"""Secondary binding sites and unintended products of a primer scheme."""

from __future__ import annotations

import random

import pytest

from vhirprimerlab.models import Alignment, ScoringParams
from vhirprimerlab.specificity import Oligo, check_specificity, reverse_complement

P = ScoringParams(max_mismatches=2, window=5)
FWD = "ACGTTGCAAGCTTCGGATCC"      # 20-mer
REV = "TTGACCGTAGGCATCGATGC"      # 20-mer, 5'->3'


def genome(length: int = 1500, seed: int = 7, avoid: tuple[str, ...] = ()) -> str:
    """A random sequence that holds none of the given oligos or their reverse complements."""
    rng = random.Random(seed)
    while True:
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        if not any(o in seq or reverse_complement(o) in seq for o in avoid):
            return seq


def place(base: str, at: int, text: str) -> str:
    return base[:at] + text + base[at + len(text):]


def scheme(row: str) -> tuple[list[Oligo], Alignment]:
    """One amplicon, forward primer at 100 and reverse primer ending at 700, in a gapless row."""
    row = place(place(row, 100, FWD), 680, reverse_complement(REV))
    aln = Alignment(["s1"], [row])
    return [Oligo("A_LEFT", 1, FWD, "+", 100, 119), Oligo("A_RIGHT", 1, REV, "-", 680, 699)], aln


def clean() -> tuple[list[Oligo], Alignment]:
    return scheme(genome(avoid=(FWD, REV)))


def sites(rep, name: str) -> int:
    return {p.name: p for p in rep.primers}[name].off_target_sites


def test_reverse_complement_handles_ambiguity_codes() -> None:
    assert reverse_complement("ACGTRYN") == "NRYACGT"


def test_a_clean_scheme_has_no_secondary_sites() -> None:
    oligos, aln = clean()
    rep = check_specificity(oligos, aln, P, max_product=900)
    assert all(p.off_target_sites == 0 for p in rep.primers) and rep.products == []


def test_a_second_exact_site_is_reported_for_that_primer_only() -> None:
    oligos, aln = clean()
    aln = Alignment(aln.ids, [place(aln.rows[0], 1100, FWD)])     # FWD again, far from its own site
    rep = check_specificity(oligos, aln, P, max_product=900)
    left = {p.name: p for p in rep.primers}["A_LEFT"]
    assert left.off_target_sites == 1 and left.sequences == 1 and left.fewest_mismatches == 0
    assert sites(rep, "A_RIGHT") == 0


def test_a_site_with_two_mismatches_counts_and_one_with_three_does_not() -> None:
    oligos, aln = clean()
    two = place(aln.rows[0], 1100, FWD[:3] + "T" + FWD[4:9] + "A" + FWD[10:])
    assert sites(check_specificity(oligos, Alignment(["s1"], [two]), P, max_product=900), "A_LEFT") == 1
    three = place(aln.rows[0], 1100, "TTT" + FWD[3:9] + "A" + FWD[10:])
    assert sites(check_specificity(oligos, Alignment(["s1"], [three]), P, max_product=900), "A_LEFT") == 0


def test_a_mismatch_in_the_3_prime_window_stops_a_site_from_counting() -> None:
    oligos, aln = clean()
    row = place(aln.rows[0], 1100, FWD[:-2] + ("A" if FWD[-2] != "A" else "C") + FWD[-1])
    assert sites(check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900), "A_LEFT") == 0


def test_a_site_on_the_other_strand_is_found() -> None:
    oligos, aln = clean()
    row = place(aln.rows[0], 1100, reverse_complement(FWD))
    assert sites(check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900), "A_LEFT") == 1


def test_a_degenerate_primer_finds_every_variant() -> None:
    oligos, aln = clean()
    oligos[0] = Oligo("A_LEFT", 1, FWD[:-1] + "R", "+", 100, 119)  # ends in A or G
    row = place(aln.rows[0], 1100, FWD[:-1] + "G")
    assert sites(check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900), "A_LEFT") == 1


def test_a_secondary_site_facing_a_primer_makes_an_unintended_product() -> None:
    oligos, aln = clean()
    row = place(aln.rows[0], 450, FWD)                            # a second FWD site before REV's own
    rep = check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900)
    [product] = rep.products
    assert (product.fwd, product.rev) == ("A_LEFT", "A_RIGHT")
    assert product.fwd_secondary and not product.rev_secondary
    assert product.length_min == product.length_max == 700 - 450 and product.sequences == 1


def test_a_product_longer_than_max_product_is_not_reported() -> None:
    oligos, aln = clean()
    row = place(aln.rows[0], 450, FWD)
    assert check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=200).products == []


def test_primers_in_different_pools_do_not_make_products() -> None:
    oligos, aln = clean()
    oligos[1] = Oligo("A_RIGHT", 2, REV, "-", 680, 699)
    row = place(aln.rows[0], 450, FWD)
    assert check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900).products == []


def test_gaps_in_the_row_do_not_shift_the_designed_site() -> None:
    oligos, aln = clean()
    gapped = "--" + aln.rows[0][:50] + "--" + aln.rows[0][50:]
    shifted = [Oligo(o.name, o.pool, o.seq, o.strand, o.first_col + 4, o.last_col + 4) for o in oligos]
    rep = check_specificity(shifted, Alignment(["s1"], [gapped]), P, max_product=900)
    assert all(p.off_target_sites == 0 for p in rep.primers)


def test_background_sequences_only_count_as_secondary_sites() -> None:
    oligos, aln = clean()
    host = place(genome(800, seed=3, avoid=(FWD, REV)), 300, FWD)
    rep = check_specificity(oligos, aln, P, max_product=900, background=[("host", host)])
    by_name = {p.name: p for p in rep.primers}
    assert by_name["A_LEFT"].background_sites == 1 and by_name["A_LEFT"].off_target_sites == 0
    assert by_name["A_RIGHT"].background_sites == 0


def test_every_sequence_is_counted_even_when_identical() -> None:
    oligos, aln = clean()
    row = place(aln.rows[0], 1100, FWD)
    rep = check_specificity(oligos, Alignment(["s1", "s2", "s3"], [row] * 3), P, max_product=900)
    left = {p.name: p for p in rep.primers}["A_LEFT"]
    assert left.off_target_sites == 3 and left.sequences == 3


def test_n_in_a_sequence_never_forms_a_site() -> None:
    oligos, aln = clean()
    row = place(aln.rows[0], 1100, "N" * 30)
    assert all(p.off_target_sites == 0
               for p in check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900).primers)


@pytest.mark.parametrize("bad", ["", "ACGT" * 3])
def test_empty_or_too_short_oligos_are_refused(bad: str) -> None:
    aln = Alignment(["s1"], ["ACGT" * 100])
    with pytest.raises(ValueError):
        check_specificity([Oligo("x", 1, bad, "+", 0, 11)], aln, P, max_product=900)


@pytest.mark.parametrize("tail, path", [("RRYYY", "regex"), ("NNNNN", "every position")])
def test_a_very_degenerate_3_prime_end_still_finds_its_sites(tail: str, path: str) -> None:
    """More than 16 seed variants are looked up by regex, more than 256 at every position."""
    oligos, aln = clean()
    degenerate = FWD[:-5] + tail
    oligos[0] = Oligo("A_LEFT", 1, degenerate, "+", 100, 119)
    row = place(aln.rows[0], 1100, FWD)
    assert sites(check_specificity(oligos, Alignment(["s1"], [row]), P, max_product=900), "A_LEFT") == 1
    assert sites(check_specificity(oligos, aln, P, max_product=900), "A_LEFT") == 0
