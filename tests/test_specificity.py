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


# --- the vectorised scan agrees with the base-by-base one -----------------------------

@pytest.mark.parametrize("seed", range(150))
def test_the_numpy_scan_finds_what_the_plain_scan_finds(seed: int) -> None:
    from vhirprimerlab.specificity import _Coded, _Scanner

    rng = random.Random(seed)
    length = rng.randint(15, 30)
    oligo = "".join(rng.choice("ACGT") for _ in range(length))
    for _ in range(rng.choice([0, 0, 1, 3])):                         # some degenerate positions
        p = rng.randrange(length)
        oligo = oligo[:p] + rng.choice("RYKMSWBDHVN") + oligo[p + 1:]
    window, mismatches = rng.randint(2, 12), rng.randint(0, 3)
    text = "".join(rng.choice("ACGT") for _ in range(rng.choice([5, length, 400, 3000])))
    sites = []
    for _ in range(rng.randint(0, 6)):                                # plant near-copies of the oligo
        concrete = "".join(rng.choice(sorted(IUPAC_BASES[c])) for c in oligo)
        for _ in range(rng.randint(0, 4)):
            q = rng.randrange(length)
            concrete = concrete[:q] + rng.choice("ACGT") + concrete[q + 1:]
        at = rng.randrange(0, max(1, len(text) - length + 1))
        sites.append(at)
        text = text[:at] + concrete + text[at + length:]
    if rng.random() < 0.5:                                            # N and gaps-turned-N in the text
        for _ in range(rng.randint(1, 20)):
            q = rng.randrange(len(text))
            text = text[:q] + rng.choice("NRY") + text[q + 1:]
    scanner = _Scanner(oligo, window, mismatches)
    assert sorted(scanner.scan(_Coded(text))) == sorted(scanner.scan_text(text))


IUPAC_BASES = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "K": "GT", "M": "AC",
               "S": "CG", "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def test_the_numpy_scan_actually_finds_planted_sites() -> None:
    from vhirprimerlab.specificity import _Coded, _Scanner

    text = genome(2000, seed=11, avoid=(FWD, REV))
    text = place(place(text, 300, FWD), 1500, FWD)
    found = _Scanner(FWD, 5, 2).scan(_Coded(text))
    assert (300, 0) in found and (1500, 0) in found


# --- check mode -----------------------------------------------------------------------

def _check_setup(tmp_path):
    """The synthetic genome family, a tiling of it, and its primers as pairs to check."""
    import csv

    from tests.synthetic import make_genome_family, write_fasta
    from vhirprimerlab.cli import main

    rows = make_genome_family()
    aln = write_fasta(tmp_path / "aln.fasta", rows)
    assert main(["--aligned", str(aln), "-n", "5", "-o", str(tmp_path / "design")]) in (0, 2)
    with (tmp_path / "design" / "amplicons.tsv").open() as fh:
        amplicons = list(csv.DictReader(fh, delimiter="\t"))
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("".join(f"a{r['amplicon']}\t{r['forward_seq']}\t{r['reverse_seq']}\n" for r in amplicons))
    return rows, pairs, amplicons


def _run_check(tmp_path, aln, pairs, name="chk", *flags):
    from vhirprimerlab.cli import main

    out = tmp_path / name
    assert main(["check", "--aligned", str(aln), "--primers", str(pairs), "-o", str(out), *flags]) in (0, 2)
    return (out / "report.txt").read_text()


def test_check_mode_scans_for_secondary_sites_and_a_clean_scheme_says_so(tmp_path) -> None:
    rows, pairs, _ = _check_setup(tmp_path)
    text = _run_check(tmp_path, tmp_path / "aln.fasta", pairs)
    assert "Specificity\n-----------" in text
    assert "No secondary binding sites found." in text


def test_check_mode_reports_a_second_site_of_a_users_primer(tmp_path) -> None:
    from tests.synthetic import write_fasta

    rows, pairs, amplicons = _check_setup(tmp_path)
    forward = amplicons[0]["forward_seq"]
    # a near-copy (one mismatch) far from its own site, in every sequence; an exact copy would make
    # check call the primer's site ambiguous and leave the pair unlocated (see the next test)
    near = forward[:8] + ("A" if forward[8] != "A" else "C") + forward[9:]
    far = len(rows[0][1]) - 400
    planted = [(sid, seq[:far] + near + seq[far + len(near):]) for sid, seq in rows]
    aln = write_fasta(tmp_path / "planted.fasta", planted)
    text = _run_check(tmp_path, aln, pairs, "chk2")
    assert f"a1 forward: {len(rows)} secondary site(s) in {len(rows)} sequence(s)" in text
    assert "No secondary binding sites found." not in text


def test_check_mode_reports_the_product_a_second_site_could_make(tmp_path) -> None:
    from tests.synthetic import write_fasta
    from vhirprimerlab.specificity import reverse_complement

    rows, pairs, amplicons = _check_setup(tmp_path)
    forward, reverse = amplicons[0]["forward_seq"], amplicons[0]["reverse_seq"]
    sid, seq = rows[0]
    at = seq.find(reverse_complement(reverse)) - 250        # a second forward site facing the reverse primer
    assert at > 0
    near = forward[:8] + ("A" if forward[8] != "A" else "C") + forward[9:]
    planted = [(s, q[:at] + near + q[at + len(near):]) for s, q in rows]
    aln = write_fasta(tmp_path / "planted.fasta", planted)
    text = _run_check(tmp_path, aln, pairs, "chk3")
    assert "Unintended products" in text and "a1 forward × a1 reverse" in text


def test_check_mode_scans_a_background(tmp_path) -> None:
    from tests.synthetic import write_fasta

    rows, pairs, amplicons = _check_setup(tmp_path)
    host = write_fasta(tmp_path / "host.fasta", [("host", "ACGT" * 50 + amplicons[1]["reverse_seq"] + "TTGA" * 50)])
    text = _run_check(tmp_path, tmp_path / "aln.fasta", pairs, "chk4", "--background", str(host))
    assert "background sequence(s)" in text
    assert "a2 reverse: 1 site(s) in the background sequences" in text


def test_check_mode_says_which_pairs_it_could_not_scan(tmp_path) -> None:
    rows, pairs, _ = _check_setup(tmp_path)
    pairs.write_text(pairs.read_text() + "nowhere\tGATTACAGATTACAGATTACA\tCCCGGGAAATTTCCCGGGAAA\n"
                     "short\tACGTACGTAC\tTTGGCCAATT\n")
    text = _run_check(tmp_path, tmp_path / "aln.fasta", pairs, "chk5")
    spec = text.split("Specificity")[1]
    assert "not scanned: nowhere" in spec and "not scanned: short" in spec


def test_check_mode_without_any_located_pair_still_finishes(tmp_path) -> None:
    rows, _, _ = _check_setup(tmp_path)
    nowhere = tmp_path / "nowhere.tsv"
    nowhere.write_text("x\tGATTACAGATTACAGATTACA\tCCCGGGAAATTTCCCGGGAAA\n")
    text = _run_check(tmp_path, tmp_path / "aln.fasta", nowhere, "chk6")
    assert "not scanned: x" in text


def test_an_exact_second_site_leaves_the_pair_unlocated_and_says_so(tmp_path) -> None:
    from tests.synthetic import write_fasta

    rows, pairs, amplicons = _check_setup(tmp_path)
    forward = amplicons[0]["forward_seq"]
    far = len(rows[0][1]) - 400
    aln = write_fasta(tmp_path / "twice.fasta",
                      [(sid, seq[:far] + forward + seq[far + len(forward):]) for sid, seq in rows])
    text = _run_check(tmp_path, aln, pairs, "chk7")
    assert "not scanned: a1 (ambiguous binding site for fwd primer" in text
