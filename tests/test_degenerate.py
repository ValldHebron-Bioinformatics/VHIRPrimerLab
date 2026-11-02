from __future__ import annotations

from viraprime.degenerate import degenerate_pair, fold, tm_range
from viraprime.models import Alignment, ScoringParams

STRICT = ScoringParams(max_mismatches=0, window=5)
FWD = "ACGTTGCAAGCTTCGGATCC"
REV = "TTGACCGTAGGCATCGATGC"
FCOLS = list(range(20))
RCOLS = list(range(40, 60))


def row(fwd: str = FWD, rev: str = REV) -> str:
    return fwd + "A" * 20 + rev


def sub(seq: str, pos: int, base: str) -> str:
    return seq[:pos] + base + seq[pos + 1:]


def test_snp_becomes_iupac() -> None:
    aln = Alignment(["a", "b", "c", "d"], [row(), row(), row(sub(FWD, 3, "C")), row(sub(FWD, 3, "C"))])
    # FWD[3] is T; adding C gives Y
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f[3] == "Y" and f[:3] + f[4:] == FWD[:3] + FWD[4:]
    assert r == REV
    assert cov.coverage_pct == 100.0


def test_reverse_site_degenerated() -> None:
    aln = Alignment(["a", "b"], [row(), row(rev=sub(REV, 10, "A"))])  # REV[10] is G -> R
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert r[10] == "R" and f == FWD and cov.coverage_pct == 100.0


def test_fold_cap_respected() -> None:
    variants = [row()] + [row(sub(FWD, p, "A" if FWD[p] != "A" else "C")) for p in (1, 4, 7, 10)]
    aln = Alignment([f"s{i}" for i in range(5)], variants)
    f, _, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=4)
    assert fold(f) <= 4
    assert cov.coverage_pct < 100.0
    assert len(cov.covered) == 3


def test_three_prime_terminal_never_degenerate() -> None:
    aln = Alignment(["a", "b"], [row(), row(sub(FWD, 19, "A"))])
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f == FWD and r == REV
    assert cov.covered == ["a"]


def test_reverse_three_prime_terminal_never_degenerate() -> None:
    aln = Alignment(["a", "b"], [row(), row(rev=sub(REV, 0, "G"))])
    _, r, _ = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert r == REV


def test_fold() -> None:
    assert fold("ARYN") == 16
    assert fold("ACGT") == 1


def test_tm_range_orders() -> None:
    lo, hi = tm_range("ACGTRACGTACGTACGTAC")
    assert lo < hi
    lo2, hi2 = tm_range("ACGTAACGTACGTACGTAC")
    assert lo2 == hi2


def test_co_occurring_mismatches_rescued() -> None:
    # seq b differs at one forward and one reverse position: no single change rescues it
    b = row(sub(FWD, 3, "C"), sub(REV, 10, "A"))
    aln = Alignment(["a", "b"], [row(), b])
    f, r, cov = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f[3] == "Y" and r[10] == "R"
    assert cov.coverage_pct == 100.0


def test_no_useless_degeneracy_when_rescue_impossible() -> None:
    # seq b also has a 3'-terminal mismatch, so it can never be rescued
    b = row(sub(sub(FWD, 3, "C"), 19, "A"))
    aln = Alignment(["a", "b"], [row(), b])
    f, r, _ = degenerate_pair(FWD, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f == FWD and r == REV


def test_degenerate_bases_counts_non_acgt() -> None:
    from viraprime.degenerate import degenerate_bases
    assert degenerate_bases("ACGTACGT") == 0
    assert degenerate_bases("ACRTYCGN") == 3


def test_gc_percent_is_mean_over_variants() -> None:
    from viraprime.degenerate import gc_percent
    assert gc_percent("GGCC") == 100.0
    assert gc_percent("AATT") == 0.0
    assert gc_percent("ASAT") == 25.0  # S is always G or C
    assert gc_percent("ARAT") == 12.5  # R is A or G: half a G/C on average


def test_degenerate_pair_keeps_gc_within_limits() -> None:
    from viraprime.degenerate import gc_percent
    fwd = "AATTGCAGCTGCATGCATAT"  # 40 % GC; a degenerate R at position 10 would drop it to 37.5 %
    aln = Alignment(["a", "b"], [row(fwd=fwd), row(fwd=sub(fwd, 10, "A"))])
    f, _, _ = degenerate_pair(fwd, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8)
    assert f != fwd
    f, _, cov = degenerate_pair(fwd, FCOLS, REV, RCOLS, aln, STRICT, max_fold=8, gc_limits=(40, 60))
    assert f == fwd and gc_percent(f) == 40.0 and cov.not_covered == {"b": cov.not_covered["b"]}
