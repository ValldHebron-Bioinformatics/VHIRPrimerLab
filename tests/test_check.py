from __future__ import annotations

from pathlib import Path

import pytest

from vhirprimerlab.check import load_assays
from vhirprimerlab.models import Assay, PipelineError

F = "ACGTACGTACGTACGTAC"
R = "TTGGCCAATTGGCCAATT"


def _write(tmp_path: Path, data: bytes | str) -> Path:
    p = tmp_path / "assays.tsv"
    if isinstance(data, str):
        data = data.encode()
    p.write_bytes(data)
    return p


# ---------- load_assays ----------

def test_load_assays_header_comments_blank(tmp_path: Path) -> None:
    p = _write(tmp_path, f"name\tforward\treverse\n# c\n\nA1\t{F.lower()}\t{R}\n")
    assert load_assays(p, 5) == [Assay("A1", F, R)]


def test_load_assays_u_to_t(tmp_path: Path) -> None:
    p = _write(tmp_path, f"A1\tACGUACGUACGUACGUAC\t{R}\n")
    assert load_assays(p, 5)[0].forward == F


def test_load_assays_excel_export(tmp_path: Path) -> None:
    data = b"\xef\xbb\xbfname\tforward\treverse\r\nA1\t" + F.encode() + b"\t" + R.encode() + b"\t\t\r\n"
    assert load_assays(_write(tmp_path, data), 5) == [Assay("A1", F, R)]


def test_load_assays_strips_field_spaces(tmp_path: Path) -> None:
    p = _write(tmp_path, f" A1 \t {F} \t{R}\n")
    assert load_assays(p, 5) == [Assay("A1", F, R)]


def test_load_assays_missing_file(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="not found"):
        load_assays(tmp_path / "nope.tsv", 5)


def test_load_assays_header_only(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="no primer pairs"):
        load_assays(_write(tmp_path, "name\tforward\treverse\n"), 5)


def test_load_assays_wrong_field_count(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="line 1: expected 3 tab-separated columns"):
        load_assays(_write(tmp_path, f"A1\t{F}\n"), 5)


def test_load_assays_invalid_name(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="invalid name"):
        load_assays(_write(tmp_path, f"a b\t{F}\t{R}\n"), 5)


def test_load_assays_duplicate_names(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="duplicate.*A1"):
        load_assays(_write(tmp_path, f"A1\t{F}\t{R}\nA1\t{F}\t{R}\n"), 5)


def test_load_assays_invalid_character(tmp_path: Path) -> None:
    text = f"name\tforward\treverse\nA1\tACGTXCGTACGTACGTAC\t{R}\n"
    with pytest.raises(PipelineError, match=r"forward.*position 5.*'X'"):
        load_assays(_write(tmp_path, text), 5)


def test_load_assays_internal_space(tmp_path: Path) -> None:
    text = f"name\tforward\treverse\nA1\tACGT ACGTACGTACGTAC\t{R}\n"
    with pytest.raises(PipelineError, match="position 5"):
        load_assays(_write(tmp_path, text), 5)


def test_load_assays_too_short_for_window(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match="shorter than --three-prime-window"):
        load_assays(_write(tmp_path, f"A1\tACGT\t{R}\n"), 5)


def test_load_assays_typo_on_first_line_is_not_a_header(tmp_path: Path) -> None:
    with pytest.raises(PipelineError, match=r"line 1: forward.*position 5.*'X'"):
        load_assays(_write(tmp_path, f"A1\tACGTXCGTACGTACGTAC\t{R}\n"), 5)


# ---------- locate / check_assay ----------

import logging
import random
import time
from itertools import product

import primer3
from Bio.Seq import reverse_complement

from vhirprimerlab.check import TM_FOLD_LIMIT, check_assay, locate
from vhirprimerlab.conservation import build_consensus, score_columns
from vhirprimerlab.degenerate import fold, tm_range
from vhirprimerlab.models import Alignment, Consensus, Location, ScoringParams
from vhirprimerlab.scoring import IUPAC, score_pair

PARAMS = ScoringParams(max_mismatches=2, window=5)


def _family(mutate: dict[int, dict[int, str]] | None = None) -> tuple[Alignment, Consensus]:
    """6 rows over a seeded 140 bp base, SNPs at 60 and 120, plus one mostly-gap column at 70."""
    rng = random.Random(3)
    base = [rng.choice("ACGT") for _ in range(140)]
    rows = []
    for k in range(6):
        seq = list(base)
        if k == 1:
            seq[60] = "A" if seq[60] != "A" else "C"
        if k == 2:
            seq[120] = "G" if seq[120] != "G" else "T"
        for pos, b in (mutate or {}).get(k, {}).items():
            seq[pos] = b
        seq.insert(70, "T" if k == 0 else "-")
        rows.append("".join(seq))
    aln = Alignment([f"s{k}" for k in range(6)], rows)
    return aln, build_consensus(score_columns(aln))


def _gc(seq: str) -> float:
    return 100 * sum(b in "GC" for b in seq) / len(seq)


def test_locate_exact() -> None:
    _, cons = _family()
    assert locate(cons.seq[30:50], cons.seq) == Location(30, 0)


def test_locate_best_with_mismatches() -> None:
    _, cons = _family()
    site = list(cons.seq[30:50])
    site[5] = "A" if site[5] != "A" else "C"
    assert locate("".join(site), cons.seq) == Location(30, 1)


def test_locate_degenerate() -> None:
    _, cons = _family()
    site = list(cons.seq[30:50])
    i = next(i for i, b in enumerate(site) if b in "AG")
    site[i] = "R"
    assert locate("".join(site), cons.seq) == Location(30, 0)


def test_locate_ambiguous() -> None:
    assert locate("ACGTAC", "TTACGTACGGACGTACTT") is None


def test_locate_too_long() -> None:
    assert locate("A" * 30, "ACGT") is None


def test_check_assay_matches_score_pair() -> None:
    aln, cons = _family()
    fwd, rev_site = cons.seq[10:30], cons.seq[80:100]
    res = check_assay(Assay("P", fwd, reverse_complement(rev_site)), aln, cons, PARAMS, 95)
    assert res.located is True
    assert res.coverage == score_pair(
        fwd, cons.cons_to_col[10:30], rev_site, cons.cons_to_col[80:100], aln, PARAMS
    )
    assert res.fwd == Location(10, 0) and res.rev == Location(80, 0)
    assert res.fwd_tm == pytest.approx(primer3.calc_tm(fwd)) and res.fwd_gc == pytest.approx(_gc(fwd))
    assert res.below_threshold is False


def test_check_assay_wrong_order() -> None:
    aln, cons = _family()
    res = check_assay(Assay("P", cons.seq[80:100], reverse_complement(cons.seq[10:30])), aln, cons, PARAMS, 95)
    assert res.located is False
    assert res.note == "reverse site is not downstream of forward site"
    assert res.coverage.coverage_pct == 0.0
    assert res.coverage.not_covered == {sid: "not located" for sid in aln.ids}
    assert res.below_threshold is True


def test_check_assay_overlap() -> None:
    aln, cons = _family()
    res = check_assay(Assay("P", cons.seq[10:30], reverse_complement(cons.seq[29:49])), aln, cons, PARAMS, 95)
    assert res.located is False
    assert res.note == "reverse site is not downstream of forward site"


def test_check_assay_degenerate_tm() -> None:
    aln, cons = _family()
    fwd = list(cons.seq[10:30])
    i = next(i for i, b in enumerate(fwd) if b in "AG")
    fwd[i] = "R"
    fwd = "".join(fwd)
    res = check_assay(Assay("P", fwd, reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 95)
    assert res.degenerate is True
    assert res.fwd_tm_range == tm_range(fwd)
    assert res.fwd_tm == pytest.approx(sum(res.fwd_tm_range) / 2)
    expansions = ["".join(v) for v in product(*(sorted(IUPAC[b]) for b in fwd))]
    assert res.fwd_gc == pytest.approx(sum(_gc(e) for e in expansions) / len(expansions))


def test_check_assay_high_fold_skips_tm() -> None:
    aln, cons = _family()
    fwd = "N" * 12 + cons.seq[22:30]
    t0 = time.monotonic()
    res = check_assay(Assay("P", fwd, reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 95)
    assert time.monotonic() - t0 < 5
    assert TM_FOLD_LIMIT == 256
    assert res.fwd_tm is None and res.fwd_tm_range is None
    assert "Tm not computed: fwd fold 16777216 > 256" in res.note


def test_check_assay_wrong_orientation_warns(caplog: pytest.LogCaptureFixture) -> None:
    aln, cons = _family()
    with caplog.at_level(logging.WARNING, logger="vhirprimerlab"):
        check_assay(Assay("P", cons.seq[10:30], cons.seq[80:100]), aln, cons, PARAMS, 95)
    assert any("wrong orientation" in r.getMessage() for r in caplog.records)


def test_check_assay_below_threshold() -> None:
    _, cons0 = _family()
    three_prime = cons0.seq[29]
    aln, cons = _family({3: {29: "A" if three_prime != "A" else "C"}})
    res = check_assay(Assay("P", cons.seq[10:30], reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 100)
    assert res.below_threshold is True
    assert "s3" in res.coverage.not_covered


# ---------- outputs ----------

import csv

from vhirprimerlab.coords import CoordMap
from vhirprimerlab.io import CHECK_TSV_COLUMNS, write_check_report, write_check_tsv
from vhirprimerlab.models import CheckReport, LengthFilterResult, Record
from vhirprimerlab.ranking import count_reasons


def _results(aln: Alignment, cons: Consensus) -> list:
    ok = check_assay(Assay("A1", cons.seq[10:30], reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 95)
    bad = check_assay(Assay("A2", cons.seq[80:100], reverse_complement(cons.seq[10:30])), aln, cons, PARAMS, 95)
    return [ok, bad]


def _read(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_count_reasons() -> None:
    nc = {"s1": "fwd 3' mismatch; rev too many mismatches", "s2": "fwd 3' mismatch"}
    assert count_reasons(nc) == [("fwd 3' mismatch", 2), ("rev too many mismatches", 1)]


def test_write_check_tsv_columns(tmp_path: Path) -> None:
    aln, cons = _family()
    cmap = CoordMap(aln, "s0", cons)
    ok, bad = _results(aln, cons)
    out = tmp_path / "check.tsv"
    write_check_tsv([ok, bad], cmap, out)
    assert out.read_text().splitlines()[0].split("\t") == CHECK_TSV_COLUMNS
    r_ok, r_bad = _read(out)
    start = cmap.col_to_ref(ok.fwd_cols[0], "start")
    end = cmap.col_to_ref(ok.rev_cols[-1], "end")
    assert (r_ok["amplicon_start"], r_ok["amplicon_end"], r_ok["amplicon_len"]) == (
        str(start), str(end), str(end - start + 1))
    assert r_ok["located"] == "True" and r_ok["coverage_pct"] == "100.0" and r_ok["n_total"] == "6"
    assert r_bad["located"] == "False" and r_bad["note"] == NOT_DOWNSTREAM_TEXT
    assert r_bad["amplicon_start"] == r_bad["amplicon_len"] == r_bad["forward_site_mismatches"] == ""
    assert r_bad["not_covered"].startswith("s0(not located)")
    for r, res in ((r_ok, ok), (r_bad, bad)):
        for side, seq in (("forward", res.assay.forward), ("reverse", res.assay.reverse)):
            assert r[f"{side}_degenerate_bases"] == str(sum(b not in "ACGT" for b in seq))
            assert r[f"{side}_variants"] == str(fold(seq))


NOT_DOWNSTREAM_TEXT = "reverse site is not downstream of forward site"


def test_write_check_tsv_reference_gap(tmp_path: Path) -> None:
    aln, cons = _family({0: {10: "-"}})
    cmap = CoordMap(aln, "s0", cons)
    ok = _results(aln, cons)[0]
    write_check_tsv([ok], cmap, tmp_path / "c.tsv")
    assert _read(tmp_path / "c.tsv")[0]["amplicon_start"] == "11"


def _report(aln: Alignment, results: list, min_coverage: float = 95.0) -> CheckReport:
    recs = [Record(i, r.replace("-", "")) for i, r in zip(aln.ids, aln.rows)]
    return CheckReport(
        command="vhirprimerlab check ...", params={"min_coverage": 95.0}, n_input=6,
        median_len=140.0, min_len=140, max_len=141,
        length=LengthFilterResult(recs, [], 140.0, "drop"), reference_id="s0",
        aln_length=aln.length, cons_length=140, n_conserved=130, mean_conservation=0.99,
        min_coverage=min_coverage, results=results, log_lines=["INFO: done"],
    )


def test_write_check_report_sections(tmp_path: Path) -> None:
    aln, cons = _family({3: {29: "A" if _family()[1].seq[29] != "A" else "C"}})
    cmap = CoordMap(aln, "s0", cons)
    low = check_assay(Assay("A1", cons.seq[10:30], reverse_complement(cons.seq[80:100])), aln, cons, PARAMS, 100)
    bad = check_assay(Assay("A2", cons.seq[80:100], reverse_complement(cons.seq[10:30])), aln, cons, PARAMS, 100)
    out = tmp_path / "report.txt"
    write_check_report(_report(aln, [low, bad], 100.0), cmap, out)
    text = out.read_text()
    assert "VHIRPrimerLab check report" in text
    assert "Pair A1" in text and "Pair A2" in text
    assert "Consensus positions: forward 11-30, reverse 81-100" in text
    assert f"Forward: {low.assay.forward}  20 nt, Tm {low.fwd_tm:.1f} °C, GC {low.fwd_gc:.1f} %" in text
    assert f"Reverse: {low.assay.reverse}  20 nt, Tm {low.rev_tm:.1f} °C" in text
    assert "Why it misses sequences" in text and "fwd 3' mismatch: 1" in text
    assert "not located" in text and NOT_DOWNSTREAM_TEXT in text
    assert "Pairs below --min-coverage (100%): A1, A2" in text
    assert text.rstrip().endswith("INFO: done")


def test_write_check_report_all_pass(tmp_path: Path) -> None:
    aln, cons = _family()
    cmap = CoordMap(aln, "s0", cons)
    ok = _results(aln, cons)[0]
    write_check_report(_report(aln, [ok]), cmap, tmp_path / "r.txt")
    assert "Every pair meets --min-coverage." in (tmp_path / "r.txt").read_text()


# ---------- final-review fixes ----------

from tests.synthetic import make_aligned_family


def _aligned_family() -> tuple[Alignment, Consensus]:
    rows = make_aligned_family()
    aln = Alignment([i for i, _ in rows], [s for _, s in rows])
    return aln, build_consensus(score_columns(aln))


def test_amplicon_outside_reference_does_not_crash(tmp_path: Path) -> None:
    aln, cons = _family({0: {p: "-" for p in range(100, 140)}})  # reference s0 ends at base 100
    cmap = CoordMap(aln, "s0", cons)
    res = check_assay(Assay("A1", cons.seq[102:120], reverse_complement(cons.seq[121:139])), aln, cons, PARAMS, 95)
    assert res.located
    write_check_tsv([res], cmap, tmp_path / "c.tsv")
    row = _read(tmp_path / "c.tsv")[0]
    assert row["amplicon_start"] == row["amplicon_len"] == ""
    assert "amplicon not on reference s0" in row["note"]
    write_check_report(_report(aln, [res]), cmap, tmp_path / "r.txt")
    assert "amplicon not on reference s0" in (tmp_path / "r.txt").read_text()


def test_reverse_pasted_as_its_complement_gets_orientation_note() -> None:
    aln, cons = _aligned_family()
    fwd, rev_site = cons.seq[100:121], cons.seq[300:321]
    res = check_assay(Assay("P", fwd, rev_site), aln, cons, PARAMS, 95)  # should be reverse_complement(rev_site)
    assert res.below_threshold is True
    assert "rev primer may be written in the wrong orientation" in res.note
    assert "ambiguous" not in res.note


def test_no_site_within_mismatches_is_not_called_ambiguous() -> None:
    aln, cons = _aligned_family()
    alien = "GATTACAGATTACAGATTACA"
    res = check_assay(Assay("P", cons.seq[100:121], alien), aln, cons, PARAMS, 95)
    assert "ambiguous" not in res.note
    assert "no binding site for rev primer within 2 mismatches" in res.note


def test_not_located_is_below_threshold_even_at_zero_coverage() -> None:
    aln, cons = _family()
    res = check_assay(Assay("P", cons.seq[80:100], reverse_complement(cons.seq[10:30])), aln, cons, PARAMS, 0)
    assert res.located is False and res.below_threshold is True
