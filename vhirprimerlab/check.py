"""Checking primer pairs the user already has against every aligned sequence."""

from __future__ import annotations

import csv
import itertools
import logging
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree

from Bio.Seq import reverse_complement

from vhirprimerlab.degenerate import fold, gc_percent, tm_range
from vhirprimerlab.models import (
    Alignment,
    Assay,
    CheckResult,
    Consensus,
    Location,
    PairCoverage,
    PipelineError,
    ScoringParams,
)
from vhirprimerlab.regions import NAME_RE
from vhirprimerlab.scoring import IUPAC, score_pair
from vhirprimerlab.structure import dimer_dg, structure_report
from vhirprimerlab.thermo import DEFAULT_THERMO, Thermo

log = logging.getLogger("vhirprimerlab")

TM_FOLD_LIMIT = 256  # above this, enumerating every expansion for a Tm range is too slow
STRUCTURE_FOLD_LIMIT = 64  # above this, checking every pair of variants is too slow
NOT_DOWNSTREAM = "reverse site is not downstream of forward site"
FORWARD_TAGS = {"LEFT", "F", "FWD", "FORWARD"}
SIDE_RE = re.compile(r"(?P<pair>.+?)[_.-](?P<side>LEFT|RIGHT|FWD|REV|FORWARD|REVERSE|F|R)", re.I)


def _clean(field: str) -> str:
    return field.strip().upper().replace("U", "T")


def _looks_like_header(fwd: str) -> bool:
    """A column title like 'forward' is mostly not A/C/G/T; a primer, even with a typo, is."""
    return sum(ch in "ACGT" for ch in fwd) * 2 < len(fwd)


def _check_oligo(seq: str, what: str, where: str, window: int) -> None:
    for pos, ch in enumerate(seq, start=1):
        if ch not in IUPAC:
            raise PipelineError(f"{where}: {what} primer has an invalid character at position {pos}: '{ch}'")
    if len(seq) < window:
        raise PipelineError(
            f"{where}: {what} primer ({len(seq)} bases) is shorter than --three-prime-window ({window})"
        )


def _split_side(name: str, where: str) -> tuple[str, str]:
    """('pair1', 'forward') from an order-sheet name such as pair1_LEFT or pair1_R."""
    m = SIDE_RE.fullmatch(name)
    if not m:
        raise PipelineError(
            f"{where}: primer name '{name}' does not say which side it is; end it in "
            "_LEFT or _RIGHT (or _F/_R, _FWD/_REV, _FORWARD/_REVERSE)"
        )
    return m["pair"], "forward" if m["side"].upper() in FORWARD_TAGS else "reverse"


_XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _column(ref: str) -> int:
    """0-based column of a cell reference such as 'C7'."""
    n = 0
    for ch in ref:
        if not ch.isalpha():
            break
        n = n * 26 + ord(ch.upper()) - 64
    return n - 1


def _xlsx_rows(path: Path) -> list[tuple[str, list[str]]]:
    """(label, cell texts) for each row of the first worksheet, using only the standard library."""
    try:
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            shared: list[str] = []
            if "xl/sharedStrings.xml" in names:
                root = ElementTree.fromstring(z.read("xl/sharedStrings.xml"))
                shared = ["".join(t.text or "" for t in si.iter(f"{_XLSX_NS}t"))
                          for si in root.iter(f"{_XLSX_NS}si")]
            book = ElementTree.fromstring(z.read("xl/workbook.xml"))
            first = next(book.iter(f"{_XLSX_NS}sheet")).get(f"{_REL_NS}id")
            rels = ElementTree.fromstring(z.read("xl/_rels/workbook.xml.rels"))
            target = next(r.get("Target") for r in rels if r.get("Id") == first)
            target = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
            sheet = ElementTree.fromstring(z.read(target))
    except (zipfile.BadZipFile, KeyError, StopIteration, ElementTree.ParseError):
        raise PipelineError(f"{path.name} is not a readable .xlsx file; save it as .xlsx, CSV or TSV") from None
    rows: list[tuple[str, list[str]]] = []
    for i, row in enumerate(sheet.iter(f"{_XLSX_NS}row"), start=1):
        cells: dict[int, str] = {}
        for c in row.iter(f"{_XLSX_NS}c"):
            kind, v = c.get("t"), c.find(f"{_XLSX_NS}v")
            if kind == "s" and v is not None:
                text = shared[int(v.text or 0)]
            elif kind == "inlineStr":
                text = "".join(t.text or "" for t in c.iter(f"{_XLSX_NS}t"))
            else:
                text = v.text or "" if v is not None else ""
            cells[_column(c.get("r", "")) if c.get("r") else len(cells)] = text
        width = max(cells, default=-1) + 1
        rows.append((f"row {row.get('r', i)}", [cells.get(j, "") for j in range(width)]))
    return rows


def _text_rows(path: Path) -> list[tuple[str, list[str]]]:
    """(label, fields) per line of a tab-, semicolon- or comma-separated file."""
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")  # Excel's "CSV" on Windows
    lines = text.splitlines()
    first = next((ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")), "")
    delimiter = next((d for d in "\t;," if d in first), "\t")
    return [(f"line {i}", fields) for i, fields in enumerate(csv.reader(lines, delimiter=delimiter), start=1)]


def _rows(path: Path) -> list[tuple[str, list[str]]]:
    if path.suffix.lower() == ".xls":
        raise PipelineError(f"{path.name}: old .xls files are not supported; save it as .xlsx, CSV or TSV")
    if path.suffix.lower() == ".xlsx" or path.read_bytes()[:2] == b"PK":
        return _xlsx_rows(path)
    return _text_rows(path)


def load_assays(path: Path, window: int) -> list[Assay]:
    """Read primer pairs from a TSV, CSV (comma or semicolon) or .xlsx file.

    Either one pair per row (name, forward, reverse) or one primer per row (name,
    sequence) as in the primers.tsv order sheet, where names ending in _LEFT/_RIGHT
    (or _F/_R, ...) are joined into pairs. Optional header row, # comments."""
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"Primer file not found: {path}")
    pairs: dict[str, dict[str, str]] = {}  # pair name -> {"forward": seq, "reverse": seq}
    dups: list[str] = []
    seen_data = False
    for label, row in _rows(path):
        fields = [f.strip() for f in row]
        if not any(fields) or fields[0].startswith("#"):
            continue
        while fields and not fields[-1]:
            fields.pop()
        where = f"{path.name} {label}"
        if len(fields) not in (2, 3):
            raise PipelineError(
                f"{where}: expected 3 columns (name, forward, reverse) or 2 (name, sequence)"
            )
        if not seen_data:
            seen_data = True
            if _looks_like_header(_clean(fields[1])):
                continue  # header line
        if len(fields) == 3:
            name, halves = fields[0], {"forward": _clean(fields[1]), "reverse": _clean(fields[2])}
        else:
            name, side = _split_side(fields[0], where)
            halves = {side: _clean(fields[1])}
        if not NAME_RE.fullmatch(name):
            raise PipelineError(f"{where}: invalid name '{name}' (use letters, digits, _ . -)")
        for side, seq in halves.items():
            _check_oligo(seq, side, where, window)
        have = pairs.setdefault(name, {})
        if have.keys() & halves.keys():
            if len(halves) == 1 and len(have) == 1:
                raise PipelineError(f"{where}: pair {name} has two {next(iter(halves))} primers")
            dups.append(name)
        have.update(halves)
    if dups:
        raise PipelineError(f"{path.name}: duplicate pair names: {', '.join(sorted(set(dups)))}")
    if not pairs:
        raise PipelineError(f"Primer file contains no primer pairs: {path}")
    for name, have in pairs.items():
        missing = {"forward", "reverse"} - have.keys()
        if missing:
            raise PipelineError(f"{path.name}: pair {name} has no {missing.pop()} primer")
    return [Assay(name, have["forward"], have["reverse"]) for name, have in pairs.items()]


def _best_sites(site: str, consensus_seq: str) -> tuple[int, list[int]]:
    """Fewest mismatches of site against the consensus, and every start that reaches it."""
    n = len(site)
    allowed = [IUPAC[b] for b in site]
    best, best_starts = n + 1, []
    for start in range(len(consensus_seq) - n + 1):
        window = consensus_seq[start:start + n]
        mm = sum(c not in ok for c, ok in zip(window, allowed))
        if mm < best:
            best, best_starts = mm, [start]
        elif mm == best:
            best_starts.append(start)
    return best, best_starts


def locate(site: str, consensus_seq: str) -> Location | None:
    """Start on the consensus with the fewest mismatches; None if tied or the site is too long."""
    if len(site) > len(consensus_seq):
        return None
    best, starts = _best_sites(site, consensus_seq)
    return Location(starts[0], best) if len(starts) == 1 else None


def _tm(primer: str, label: str, notes: list[str],
        thermo: Thermo) -> tuple[float | None, tuple[float, float] | None]:
    n = fold(primer)
    if n == 1:
        return thermo.tm(primer), None
    if n > TM_FOLD_LIMIT:
        notes.append(f"Tm not computed: {label} fold {n} > {TM_FOLD_LIMIT}")
        return None, None
    lo, hi = tm_range(primer, thermo)
    return (lo + hi) / 2, (lo, hi)


def _structure(primer: str, thermo: Thermo) -> list[str]:
    n = fold(primer)
    if n > STRUCTURE_FOLD_LIMIT:
        return [f"hairpins and self-dimers not checked: {n} variants > {STRUCTURE_FOLD_LIMIT}"]
    return structure_report(primer, thermo)


def _dimer_dg(a: str, b: str, thermo: Thermo) -> float | None:
    if fold(a) * fold(b) > STRUCTURE_FOLD_LIMIT ** 2:
        return None
    return dimer_dg(a, b, thermo)


def cross_dimers(assays: list[Assay], thermo: Thermo = DEFAULT_THERMO) -> list[tuple[str, str, float]]:
    """Primers of different pairs that stick to each other (ΔG below thermo.max_dimer_dg, cal/mol).

    This matters when the pairs share a tube; forward x reverse of one pair is checked
    by check_assay."""
    primers = [(f"{a.name} {side}", a.name, seq)
               for a in assays for side, seq in (("forward", a.forward), ("reverse", a.reverse))]
    out = []
    for (na, pa, a), (nb, pb, b) in itertools.combinations(primers, 2):
        if pa == pb:
            continue
        dg = _dimer_dg(a, b, thermo)
        if dg is not None and dg < thermo.max_dimer_dg:
            out.append((na, nb, dg))
    return sorted(out, key=lambda d: d[2])


def check_assay(
    assay: Assay, aln: Alignment, cons: Consensus, params: ScoringParams, min_coverage: float,
    thermo: Thermo = DEFAULT_THERMO,
) -> CheckResult:
    """Locate both primers on the consensus and score the pair on every sequence."""
    notes: list[str] = []
    fwd_site, rev_site = assay.forward, reverse_complement(assay.reverse)
    locs: dict[str, Location | None] = {}
    for label, site in (("fwd", fwd_site), ("rev", rev_site)):
        if len(site) > len(cons.seq):
            notes.append(f"{label} primer is longer than the consensus")
            locs[label] = None
            continue
        best, starts = _best_sites(site, cons.seq)
        loc = Location(starts[0], best) if len(starts) == 1 else None
        if best > params.max_mismatches:
            # a unique poor site is still scored (it usually covers nothing); ties are not located
            found = [f"no binding site for {label} primer within {params.max_mismatches} "
                     f"mismatches (best {best})"]
            flipped, _ = _best_sites(reverse_complement(site), cons.seq)
            if flipped <= params.max_mismatches:
                found.append(f"{label} primer may be written in the wrong orientation")
            notes += found
            log.warning("Pair %s: %s", assay.name, "; ".join(found))
        elif loc is None:
            notes.append(f"ambiguous binding site for {label} primer")
        locs[label] = loc
    fwd, rev = locs["fwd"], locs["rev"]
    located = fwd is not None and rev is not None
    if located and rev.start <= fwd.start + len(fwd_site) - 1:
        notes.append(NOT_DOWNSTREAM)
        located = False

    fwd_cols = rev_cols = []
    if located:
        fwd_cols = cons.cons_to_col[fwd.start:fwd.start + len(fwd_site)]
        rev_cols = cons.cons_to_col[rev.start:rev.start + len(rev_site)]
        coverage = score_pair(fwd_site, fwd_cols, rev_site, rev_cols, aln, params)
    else:
        coverage = PairCoverage([], {sid: "not located" for sid in aln.ids}, 0.0)

    fwd_tm, fwd_range = _tm(assay.forward, "fwd", notes, thermo)
    rev_tm, rev_range = _tm(assay.reverse, "rev", notes, thermo)
    return CheckResult(
        assay=assay, located=located, note="; ".join(notes), fwd=fwd, rev=rev,
        fwd_cols=list(fwd_cols), rev_cols=list(rev_cols), coverage=coverage,
        fwd_tm=fwd_tm, rev_tm=rev_tm, fwd_tm_range=fwd_range, rev_tm_range=rev_range,
        fwd_gc=gc_percent(assay.forward), rev_gc=gc_percent(assay.reverse),
        degenerate=fold(assay.forward) > 1 or fold(assay.reverse) > 1,
        below_threshold=not located or coverage.coverage_pct < min_coverage,
        fwd_structure=_structure(assay.forward, thermo), rev_structure=_structure(assay.reverse, thermo),
        pair_dimer_dg=_dimer_dg(assay.forward, assay.reverse, thermo),
    )


MAX_CROSS_WARNINGS = 10


def structure_warnings(results: list[CheckResult], dimers: list[tuple[str, str, float]],
                       thermo: Thermo = DEFAULT_THERMO) -> list[str]:
    """One line per primer that folds or sticks to itself, per pair whose primers stick to
    each other, and per primer pair across pairs (only the strongest few)."""
    out: list[str] = []
    for res in results:
        for side, problems in (("forward", res.fwd_structure), ("reverse", res.rev_structure)):
            if problems:
                out.append(f"Primer structure: {res.assay.name} {side}: {'; '.join(problems)}")
        if res.pair_dimer_dg is not None and res.pair_dimer_dg < thermo.max_dimer_dg:
            out.append(f"Primer dimer: {res.assay.name} forward × reverse, "
                       f"ΔG {res.pair_dimer_dg / 1000:.1f} kcal/mol")
    for a, b, dg in dimers[:MAX_CROSS_WARNINGS]:
        out.append(f"Primer dimer if mixed in one tube: {a} × {b}, ΔG {dg / 1000:.1f} kcal/mol")
    if len(dimers) > MAX_CROSS_WARNINGS:
        out.append(f"Primer dimer if mixed in one tube: {len(dimers) - MAX_CROSS_WARNINGS} more "
                   "in report.txt")
    return out
