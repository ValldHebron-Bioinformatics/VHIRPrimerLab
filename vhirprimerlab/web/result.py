"""Turn a finished job directory into the result the page shows (server and browser engines)."""

from __future__ import annotations

import csv
from pathlib import Path

OUTPUT_FILES = (
    "primers.tsv", "amplicons.tsv", "pairs.tsv", "scheme.bed", "check.tsv", "report.txt",
    "alignment.fasta",
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def reference(out: Path) -> tuple[str | None, int | None]:
    """Reference ID (from report.txt) and its ungapped length (from alignment.fasta)."""
    report, aln = out / "report.txt", out / "alignment.fasta"
    if not report.is_file() or not aln.is_file():
        return None, None
    ref = next((line.split(": ", 1)[1] for line in report.read_text().splitlines()
                if line.startswith("Reference: ")), None)
    current, length = None, 0
    for line in aln.read_text().splitlines():
        if line.startswith(">"):
            if current == ref:
                break
            current, length = line[1:].split()[0], 0
        elif current == ref:
            length += len(line.strip().replace("-", ""))
    return ref, length if current == ref else None


def build_result(job_dir: Path, argv: list[str], code: int) -> dict:
    """The result object for a run whose log is job_dir/log.txt and outputs job_dir/out."""
    job_dir = Path(job_dir)
    log = job_dir / "log.txt"
    lines = log.read_text().splitlines() if log.is_file() else []
    out = job_dir / "out"
    mode = "check" if argv[:1] == ["check"] else "tile" if "-n" in argv else "regions"
    table = out / {"check": "check.tsv", "tile": "amplicons.tsv", "regions": "pairs.tsv"}[mode]
    errors = [ln[len("ERROR: "):] for ln in lines if ln.startswith("ERROR: ")]
    usage = [ln.split("error: ", 1)[1] for ln in lines if "error: " in ln and ln.startswith("vhirprimerlab")]
    ref, ref_len = reference(out)
    return {
        "exit_code": code,
        "mode": mode,
        "rows": read_rows(table) if table.is_file() and code != 1 else [],
        "reference_id": ref,
        "reference_length": ref_len,
        "warnings": [ln[len("WARNING: "):] for ln in lines if ln.startswith("WARNING: ")],
        "error": (errors or usage or ["the job failed; see the log"])[-1] if code == 1 else None,
    }
