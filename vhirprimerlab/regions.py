"""Parsing of target regions from --region specs or a --regions-file."""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from vhirprimerlab.models import PipelineError, Region

NAME_RE = re.compile(r"[A-Za-z0-9_.-]+")
SPEC_RE = re.compile(r"^(?:(?P<name>[^:]*):)?(?P<start>\d+)-(?P<end>\d+)$")


def _make_region(name: str, start: int, end: int, where: str) -> Region:
    if not NAME_RE.fullmatch(name):
        raise PipelineError(f"{where}: invalid region name '{name}' (use letters, digits, _ . -)")
    if start < 1:
        raise PipelineError(f"{where}: region start must be >= 1 (1-based coordinates)")
    if start >= end:
        raise PipelineError(f"{where}: region start ({start}) must be smaller than end ({end})")
    return Region(name, start, end)


def parse_region_spec(spec: str, index: int) -> Region:
    """Parse 'START-END' or 'NAME:START-END'; unnamed regions become region<index>."""
    m = SPEC_RE.match(spec.strip())
    if not m:
        raise PipelineError(f"Invalid --region '{spec}': expected START-END or NAME:START-END")
    name = m.group("name")
    if name is None:
        name = f"region{index}"
    return _make_region(name, int(m.group("start")), int(m.group("end")), f"--region '{spec}'")


def load_regions_file(path: Path) -> list[Region]:
    """Read a tab-separated name/start/end file (optional header, # comments)."""
    path = Path(path)
    if not path.is_file():
        raise PipelineError(f"Regions file not found: {path}")
    regions: list[Region] = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.rstrip("\r").split("\t")
        where = f"{path.name} line {lineno}"
        if len(fields) != 3:
            raise PipelineError(f"{where}: expected 3 tab-separated columns (name, start, end)")
        name, start, end = (f.strip() for f in fields)
        if not regions and not start.isdigit():
            continue  # header line
        if not (start.isdigit() and end.isdigit()):
            raise PipelineError(f"{where}: start and end must be integers")
        regions.append(_make_region(name, int(start), int(end), where))
    if not regions:
        raise PipelineError(f"Regions file contains no regions: {path}")
    return regions


def resolve_regions(specs: list[str] | None, file: Path | None, ref_length: int) -> list[Region]:
    """Combine CLI region sources into a list of uniquely named regions."""
    if specs and file:
        raise PipelineError("Use either --region or --regions-file, not both")
    if file:
        regions = load_regions_file(file)
    elif specs:
        regions = [parse_region_spec(s, i) for i, s in enumerate(specs, start=1)]
    else:
        return [Region("full", 1, ref_length)]
    dupes = sorted(name for name, n in Counter(r.name for r in regions).items() if n > 1)
    if dupes:
        raise PipelineError(f"Duplicate region names: {', '.join(dupes)}")
    return regions
