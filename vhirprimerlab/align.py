"""Multiple sequence alignment via MAFFT or MUSCLE (subprocess)."""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Literal

from vhirprimerlab.io import alignment_from_records, load_fasta
from vhirprimerlab.models import Alignment, PipelineError, Record

log = logging.getLogger(__name__)

Aligner = Literal["mafft", "muscle"]
AlignFn = Callable[[list[Record]], Alignment]

# aligners provided by the host program (e.g. a WebAssembly MAFFT in the browser)
_REGISTERED: dict[str, AlignFn] = {}


def register_aligner(name: str, fn: AlignFn | None) -> None:
    """Use fn instead of the external program `name`; None removes the registration."""
    if fn is None:
        _REGISTERED.pop(name, None)
    else:
        _REGISTERED[name] = fn


def check_binary(tool: Aligner) -> str:
    """Return the path of the aligner binary, or raise with an install hint."""
    path = shutil.which(tool)
    if path is None:
        raise PipelineError(
            f"'{tool}' was not found on PATH. Install it, e.g. `conda install -c bioconda {tool}`, "
            "or pass a pre-aligned file with --aligned"
        )
    return path


def run_aligner(records: list[Record], tool: Aligner) -> Alignment:
    """Align records and return rows in the input order."""
    if tool in _REGISTERED:
        log.info("Aligning %d sequences with %s", len(records), tool)
        result = _REGISTERED[tool](records)
        rows = dict(zip(result.ids, result.rows))
        missing = [r.id for r in records if r.id not in rows]
        if missing:
            raise PipelineError(f"{tool} output is missing sequences: {', '.join(missing)}")
        return alignment_from_records([Record(r.id, rows[r.id]) for r in records])
    binary = check_binary(tool)
    with tempfile.TemporaryDirectory(prefix="vhirprimerlab_") as tmp:
        tmpdir = Path(tmp)
        in_fa = tmpdir / "in.fasta"
        out_fa = tmpdir / "out.fasta"
        in_fa.write_text("".join(f">{r.id}\n{r.seq}\n" for r in records))
        if tool == "mafft":
            cmd = [binary, "--auto", "--quiet", str(in_fa)]
        else:
            cmd = [binary, "-align", str(in_fa), "-output", str(out_fa)]
        log.info("Aligning %d sequences with %s", len(records), tool)
        log.debug("Command: %s", " ".join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            tail = "\n".join(proc.stderr.strip().splitlines()[-20:])
            raise PipelineError(f"{tool} failed with exit code {proc.returncode}:\n{tail}")
        if tool == "mafft":
            out_fa.write_text(proc.stdout)
        aligned = {r.id: r for r in load_fasta(out_fa, allow_gaps=True)}

    missing = [r.id for r in records if r.id not in aligned]
    if missing:
        raise PipelineError(f"{tool} output is missing sequences: {', '.join(missing)}")
    return alignment_from_records([aligned[r.id] for r in records])
