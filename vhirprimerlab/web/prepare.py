"""Check a web request and turn it into vhirprimerlab arguments (server and browser engines).

No web framework here: the browser worker imports this module under Pyodide.
"""

from __future__ import annotations

import contextlib
import io
import tempfile
from pathlib import Path

from vhirprimerlab.check import load_assays
from vhirprimerlab.cli import MAX_BACKGROUND_BASES, parse_args, validate_args
from vhirprimerlab.io import load_fasta
from vhirprimerlab.models import PipelineError
from vhirprimerlab.regions import parse_region_spec
from vhirprimerlab.tiling import amplicons_for, max_amplicons, target_length
from vhirprimerlab.web.params import MODES, Limits, build_argv, validate_upload


def _usage_message(stderr: str) -> str:
    lines = [ln for ln in stderr.splitlines() if "error: " in ln]
    return lines[-1].split("error: ", 1)[1] if lines else "invalid settings"


def prepare_job(
    mode: str,
    settings: object,
    fasta: bytes,
    primers_text: str,
    limits: Limits,
    aligners: list[str],
    background_text: str = "",
) -> tuple[list[str], dict[str, bytes]]:
    """(argv without -o, files for the job directory), or ValueError with a user-facing message.

    background_text is an optional FASTA of sequences the primers must not bind (tiling and check).
    """
    if mode not in MODES:
        raise ValueError(f"unknown mode: {mode}")
    if not isinstance(settings, dict):
        raise ValueError("settings must be a JSON object")
    max_bytes = int(limits.max_upload_mb * 1024 * 1024)
    if len(fasta) > max_bytes:
        raise ValueError(f"The file is larger than {limits.max_upload_mb:g} MB; "
                         "use the command line for larger inputs.")
    if not fasta.strip():
        raise ValueError("Add your sequences (FASTA) first.")

    regions = settings.get("tile_regions")
    if mode == "genes" and not (isinstance(regions, list) and any(
            isinstance(r, str) and r.strip() for r in regions)):
        raise ValueError("Pick at least one region on the genome first.")

    background = bool(background_text.strip())
    if background and mode not in ("tile", "genes", "check"):
        raise ValueError("Sequences to avoid are scanned when tiling and in check mode only.")
    if background and len(background_text.encode()) > max_bytes:
        raise ValueError(f"The background file is larger than {limits.max_upload_mb:g} MB; "
                         "use the command line for larger inputs.")

    flags = build_argv(mode, settings)
    argv = (["check"] if mode == "check" else []) + ["-i", "input.fasta"]
    argv += (["--primers", "primers.tsv"] if mode == "check" else []) + flags
    argv += ["--background", "background.fasta"] if background else []
    stderr = io.StringIO()
    try:
        with contextlib.redirect_stderr(stderr):
            args, _, _ = parse_args([*argv, "-o", "out"])  # the runner adds -o itself
        validate_args(args)
    except SystemExit:
        raise ValueError(_usage_message(stderr.getvalue())) from None
    except PipelineError as exc:
        raise ValueError(str(exc)) from None
    if args.aligner not in aligners:
        raise ValueError(f"{args.aligner.upper()} is not installed on this server; use "
                         f"{' or '.join(a.upper() for a in aligners) or 'the command line'}.")

    with tempfile.TemporaryDirectory() as tmp:
        upload = Path(tmp) / "input.fasta"
        upload.write_bytes(fasta)
        _, longest = validate_upload(upload, limits)
        files = {"input.fasta": fasta}
        if background:
            path = Path(tmp) / "background.fasta"
            path.write_text(background_text)
            try:
                records = load_fasta(path)
            except PipelineError as exc:
                raise ValueError("The background sequences: " + str(exc).replace(str(path), "your file")) from None
            bases = sum(len(r.seq) for r in records)
            if bases > MAX_BACKGROUND_BASES:
                raise ValueError(f"The background sequences hold {bases:,} bases; the web app scans at most "
                                 f"{MAX_BACKGROUND_BASES:,}. Keep only what the primers could really meet, "
                                 "or use the command line.")
            files["background.fasta"] = background_text.encode()
        if mode == "tile":
            n = args.amplicons
            if n > limits.max_amplicons:
                raise ValueError(f"The web app tiles up to {limits.max_amplicons} amplicons.")
            most = max_amplicons(longest, args.overlap, args.primer_len_max)
            if n > most:
                raise ValueError(f"N too large: amplicons would be ~"
                                 f"{target_length(longest, n, args.overlap)} bp; use at most {most}")
        if mode == "genes":
            total = 0
            for i, spec in enumerate(args.tile_region, start=1):
                try:
                    region = parse_region_spec(spec, i)
                except PipelineError as exc:
                    raise ValueError(str(exc)) from None
                if region.end > longest:
                    raise ValueError(f"Region {region.name} ends at {region.end:,}, past the end "
                                     f"of your longest sequence ({longest:,} bp).")
                total += amplicons_for(region.end - region.start + 1, args.amplicon_len, args.overlap)
            if total > limits.max_amplicons:
                raise ValueError(f"These regions need about {total:,} amplicons; the web app tiles "
                                 f"up to {limits.max_amplicons}. Pick fewer regions or lengthen the amplicons.")
        if mode == "check":
            if not primers_text.strip():
                raise ValueError("Paste your primer pairs (name, forward, reverse) first.")
            primers = Path(tmp) / "primers.tsv"
            primers.write_text(primers_text)
            try:
                load_assays(primers, args.three_prime_window)
            except PipelineError as exc:
                raise ValueError(str(exc).replace(str(primers), "primers")) from None
            files["primers.tsv"] = primers_text.encode()
    return argv, files
