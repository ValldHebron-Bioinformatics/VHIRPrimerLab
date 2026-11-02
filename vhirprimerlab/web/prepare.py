"""Check a web request and turn it into vhirprimerlab arguments (server and browser engines).

No web framework here: the browser worker imports this module under Pyodide.
"""

from __future__ import annotations

import contextlib
import io
import tempfile
from pathlib import Path

from vhirprimerlab.check import load_assays
from vhirprimerlab.cli import parse_args, validate_args
from vhirprimerlab.models import PipelineError
from vhirprimerlab.tiling import max_amplicons, target_length
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
) -> tuple[list[str], dict[str, bytes]]:
    """(argv without -o, files for the job directory), or ValueError with a user-facing message."""
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

    flags = build_argv(mode, settings)
    argv = (["check"] if mode == "check" else []) + ["-i", "input.fasta"]
    argv += (["--primers", "primers.tsv"] if mode == "check" else []) + flags
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
        if mode == "tile":
            n = args.amplicons
            if n > limits.max_amplicons:
                raise ValueError(f"The web app tiles up to {limits.max_amplicons} amplicons.")
            most = max_amplicons(longest, args.overlap, args.primer_len_max)
            if n > most:
                raise ValueError(f"N too large: amplicons would be ~"
                                 f"{target_length(longest, n, args.overlap)} bp; use at most {most}")
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
