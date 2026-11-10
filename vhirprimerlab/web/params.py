"""Web settings: limits, and translating form settings into vhirprimerlab arguments."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path

from vhirprimerlab.io import load_fasta
from vhirprimerlab.lengths import ungapped_length
from vhirprimerlab.models import PipelineError

MODES = ("tile", "regions", "check")


@dataclass(frozen=True)
class Limits:
    max_upload_mb: float = 20
    max_seqs: int = 500
    max_seq_len: int = 35_000
    max_amplicons: int = 200
    max_queue: int = 20
    workers: int = 2
    timeout_s: float = 900
    retention_s: float = 86_400

    @classmethod
    def from_env(cls) -> Limits:
        """Defaults, overridden by VHIRPRIMERLAB_<FIELD> environment variables."""
        values = {}
        for name, default in cls().__dict__.items():
            raw = os.environ.get(f"VHIRPRIMERLAB_{name.upper()}")
            if raw is not None:
                values[name] = type(default)(raw)
        return cls(**values)


# setting -> (flag, kind); kind is int, float, str, "aligner", "bool", "regions"
_SHARED = {
    "max_mismatches": ("--max-mismatches", int),
    "three_prime_window": ("--three-prime-window", int),
    "min_coverage": ("--min-coverage", float),
    "aligner": ("--aligner", "aligner"),
    "reference_id": ("--reference-id", str),
}
_DESIGN = {
    "tm_min": ("--tm-min", float),
    "tm_opt": ("--tm-opt", float),
    "tm_max": ("--tm-max", float),
    "max_tm_diff": ("--max-tm-diff", float),
    "primer_len_min": ("--primer-len-min", int),
    "primer_len_max": ("--primer-len-max", int),
    "gc_min": ("--gc-min", float),
    "gc_max": ("--gc-max", float),
    "max_degeneracy": ("--max-degeneracy", int),
    **_SHARED,
}
SETTINGS: dict[str, dict[str, tuple[str, object]]] = {
    "tile": {
        "amplicons": ("-n", int),
        "overlap": ("--overlap", int),
        "max_start_loss": ("--max-start-loss", int),
        "max_end_loss": ("--max-end-loss", int),
        "degenerate": ("--no-degenerate", "bool_off"),
        **_DESIGN,
    },
    "regions": {
        "regions": ("--region", "regions"),
        "amplicon_min": ("--amplicon-min", int),
        "amplicon_max": ("--amplicon-max", int),
        "allow_degenerate": ("--allow-degenerate", "bool_on"),
        **_DESIGN,
    },
    "check": dict(_SHARED),
}


def _number(key: str, value: object, kind: type) -> str:
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f"setting {key} must be a number")
    try:
        number = kind(value) if kind is float else int(str(value))
        if not math.isfinite(number):
            raise ValueError
    except ValueError:
        what = "a whole number" if kind is int else "a number"
        raise ValueError(f"setting {key} must be {what}") from None
    return str(number)


def build_argv(mode: str, params: dict[str, object]) -> list[str]:
    """vhirprimerlab flags for the given settings (only whitelisted settings are accepted)."""
    if mode not in SETTINGS:
        raise ValueError(f"unknown mode: {mode}")
    allowed = SETTINGS[mode]
    argv: list[str] = []
    for key, value in params.items():
        if key not in allowed:
            raise ValueError(f"unknown setting: {key}")
        if value is None or value == "":
            continue
        flag, kind = allowed[key]
        if kind in ("bool_on", "bool_off"):
            if not isinstance(value, bool):
                raise ValueError(f"setting {key} must be true or false")
            if value == (kind == "bool_on"):
                argv.append(flag)
        elif kind == "regions":
            if not isinstance(value, list) or not all(isinstance(r, str) for r in value):
                raise ValueError("setting regions must be a list of START-END or NAME:START-END")
            for region in value:
                if region.strip():
                    argv += [flag, region.strip()]
        elif kind == "aligner":
            if value not in ("mafft", "muscle"):
                raise ValueError("setting aligner must be mafft or muscle")
            argv += [flag, str(value)]
        elif kind is str:
            if not isinstance(value, str):
                raise ValueError(f"setting {key} must be text")
            argv += [flag, value]
        else:
            argv += [flag, _number(key, value, kind)]  # type: ignore[arg-type]
    return argv


def validate_upload(path: Path, limits: Limits) -> tuple[int, int]:
    """(number of sequences, longest length) of an uploaded FASTA, or ValueError."""
    try:
        records = load_fasta(Path(path))
    except PipelineError as exc:
        raise ValueError(str(exc).replace(str(path), "your file")) from None
    if len(records) > limits.max_seqs:
        raise ValueError(
            f"Too many sequences: {len(records)}. The web app accepts up to {limits.max_seqs}; "
            "use the command line for larger sets."
        )
    longest = max(records, key=lambda r: ungapped_length(r.seq))
    n = ungapped_length(longest.seq)
    if n > limits.max_seq_len:
        raise ValueError(
            f"Sequence {longest.id} is {n:,} bp. The web app accepts up to "
            f"{limits.max_seq_len:,} bp per sequence; use the command line for longer genomes."
        )
    return len(records), n
