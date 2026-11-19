"""The form and the engine must agree on every default.

The web form (index.html) and the command line each state the default of a setting. They drifted
once: the form showed 4 for "max variants per primer" in every mode while the engine used 8 for
regions, so a user who typed 4 sent nothing and silently got 8. This test reads the defaults from
both and fails when they differ, so the next drift is caught by CI, not by a user.
"""

from __future__ import annotations

from html.parser import HTMLParser
from pathlib import Path

import pytest

from vhirprimerlab.cli import parse_args
from vhirprimerlab.web.params import SETTINGS

INDEX = Path(__file__).resolve().parent.parent / "vhirprimerlab" / "web" / "static" / "index.html"

# The smallest argv that selects each mode; nothing is read from these paths.
ARGV = {
    "tile": ["-i", "x.fa", "-n", "8", "-o", "out"],
    "regions": ["-i", "x.fa", "--amplicon-min", "150", "--amplicon-max", "400", "-o", "out"],
    "genes": ["-i", "x.fa", "--tile-region", "A:1-500", "-o", "out"],
    "check": ["check", "-i", "x.fa", "--primers", "p.tsv", "-o", "out"],
}

VOID = ("input", "option", "br", "img", "meta", "link", "hr")


class _Fields(HTMLParser):
    """Every [data-setting] control of the page: (key, modes it shows in, default, per-mode defaults)."""

    def __init__(self) -> None:
        super().__init__()
        self.fields: list[dict] = []
        self._modes: list[str | None] = []
        self._select: dict | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        a = dict(attrs)
        if tag not in VOID:
            self._modes.append(a.get("data-modes"))
        if tag == "input" and "data-setting" in a:
            self.fields.append({"key": a["data-setting"], "modes": self._current_modes(),
                                "default": a.get("value"),
                                "per_mode": {k[len("data-default-"):]: v for k, v in a.items()
                                             if k.startswith("data-default-")}})
        elif tag == "select" and "data-setting" in a:
            self._select = {"key": a["data-setting"], "modes": self._current_modes(), "default": None,
                            "per_mode": {}}
            self.fields.append(self._select)
        elif tag == "option" and self._select is not None:
            if self._select["default"] is None or "selected" in a:
                self._select["default"] = a.get("value")

    def handle_endtag(self, tag: str) -> None:
        if tag == "select":
            self._select = None
        if tag not in VOID and self._modes:
            self._modes.pop()

    def _current_modes(self) -> list[str]:
        for m in reversed(self._modes):
            if m:
                return m.split()
        return []


def _fields() -> list[dict]:
    parser = _Fields()
    parser.feed(INDEX.read_text())
    return [f for f in parser.fields if f["key"] in {k for s in SETTINGS.values() for k in s}]


def _engine_default(mode: str, key: str):
    flag = SETTINGS[mode][key][0]
    args = parse_args(ARGV[mode])[0]
    dest = {"-n": "amplicons"}.get(flag, flag.lstrip("-").replace("-", "_"))
    return getattr(args, dest)


CASES = [(f["key"], mode) for f in _fields() for mode in f["modes"]
         if mode in SETTINGS and f["key"] in SETTINGS[mode] and f["default"] not in (None, "")]


@pytest.mark.parametrize("key, mode", CASES)
def test_the_form_default_is_the_engine_default(key: str, mode: str) -> None:
    field = next(f for f in _fields() if f["key"] == key and mode in f["modes"])
    shown = field["per_mode"].get(mode, field["default"])
    engine = _engine_default(mode, key)
    if isinstance(engine, bool) or engine is None:
        pytest.skip(f"{key} has no plain default in {mode} mode")
    try:
        assert float(shown) == pytest.approx(float(engine)), f"{key} in {mode}: form {shown!r}, engine {engine!r}"
    except ValueError:  # a select: both are words
        assert str(shown) == str(engine), f"{key} in {mode}: form {shown!r}, engine {engine!r}"


def test_the_parity_test_sees_the_fields_it_is_meant_to_guard() -> None:
    keys = {key for key, _ in CASES}
    assert {"max_degeneracy", "tm_min", "min_coverage", "overlap", "aligner"} <= keys
    assert ("max_degeneracy", "regions") in CASES and ("max_degeneracy", "tile") in CASES
