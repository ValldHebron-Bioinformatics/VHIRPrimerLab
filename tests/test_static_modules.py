"""Runs the Node tests of the plain-JavaScript modules in web/static (no Pyodide, no built site).

Needs Node 20.19 / 22.7 or newer, which decides a .js file is ESM by its syntax. A
`{"type": "module"}` package.json in web/static would say so outright, but that folder also
holds the classic scripts the page's worker loads with importScripts() and the Node tests with
require() (engine.js, engine-worker.js, mafft-runner.js, and the MAFFT WebAssembly modules the
site build puts beside them), and "type" applies to a whole subtree: declaring it there makes
require() of those return an empty namespace, and tests/browser/mafft.test.mjs and
engine.test.mjs fail. Separating the two kinds would mean renaming files the browser loads.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent / "browser"

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="needs node")


def test_static_module_node_tests() -> None:
    proc = subprocess.run(
        ["node", "--test", "anchors.test.mjs", "bestmatch.test.mjs", "annotation.test.mjs", "blocks.test.mjs",
         "app.test.mjs"], cwd=HERE,
        capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-2000:]
