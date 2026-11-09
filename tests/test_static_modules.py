"""Runs the Node tests of the plain-JavaScript modules in web/static (no Pyodide, no built site)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent / "browser"

pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="needs node")


def test_static_module_node_tests() -> None:
    proc = subprocess.run(
        ["node", "--test", "anchors.test.mjs"], cwd=HERE,
        capture_output=True, text=True, timeout=300,
    )
    assert proc.returncode == 0, proc.stdout[-4000:] + proc.stderr[-2000:]
