"""In CI the Node suites must run, never skip.

test_browser_engine.py and test_static_modules.py skip when Node, Pyodide or the MAFFT WebAssembly
files are missing, which is right on a laptop and wrong in CI, where a skip looks like a pass.
CI sets REQUIRE_BROWSER_TESTS=1; these tests then fail loudly when a prerequisite is absent.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

pytestmark = pytest.mark.skipif(
    not os.environ.get("REQUIRE_BROWSER_TESTS"), reason="only enforced in CI (REQUIRE_BROWSER_TESTS=1)")


def test_node_is_installed() -> None:
    assert shutil.which("node"), "node is missing: the browser tests would be skipped"


def test_pyodide_is_installed_for_the_node_tests() -> None:
    assert (ROOT / "tests" / "browser" / "node_modules" / "pyodide").is_dir(), \
        "run `npm ci` in tests/browser: the engine tests would be skipped"


@pytest.mark.parametrize("name", ["disttbfast", "tbfast", "dndpre", "dvtditr"])
def test_mafft_wasm_is_in_the_repository(name: str) -> None:
    assert (ROOT / "site" / "app" / "engine" / f"mafft-{name}.wasm").is_file()
