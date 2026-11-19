from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from tests.synthetic import make_family, write_fasta

from vhirprimerlab import thermo

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")

# The server-mode tests import FastAPI at module level, so deselecting them is not enough on a
# machine without the web extra: they must not be collected at all unless asked for (-m server).
SERVER_TESTS = {"test_web_api.py", "test_web_jobs.py"}


def pytest_ignore_collect(collection_path: Path, config: pytest.Config) -> bool | None:
    if collection_path.name in SERVER_TESTS:
        return "not server" in (config.getoption("markexpr") or "")
    return None


@pytest.fixture
def family_base() -> str:
    return make_family()[0]


@pytest.fixture
def family(tmp_path: Path) -> Path:
    return write_fasta(tmp_path / "family.fasta", make_family()[1])


@pytest.fixture(autouse=True)
def _fresh_thermo_cache():
    """Tests replace primer3 with stand-ins; an answer cached from a real call must not leak."""
    thermo.clear_caches()
    yield
    thermo.clear_caches()
