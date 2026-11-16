from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from tests.synthetic import make_family, write_fasta

from vhirprimerlab import thermo

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")


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
