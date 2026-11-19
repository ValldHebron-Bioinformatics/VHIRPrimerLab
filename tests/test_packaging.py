"""What `pip install` ships: the version, every module, the web assets and the type marker.

The browser app installs the wheel this builds, so a file missing here is missing on the site.
"""

from __future__ import annotations

import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

import vhirprimerlab

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def wheel(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("wheel")
    subprocess.run([sys.executable, "-m", "pip", "wheel", "--quiet", "--no-deps", "-w", str(out), str(ROOT)],
                   check=True, capture_output=True)
    (built,) = out.glob("vhirprimerlab-*.whl")
    return built


def test_the_wheel_carries_the_version_of_the_package(wheel: Path) -> None:
    assert wheel.name.startswith(f"vhirprimerlab-{vhirprimerlab.__version__}-")


def test_the_wheel_has_every_module_and_the_web_assets(wheel: Path) -> None:
    names = set(zipfile.ZipFile(wheel).namelist())
    modules = {f"vhirprimerlab/{p.name}" for p in (ROOT / "vhirprimerlab").glob("*.py")}
    modules |= {f"vhirprimerlab/web/{p.name}" for p in (ROOT / "vhirprimerlab" / "web").glob("*.py")}
    assert modules <= names, sorted(modules - names)
    static = {f"vhirprimerlab/web/static/{p.relative_to(ROOT / 'vhirprimerlab/web/static').as_posix()}"
              for p in (ROOT / "vhirprimerlab/web/static").rglob("*") if p.is_file()}
    assert static <= names, sorted(static - names)


def test_the_wheel_marks_the_package_as_typed(wheel: Path) -> None:
    assert "vhirprimerlab/py.typed" in zipfile.ZipFile(wheel).namelist()
