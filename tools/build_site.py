"""Build the GitHub Pages site: the landing page plus the browser app under app/.

    python tools/build_site.py OUT_DIR

OUT_DIR receives site/ as is, vhirprimerlab/web/static/ under app/, and in app/engine/ the
committed engine files (primer3 wheel, MAFFT module), a freshly built vhirprimerlab wheel and
manifest.json, which tells the worker which wheels to load.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PYODIDE = "0.27.7"


def build(out: Path) -> None:
    out = out.resolve()
    if out == ROOT or ROOT / "site" in (out, *out.parents) or out in ROOT.parents:
        raise SystemExit(f"refusing to build into {out}")
    if out.exists() and any(out.iterdir()):
        if not (out / "app" / "engine" / "manifest.json").is_file():
            raise SystemExit(f"{out} is not a previous site build; choose a new or empty folder")
        shutil.rmtree(out)
    elif out.exists():
        out.rmdir()
    shutil.copytree(ROOT / "site", out)
    app = out / "app"
    shutil.copytree(ROOT / "vhirprimerlab" / "web" / "static", app, dirs_exist_ok=True)
    engine = app / "engine"
    subprocess.run([sys.executable, "-m", "pip", "wheel", "--quiet", "--no-deps",
                    "-w", str(engine), str(ROOT)], check=True)
    vhirprimerlab = sorted(engine.glob("vhirprimerlab-*.whl"))
    primer3 = sorted(engine.glob("primer3_py-*.whl"))
    if len(vhirprimerlab) != 1 or len(primer3) != 1:
        raise SystemExit(f"expected one vhirprimerlab and one primer3 wheel in {engine}")
    manifest = {"vhirprimerlab": vhirprimerlab[0].name, "primer3": primer3[0].name, "pyodide": PYODIDE}
    (engine / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    build(Path(sys.argv[1]))
