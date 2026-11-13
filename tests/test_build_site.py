from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_build_site(tmp_path: Path) -> None:
    out = tmp_path / "_site"
    subprocess.run([sys.executable, str(ROOT / "tools" / "build_site.py"), str(out)],
                   check=True, capture_output=True)
    assert (out / "index.html").is_file()                      # the landing page
    app = out / "app"
    for name in ("index.html", "app.js", "app.css", "engine-worker.js", "mafft-runner.js",
                 "anchors.js", "annotation.js", "blocks.js", "organisms.json", "example.fasta"):
        assert (app / name).is_file(), name
    engine = app / "engine"
    manifest = json.loads((engine / "manifest.json").read_text())
    assert manifest["pyodide"] == "0.27.7"
    for key in ("vhirprimerlab", "primer3"):
        assert (engine / manifest[key]).is_file(), key
    assert manifest["vhirprimerlab"].startswith("vhirprimerlab-") and manifest["vhirprimerlab"].endswith(".whl")
    for name in ("disttbfast", "tbfast", "dndpre", "dvtditr"):    # what `mafft --auto` runs
        assert (engine / f"mafft-{name}.js").is_file() and (engine / f"mafft-{name}.wasm").is_file(), name


def test_build_site_replaces_an_old_build(tmp_path: Path) -> None:
    out = tmp_path / "_site"
    (out / "app" / "engine").mkdir(parents=True)
    (out / "app" / "engine" / "vhirprimerlab-0.0.1-py3-none-any.whl").write_text("old")
    (out / "app" / "engine" / "manifest.json").write_text("{}")  # marks a previous build
    subprocess.run([sys.executable, str(ROOT / "tools" / "build_site.py"), str(out)],
                   check=True, capture_output=True)
    wheels = sorted(p.name for p in (out / "app" / "engine").glob("vhirprimerlab-*.whl"))
    assert wheels == [json.loads((out / "app" / "engine" / "manifest.json").read_text())["vhirprimerlab"]]


def test_build_site_refuses_to_replace_a_folder_that_is_not_a_build(tmp_path: Path) -> None:
    precious = tmp_path / "docs"
    precious.mkdir()
    (precious / "notes.md").write_text("keep me")
    proc = subprocess.run([sys.executable, str(ROOT / "tools" / "build_site.py"), str(precious)],
                          capture_output=True, text=True)
    assert proc.returncode != 0 and "not a previous site build" in proc.stderr
    assert (precious / "notes.md").read_text() == "keep me"


def test_mafft_licence_ships_with_the_module(tmp_path: Path) -> None:
    out = tmp_path / "_site"
    subprocess.run([sys.executable, str(ROOT / "tools" / "build_site.py"), str(out)],
                   check=True, capture_output=True)
    text = (out / "app" / "engine" / "MAFFT-LICENSE.txt").read_text()
    assert "Kazutaka Katoh" in text and "Redistribution" in text


def test_pyodide_is_served_from_the_site_and_matches_the_pinned_hashes(tmp_path: Path) -> None:
    import hashlib

    out = tmp_path / "_site"
    subprocess.run([sys.executable, str(ROOT / "tools" / "build_site.py"), str(out)],
                   check=True, capture_output=True)
    pins = json.loads((ROOT / "tools" / "pyodide-0.27.7.json").read_text())["sha256"]
    assert "pyodide.asm.wasm" in pins and len(pins) >= 9
    for name, digest in pins.items():
        data = (out / "app" / "engine" / "pyodide" / name).read_bytes()
        assert hashlib.sha256(data).hexdigest() == digest, name


def test_no_page_or_script_loads_anything_from_a_third_party() -> None:
    """The pages promise that nothing leaves the browser: no CDN, no font service.

    Only resource loaders are checked (tags, CSS url(), importScripts), not links the reader may click.
    """
    import re

    static = ROOT / "vhirprimerlab" / "web" / "static"
    pages = [ROOT / "site" / "index.html", static / "index.html"]
    loaders = (r'<(?:script|link|img|source|iframe|video)\b[^>]*\b(?:src|href)="https?://([\w.-]+)',
               r"url\(\s*['\"]?https?://([\w.-]+)", r"importScripts\([^)]*https?://([\w.-]+)")
    sources = [*pages, *(ROOT / "site").glob("*.css"), *static.glob("*.css"), static / "engine-worker.js"]
    for path in sources:
        for pattern in loaders:
            for host in re.findall(pattern, path.read_text()):
                raise AssertionError(f"{path.name} loads from {host}")
    assert "cdn.jsdelivr.net" not in (static / "engine-worker.js").read_text()
    for html in pages:
        assert "Content-Security-Policy" in html.read_text(), html
