# Browser engine files

Compiled parts of the vhirprimerlab browser engine, loaded by `app/engine-worker.js`.
Rebuild them with `bash tools/build_wasm.sh` (from the repository root).

| File | What it is | Source and licence |
|---|---|---|
| `primer3_py-2.2.0-cp312-cp312-pyodide_2024_0_wasm32.whl` | primer3-py 2.2.0 built for Pyodide 0.27 (Emscripten 3.1.58) | [primer3-py](https://github.com/libnano/primer3-py), GPL-2.0 |
| `mafft-<program>.js`, `mafft-<program>.wasm` | MAFFT 7.526 `disttbfast`, `tbfast`, `dndpre` and `dvtditr`, the programs `mafft --auto` runs, as WebAssembly modules; factories `createMafft_<program>` | [MAFFT](https://mafft.cbrc.jp/alignment/software/), BSD (see `MAFFT-LICENSE.txt`) |

The vhirprimerlab wheel and `manifest.json` are added at deploy time by `tools/build_site.py`.
Pyodide, Biopython, NumPy and micropip are downloaded at build time from the jsDelivr CDN, verified against the SHA-256 hashes in `tools/pyodide-0.27.7.json` and served from `app/engine/pyodide/`.
