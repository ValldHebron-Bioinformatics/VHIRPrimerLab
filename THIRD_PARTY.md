# Third-party software and data

VHIRPrimerLab itself is under the GNU General Public License, version 3 or later (`LICENSE`).
The web app on GitHub Pages also ships, or downloads when the site is built, the components below.
Each keeps its own licence. "Checked" means read from the package's own metadata or licence file
in this repository; "not checked" means I could not confirm it from here, so confirm it before
relying on it (for example before redistributing the site elsewhere).

| Component | Used for | Licence | Checked |
|---|---|---|---|
| [primer3-py](https://github.com/libnano/primer3-py) (wraps [primer3](https://github.com/primer3-org/primer3)) | Tm, hairpin and dimer calculations; the WebAssembly build in `site/app/engine/` is version 2.2.0 | GPL-2.0-or-later | Checked in the bundled 2.2.0 wheel: 39 of its source files say "either version 2 of the License, or (at your option) any later version", so it can be combined with this project under GPL-3.0. Its package metadata only says "GPLv2", which is a loose label, not a restriction to version 2 only. The 2.3.1 metadata (used by the command line) says GPL-2.0-or-later |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) 7.526 | Alignment in the browser (WebAssembly) and on the command line (installed separately) | BSD-style | Checked: `site/app/engine/MAFFT-LICENSE.txt` |
| [Pyodide](https://pyodide.org/) 0.27.7 | Python in the browser | MPL-2.0 | Not checked |
| [Biopython](https://biopython.org/) | Sequence input and output | Biopython License Agreement (BSD-style) | Package metadata: `LicenseRef-Biopython-License-Agreement` |
| [NumPy](https://numpy.org/) | Specificity scan | BSD-3-Clause (with some bundled permissive components) | Package metadata |
| [FastAPI](https://fastapi.tiangolo.com/), [Uvicorn](https://www.uvicorn.org/) | Optional local server (`vhirprimerlab web`) | MIT; BSD-3-Clause | Package metadata |
| [Raleway](https://github.com/impallari/Raleway) and [Geist Mono](https://github.com/vercel/geist-font) fonts | Page typography, served from `site/fonts/` and `vhirprimerlab/web/static/fonts/` | SIL Open Font License 1.1 | Stated in `fonts.css`; the licence text is not included in this repository |

The Pyodide files are downloaded at build time and verified against the SHA-256 hashes in
`tools/pyodide-0.27.7.json` (see `site/app/engine/README.md`).
