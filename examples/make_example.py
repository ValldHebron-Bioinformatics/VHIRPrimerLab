"""Write the synthetic example families.

* examples/family.fasta: 10 short sequences (region design and check examples)
* examples/genomes.fasta and vhirprimerlab/web/static/example.fasta: 24 genomes of 9 kb
  (tiling example and the web app's "Try the example")
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from tests.synthetic import make_family, make_genome_family, write_fasta  # noqa: E402

GENOMES = dict(length=9000, n_seqs=24, seed=11, hyper=(4200, 4500))

if __name__ == "__main__":
    outs = [write_fasta(ROOT / "examples" / "family.fasta", make_family()[1])]
    genomes = make_genome_family(**GENOMES)
    for path in (ROOT / "examples" / "genomes.fasta", ROOT / "vhirprimerlab" / "web" / "static" / "example.fasta"):
        path.parent.mkdir(parents=True, exist_ok=True)
        outs.append(write_fasta(path, genomes))
    for out in outs:
        print(f"Wrote {out}")
