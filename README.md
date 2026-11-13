# VHIRPrimerLab

Design PCR primers that work across a family of related sequences (strains of a virus,
variants of a gene, …). Give it your genomes and the number of amplicons you want: it tiles
the whole genome with that many overlapping amplicons in two pools, choosing primer sites that
match as many of your sequences as possible. Every primer is checked against every sequence,
not just a consensus.

```bash
vhirprimerlab -i genomes.fasta -n 10 -o out
```

The command is `vhirprimerlab` (lower case).

Prefer not to install anything? See [Run it from the web](#run-it-from-the-web).

## Installation

```bash
conda create -n primers -c conda-forge -c bioconda python=3.11 biopython "primer3-py>=2" pytest mafft muscle
conda activate primers
pip install -e .
```

Python ≥ 3.10 is required. Python dependencies are listed in `requirements.txt`; MAFFT (default)
or MUSCLE 5 must be on `PATH` unless you pass a pre-aligned file with `--aligned`.

## Quick start: tile a genome

```bash
python examples/make_example.py                # writes examples/genomes.fasta (24 genomes, 9 kb)
vhirprimerlab -i examples/genomes.fasta -n 8 -o out
```

This aligns the genomes, then picks 8 overlapping amplicons (~1.2 kb each, at least 50 bp
overlap) covering the genome end to end. Odd amplicons go in pool 1 and even ones in pool 2,
so overlapping products never share a tube. It writes:

- `primers.tsv`: the primers to order, one row per oligo with just `name` and `sequence`
  (`VHIRPrimerLab_1_LEFT`, `VHIRPrimerLab_1_RIGHT`, …), ready for a supplier's bulk upload.
- `amplicons.tsv`: one row per amplicon with every detail: pool, primer names and
  sequences, and for each primer its length, Tm, GC %, number of degenerate bases and number of
  variants (for degenerate primers also the Tm range over the variants); then reference
  coordinates, overlap with the next amplicon, coverage and the sequences each amplicon misses,
  with the reason.
- `scheme.bed`: the primer scheme in ARTIC format (`chrom start end name pool strand
  sequence`), ready for amplicon-sequencing pipelines.
- `report.txt`: settings, input and alignment statistics, the tiling (each amplicon with both
  primers: sequence, length, Tm, GC % and degenerate bases), bases left uncovered at the genome
  ends, possible primer dimers within a pool, and the log.
- `alignment.fasta`: the alignment used.

### How the amplicons are chosen

1. Every position of the consensus is tried as a primer 3′ end; primers outside the Tm, GC and
   length limits, or with hairpins, self-dimers or a G/C-heavy 3′ end, are dropped. Within
   `--max-start-loss` of the genome start and `--max-end-loss` of its end, where the sequence
   is often AT-rich, primers may be up to 32 nt long and have any GC %, so that a primer with
   the right Tm can still reach the end.
2. The most conserved primers are kept and checked against every sequence with the matching
   rules below. Degenerate (IUPAC) bases are added where they rescue sequences, up to
   `--max-degeneracy` variants per primer, as long as the primer's GC % (averaged over its
   variants) stays within `--gc-min`/`--gc-max` and no variant forms a hairpin, a self-dimer
   (also between two variants of the mix) or a G/C-heavy 3′ end.
3. A dynamic-programming search picks exactly N amplicons. The first insert must start within
   `--max-start-loss` (50) bases of the genome start, and the last one end within
   `--max-end-loss` (50) bases of its end; the two limits are separate, since one end is often
   harder than the other. Within that, it maximises the worst amplicon's coverage, then the
   total coverage, then how many sequences each primer matches on its own, then how close the
   tiling reaches both genome ends, then primer quality. Sequences that stop short of a genome
   end count as not covered by the end amplicon, so an alignment with incomplete ends shows a
   low coverage there. If no primer reaches an end, the tool says which limit to raise and to
   what: RSV, for one, has no primer with a usable Tm in its last ~270 bases. Amplicon lengths may vary ±20 % around the target; if no tiling fits, the
   tool retries at ±35 % and ±50 % and says so.
4. Primers in the same pool are checked for primer-dimers (ΔG below −9 kcal/mol between any
   variant of one primer and any variant of the other). A primer in a dimer is swapped for
   another when the tiling can be chosen again without it and the weakest amplicon keeps
   at least as many sequences; dimers that cannot be removed this way are listed in
   `report.txt`.

### Tiling chosen spans

To tile only parts of the genome (say, two genes) give each span with `--tile-region` instead of `-n`:

```bash
vhirprimerlab -i genomes.fasta -o out \
    --tile-region F:1001-3000 --tile-region SH:4001-4600 --amplicon-len 700
```

Each span is tiled on its own with amplicons of about `--amplicon-len`, overlapping by at least
`--overlap`; a span that fits in one amplicon gets one. Amplicons are numbered and pooled across
the whole scheme and named after their span (`F_1_LEFT`, `F_2_LEFT`, …; `SH_LEFT` for a
single-amplicon span). `--max-start-loss` and `--max-end-loss` apply at the start and end of every
span, and primer dimers are checked across the whole scheme.

### Tiling options

| Option | Default | Description |
|---|---|---|
| `-n/--amplicons` | — | number of amplicons |
| `--tile-region` | — | `[NAME:]START-END`, 1-based reference coordinates; tile only this span, repeatable, one scheme for all spans (replaces `-n`; cannot be combined with `-n`, `--region` or `--regions-file`) |
| `--amplicon-len` | 400 | target amplicon length with `--tile-region` (bp); each span gets as many amplicons as it needs |
| `--overlap` | 50 | minimum overlap between neighbouring amplicons (bp, between the primers) |
| `--max-start-loss` | 50 | most genome bases left outside the first insert (bp, primer included) |
| `--max-end-loss` | 50 | most genome bases left outside the last insert (bp, primer included) |
| `--tm-min/--tm-opt/--tm-max` | 57 / 60 / 63 | primer Tm (°C) |
| `--max-tm-diff` | 2 | max Tm difference within a pair (°C) |
| `--primer-len-min/--primer-len-max` | 18 / 25 | primer length |
| `--gc-min/--gc-max` | 40 / 60 | primer GC % |
| `--max-mismatches` | 2 | max mismatches per primer |
| `--max-ambiguous` | 2 | `N`/ambiguous bases in a sequence forgiven per primer (not at the 3′ end) |
| `--three-prime-window` | 5 | 3′ bases that must match exactly |
| `--no-degenerate` | off | plain primers only |
| `--max-degeneracy` | 4 | max variants per degenerate primer |
| `--min-coverage` | 95 | coverage target (%) |
| `-i`/`--aligned`, `--aligner`, `--reference-id`, `--length-tolerance`, `--outlier-action`, `-v` | | as below |

Exit codes: `0` every amplicon reaches `--min-coverage`, `2` at least one does not (outputs are
still written), `1` error (for example too many amplicons for the genome length; the message
says the most that fit).

## Run it from the web

Open the web app at **https://valldhebron-bioinformatics.github.io/VHIRPrimerLab/app/**: drop your
FASTA, choose the number of amplicons, and get the primers with a coverage map and downloads.
A "Check my primers" tab tests pairs you already have.

The analysis runs **in your browser**: VHIRPrimerLab, primer3 and MAFFT are compiled to
WebAssembly, so your sequences are never uploaded and nothing needs installing. The first run
loads about 20 MB of tools, which the browser keeps for later visits. It gives the same
results as the command line: the design steps are the same code, and the alignment runs the
same MAFFT programs (version 7.526) with the same method choice as `mafft --auto` (L-INS-i,
FFT-NS-i or FFT-NS-2, depending on the number and length of the sequences). For very large
sets (hundreds of long genomes) the command line is faster.

### Run it on your own server

The same page can run jobs on a server instead (useful for large inputs or a shared lab
machine). On your own computer:

```bash
pip install -e ".[web]"
vhirprimerlab web            # opens http://127.0.0.1:8000
```

Or with Docker (MAFFT included):

```bash
docker build -t vhirprimerlab-web .
docker run -d --restart unless-stopped -p 7860:7860 vhirprimerlab-web   # http://<host>:7860
```

Server limits can be changed with environment variables: `VHIRPRIMERLAB_MAX_SEQS` (500),
`VHIRPRIMERLAB_MAX_SEQ_LEN` (35000), `VHIRPRIMERLAB_MAX_UPLOAD_MB` (20), `VHIRPRIMERLAB_MAX_AMPLICONS` (200),
`VHIRPRIMERLAB_WORKERS` (2), `VHIRPRIMERLAB_MAX_QUEUE` (20), `VHIRPRIMERLAB_TIMEOUT_S` (900),
`VHIRPRIMERLAB_RETENTION_S` (86400) and `VHIRPRIMERLAB_JOBS_DIR`. Server uploads are deleted after
`VHIRPRIMERLAB_RETENTION_S`.

### How the browser version is built

`tools/build_wasm.sh` rebuilds the compiled parts in `site/app/engine/` (primer3-py for
Pyodide, and the MAFFT programs `mafft --auto` runs, for WebAssembly); `tools/build_site.py`
assembles the GitHub Pages site, and `.github/workflows/pages.yml` deploys it on every push to master.

## Pick by gene

The web app's "Pick by gene" mode draws an organism's genes over your own reference, so you
click the genes you want and get primers covering them instead of counting amplicons.

1. Drop your sequences. The genes are placed on whichever sequence you choose as the reference.
2. Choose the organism: one of the usual targets in the list (the list is
   `vhirprimerlab/web/static/organisms.json`; edit it to add your own), or type a name or an
   accession in the box under it. A name such as `dengue virus 2` is searched among NCBI's RefSeq
   records and the matches are listed with their accession and length: a length near the whole
   genome is the genome, a short one is a single gene. Something that is plainly an accession
   (`NC_039199.1`) is loaded directly. Records over 500,000 nt are refused. "Load a GenBank file",
   the button under the box, takes a record you already have and needs no network at all.
3. Click the genes on the drawing; drag the handles of a picked span to adjust it. A gene the
   anchor map cannot place on your reference is greyed and cannot be clicked.
4. Press Design. Each span is tiled on its own and the whole scheme comes back as usual: the
   amplicon map, the table, `primers.tsv`, `amplicons.tsv`, `scheme.bed` and `report.txt`. Every
   span is padded around the gene you clicked, so the gene lands inside the inserts rather than
   under a primer; where it sits at the very edge of your reference there is no room for that
   padding, and the list of spans says how many of its bases may fall outside the amplicons.
   "Cover the whole genome instead", under the list, tiles the whole reference as the design page
   does.

Your browser asks NCBI for the annotation and keeps it for next time; only the name or accession
you type is sent, never your sequences, and the server makes no network call. A GenBank file needs
no network at all, and neither does a list that fails to load.

## Advanced: target regions

Instead of tiling the whole genome you can design pairs inside chosen regions
(`--region`, `--amplicon-min`, `--amplicon-max`), or check primers you already have
(`vhirprimerlab check`). Region pairs are picked by primer3 with the same 3′ rule as tiling: at
most 3 G/C in the last 5 bases.

### Example

```bash
python examples/make_example.py          # writes examples/family.fasta (10 synthetic sequences)
vhirprimerlab -i examples/family.fasta -o out \
    --region cons:61-540 --region hv:701-1100 \
    --amplicon-min 100 --amplicon-max 300 --allow-degenerate
```

The example family has a conserved block (61–540) and a hypervariable block (701–1100). The run
finds pairs with 100 % coverage in `cons`, reports the best available (30 %) pairs for `hv` with
`below_threshold=True`, explains why in `out/report.txt`, and exits with code 2.

`python -m vhirprimerlab …` works the same as `vhirprimerlab …`.

### Checking existing primers

`vhirprimerlab check` tests primer pairs you already have (a published assay, pairs from an earlier
run, pairs from another tool) against your sequences, without designing new ones. Each primer is
placed where it fits the consensus best, then the pair is scored on every sequence with exactly
the same matching rules as design mode.

List the pairs in a TSV, CSV (comma or semicolon) or Excel `.xlsx` file (first sheet), written
5′→3′ as you would order them (IUPAC codes allowed; a header row, blank rows and `#` comments are
optional), one pair per row:

```
name	forward	reverse
cons	CAGTCCCCTACGTCGCATATC	TAGTTGGGCGTCCTGAGAAAC
hv	CTCCAGTGGMGGTCTTGGTKG	GCAYKCATCCMAGTAAGCACAC
```

or one primer per row (`name, sequence`), like the `primers.tsv` order sheet that design mode
writes. Names ending in `_LEFT`/`_RIGHT` (or `_F`/`_R`, `_FWD`/`_REV`, `_FORWARD`/`_REVERSE`) are
joined into pairs, so `VHIRPrimerLab_3_LEFT` and `VHIRPrimerLab_3_RIGHT` become pair
`VHIRPrimerLab_3`. You can check a design straight from its order sheet:

```bash
vhirprimerlab check -i examples/family.fasta --primers assays.tsv -o out_check
vhirprimerlab check -i examples/family.fasta --primers out/primers.tsv -o out_check
```

Every primer is also checked for hairpins and self-dimers (every variant of a degenerate primer,
Tm above 47 °C), for more than 3 G/C in its last 5 bases, and for dimers with its partner and
with the primers of the other pairs (ΔG below −9 kcal/mol). These are warnings: they do not change
coverage or the exit code. Dimers between different pairs only matter for pairs you mix in one
tube, so they are listed in `report.txt` (the strongest also as warnings).

It writes `alignment.fasta`, `report.txt` (one section per pair: each primer's length, Tm, GC % and degenerate bases,
where it binds, and why it misses sequences)
and `check.tsv`, one row per pair: `name, forward_seq, reverse_seq, forward_tm, reverse_tm,
forward_gc, reverse_gc, forward_len, reverse_len, degenerate, forward_degenerate_bases,
reverse_degenerate_bases, forward_variants, reverse_variants, forward_tm_range,
reverse_tm_range, located, forward_site_mismatches, reverse_site_mismatches, amplicon_start,
amplicon_end, amplicon_len, coverage_pct, n_covered, n_total, below_threshold, not_covered,
forward_structure, reverse_structure, pair_dimer_dg, note` (`*_structure` is `ok` or what is
wrong; `pair_dimer_dg` is forward × reverse in kcal/mol). On the example family the two pairs above give 100 % and 30 % coverage, as in the design
run.

A pair is reported as not located (0 % coverage, reason in `note`) when a primer fits two places
equally well or the reverse primer does not bind downstream of the forward one. When a primer has
no site within `--max-mismatches`, `note` says so, and also says when the primer would fit if it
were written in the other orientation (a common copy-paste slip with reverse primers). If the
reference has no bases where the amplicon falls (for example a partial genome), the amplicon
coordinates are left empty and `note` says `amplicon not on reference <id>`. Options shared
with design mode (`-i`/`--aligned`, `--aligner`, `--reference-id`, `--length-tolerance`,
`--outlier-action`, `--max-mismatches`, `--three-prime-window`, `--min-coverage`, `-v`) work
the same; design-only options such as `--region` are rejected. Exit codes: `0` every pair meets
`--min-coverage`, `2` at least one does not (outputs are still written), `1` error.

### Pipeline (region design)

1. **Load & validate** the FASTA: errors on missing/empty files, empty records, duplicate IDs
   and non-IUPAC characters. Lowercase and `U` are accepted (converted to upper-case / `T`).
2. **Length check**: sequences whose length differs from the median by more than
   `--length-tolerance` × median are dropped, kept with a warning, or abort the run
   (`--outlier-action`).
3. **Align** with MAFFT or MUSCLE, or use `--aligned` to skip this step.
4. **Consensus**: majority base per column; columns that are gaps in ≥ 50 % of sequences are left
   out of the consensus.
5. **Design** candidate pairs with primer3 on the consensus, with the whole amplicon inside each
   region.
6. **Re-score** each pair against every sequence (see *Matching rules*), optionally adding IUPAC
   degenerate bases (`--allow-degenerate`).
7. **Rank** by coverage, then primer3 penalty, and report the top `--num-pairs` per region.

For degenerate primers, GC % is the average over all variants (degenerate bases are only added
while it stays within `--gc-min`/`--gc-max`), and `*_tm_range` gives the
lowest and highest Tm over all variants. `*_degenerate_bases` counts the
IUPAC positions and `*_variants` the distinct oligos the primer stands for.

### Matching rules (all modes)

A primer matches a sequence when it has at most `--max-mismatches` mismatches **and** none in the
last `--three-prime-window` bases of its 3′ end. A pair covers a sequence when both primers match.

- A gap counts as a mismatch.
- An `N` or other ambiguous base in the sequence is *unknown*, not different: when it could be the
  primer's base (`N`, or `R` against `A`) it is forgiven, up to `--max-ambiguous` of them per primer
  (default 2); further ones count as mismatches, so a run of `N` still fails. One that cannot be
  the primer's base (`Y` against `A`) is a mismatch. Within the last `--three-prime-window` bases an
  unknown base fails the primer (reported as "ambiguous 3' base"), because it cannot be trusted
  there. `--max-ambiguous 0` counts every one of them as a mismatch.
- Each base the sequence has in addition to the primer (an insertion) counts as one mismatch.
- A degenerate primer base matches any base it stands for (`R` matches `A` and `G`).

With `--allow-degenerate`, IUPAC bases are added one at a time, each time choosing the change
that rescues the most sequences (then the one removing the most mismatches, then the one
farthest from the 3′ end), as long as the primer's fold (number of distinct oligos) stays at or
below `--max-degeneracy`. The 3′-terminal base is never made degenerate, and changes that do not
end up rescuing a sequence are undone.

### Coordinates and regions

All coordinates are 1-based, inclusive, on the reference sequence (`--reference-id`, default the
first sequence). They are converted to alignment columns internally.

- `--region START-END` or `--region NAME:START-END`, repeatable. Unnamed regions are called
  `region1`, `region2`, … in command-line order.
- `--regions-file FILE`: tab-separated `name`, `start`, `end`; an optional header line, blank lines
  and `#` comments are allowed.
- Without either, the whole reference is one region called `full`.

Regions are designed independently and may overlap.

### Region design options

| Option | Default | Description |
|---|---|---|
| `-i/--input` | — | multi-FASTA (exactly one of `-i` / `--aligned` is required) |
| `--aligned` | — | pre-aligned FASTA; skips alignment |
| `-o/--outdir` | required | output directory (created if missing) |
| `--aligner` | `mafft` | `mafft` or `muscle` |
| `--reference-id` | first sequence | reference for coordinates |
| `--region` / `--regions-file` | whole reference | target regions (see above) |
| `--amplicon-min`, `--amplicon-max` | **required** for region design | amplicon length range (bp) |
| `--tm-min/--tm-opt/--tm-max` | 57 / 60 / 63 | primer Tm (°C) |
| `--max-tm-diff` | 2 | max Tm difference within a pair (°C) |
| `--primer-len-min/--primer-len-max` | 18 / 25 | primer length |
| `--gc-min/--gc-max` | 40 / 60 | primer GC % |
| `--length-tolerance` | 0.10 | allowed deviation from median length (fraction) |
| `--outlier-action` | `drop` | `warn`, `drop` or `abort` |
| `--max-mismatches` | 2 | max mismatches per primer |
| `--max-ambiguous` | 2 | `N`/ambiguous bases in a sequence forgiven per primer (not at the 3′ end) |
| `--three-prime-window` | 5 | 3′ bases that must match exactly |
| `--min-coverage` | 95 | minimum coverage (%) |
| `--allow-degenerate` | off | allow IUPAC bases at variable positions |
| `--max-degeneracy` | 8 | max fold per degenerate primer |
| `--candidates` | 200 | primer3 pairs re-scored per region |
| `--num-pairs` | 5 | pairs reported per region |
| `-v/--verbose` | off | debug logging and tracebacks |

### Region design outputs

- `alignment.fasta` — the alignment used (after outlier removal).
- `primers.tsv` — the primers to order: the best (rank 1) pair of each region as `name` and
  `sequence` rows (`<region>_LEFT`, `<region>_RIGHT`).
- `pairs.tsv` — one row per reported pair: `region, rank, forward_seq, reverse_seq, forward_tm,
  reverse_tm, forward_gc, reverse_gc, forward_len, reverse_len, amplicon_start, amplicon_end,
  amplicon_len, coverage_pct, n_covered, n_total, degenerate, forward_tm_range, reverse_tm_range,
  penalty, below_threshold, not_covered`. Amplicon coordinates are on the reference; `rank`
  restarts at 1 for each region; `not_covered` lists `ID(reason)` for every missed sequence. Tm
  ranges are given for degenerate primers (min–max over all expansions).
- `report.txt` — parameters, input statistics, dropped outliers, alignment statistics, one section
  per region (candidates, coverage, the best pair's primers with length, Tm, GC % and degenerate
  bases, and why it misses sequences when no pair reaches `--min-coverage`), a summary, and the full log.

### Exit codes (region design)

| Code | Meaning |
|---|---|
| 0 | every region has at least one pair meeting `--min-coverage` |
| 1 | error (bad input or parameters, aligner missing or failed) |
| 2 | finished, but at least one region has no pair meeting `--min-coverage` (or no primer3 candidates); outputs are still written |

## Tests

```bash
pytest
```

Tests use a synthetic family (`tests/synthetic.py`); those that run MAFFT are skipped if it is not
installed.
