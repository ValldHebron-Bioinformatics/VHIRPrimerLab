# Genome browser: pick amplicons by gene

**Date:** 2026-11-09
**Status:** approved design, ready for an implementation plan
**Stacks on:** PR #11 (`feat/tiling-genome-ends`), which this work extends

## Why

Today a user who wants primers over one gene must open a genome record elsewhere, read the
gene's coordinates, convert them to their own reference, and type them into the regions box.
That is the step where a virologist gives up or gets the coordinates wrong.

This feature makes the genome itself the interface: pick the organism, see its genes and
proteins drawn over your own sequences, click the protein you care about, and get the primers
that cover it.

**Who it is for:** virology and parasitology labs that design primers over specific genes
(a genotyping target, a resistance locus, one protein for a new assay), not over whole genomes.

**Success:** for human metapneumovirus, a user who has never used the regions mode can go from
a FASTA to an order sheet covering the F gene without typing a single coordinate, and the gene's
boundaries on their own reference are right to within a few bases.

## Decisions taken with the user

| Question | Decision |
|---|---|
| May sequences leave the machine? | Only on an explicit click, with a notice. Phase 1 sends no sequence at all; only an accession number. |
| Which organisms? | Viruses and parasites. For parasites the input is a single gene or locus, never a whole chromosome. |
| Segmented viruses (influenza, rotavirus)? | One segment per run, so the viewer always draws a single molecule. |
| Clicking a protein that is longer than one amplicon? | A block of overlapping amplicons covering the whole protein. |
| Where does it live? | A new mode beside tiling, regions and check. The regions mode stays as it is. |
| Automatic organism detection? | Phase 2, on top of a catalogue and a name search (phase 1). |
| Contact address for the BLAST services? | A per-user setting stored in their own browser. No address is committed to the repository. |

## What was measured before designing

These three probes decided the architecture. All were run against the user's real hMPV data
(93 public Pathoplexus genomes) and the public RefSeq record `NC_039199.1`.

1. **The annotation is available and structured.** NCBI's efetch service returns INSDSeq XML
   with the organism, the record length and every feature (`gene`, `CDS`, `mat_peptide`, UTRs)
   carrying its location, gene name and product. It answers with
   `Access-Control-Allow-Origin: *`, so a browser may fetch it directly.
2. **Coordinates can be mapped without an aligner.** Shared exact 20-mers between the RefSeq
   record and the run's reference, chained into a monotone list and interpolated between
   anchors, reproduced a MAFFT pairwise alignment to within **13 nt over 13 kb** in **0.04 s**
   (N gene exact, F gene off by 13 at one edge, L gene off by 3). A MAFFT pairwise alignment of
   the same two sequences takes 0.33 s but needs an aligner, which the browser build does not
   have for this step.
3. **Automatic detection is the fragile part.** NCBI's public BLAST API answered "Temporarily
   unavailable"; EBI's REST service refused a job without a valid contact address. This is why
   detection is phase 2 and never a prerequisite.

## Scope

**In scope (phase 1)**

- Tiling restricted to one or more spans of the reference, in the engine and on the command line.
- Fetching and reading a RefSeq annotation in the browser, with a cache.
- Mapping record coordinates onto the run's reference with shared k-mers.
- A new web mode: pick the organism, see the genome drawn, click genes, adjust blocks, design.
- A small catalogue of the lab's usual targets, plus a name search.
- Loading a GenBank file by hand, so the mode works with no network and with private annotations.

**Out of scope**

- Automatic organism detection from the sequence (phase 2, its own spec).
- Whole bacterial or parasite chromosomes; records are capped at 500 kb.
- Several segments in one run.
- Zoom and pan: the view is cropped to the part the user's sequences cover.
- Changing the regions mode, the tiling mode or the check mode.

## Architecture

Four units with narrow interfaces. Only the first is Python; the rest are browser-side, which is
what keeps the server free of any network code and makes the mode work in the GitHub Pages build.

```
browser                                                   engine (Python)
┌──────────────┐  accession   ┌───────────────┐
│ catalogue /  ├─────────────►│ annotation.js │  (fetch + parse INSDSeq XML, cache)
│ name search  │              └───────┬───────┘
└──────────────┘                      │ features in record coordinates
                                      ▼
┌──────────────┐  reference    ┌──────────────┐
│ uploaded     ├──────────────►│ anchors.js   │  (shared 20-mers → position map)
│ sequences    │               └──────┬───────┘
└──────────────┘                      │ features in reference coordinates
                                      ▼
                              ┌──────────────┐   spans in reference
                              │ blocks.js    │   coordinates          ┌──────────────┐
                              │ + the SVG    ├───────────────────────►│ span tiling  │
                              └──────────────┘   --tile-region        └──────────────┘
```

### Unit 1: span tiling (Python)

**What it does.** Runs the existing tiling search inside a given span of the reference instead of
over the whole genome, for one or more spans, and merges the results into one scheme.

**Why in the engine.** Overlapping amplicons that actually overlap, two pools that alternate, and
a search that maximises the weakest amplicon's coverage are all properties of the tiling search.
The regions mode cannot give them: it returns one best pair per region, chosen independently, so
neighbouring pairs need not overlap and get no pool.

**Command line.**

```
vhirprimerlab -i in.fasta --tile-region F:3012-4631 --tile-region L:7074-13090 \
              --amplicon-len 400 [--overlap 50] -o out
```

- `--tile-region NAME:START-END` is repeatable, in 1-based reference coordinates, like `--region`.
- `--amplicon-len` is the target amplicon length, default 400 bp, with the existing ±20/35/50 %
  tolerance. The number of amplicons per span follows from it:
  `n = max(1, ceil((span_len - overlap) / (amplicon_len - overlap)))`.
- `--max-start-loss` and `--max-end-loss` (from PR #11) apply to the span's edges, so the span is
  covered end to end but for that allowance.
- `--tile-region` is mutually exclusive with `-n`, `--region` and `--regions-file`.

**Interface inside the engine.**

- `tiling.amplicons_for(span_len, amplicon_len, overlap) -> int`, a pure function.
- `tiling.tile(..., span: tuple[int, int] | None)`, where `span` is a consensus interval,
  0-based inclusive. With a span, candidates are restricted to it, the geometry and the end
  rules treat the span as the genome, and `end_loss` is measured from the span's edges.
  Amplicon coordinates stay absolute, so nothing downstream changes.
- A new `tiling.tile_spans(...)` runs `tile` per span and merges: amplicons numbered across the
  whole scheme, pools alternating, and the dimer check run once over the merged set, because in
  a multiplex every primer of a pool shares a tube whatever span it came from.
- Primer names come from the span: `F_LEFT`/`F_RIGHT` for a single amplicon, `L_1_LEFT`,
  `L_2_LEFT`, … when a span needs several. The `_LEFT`/`_RIGHT` suffixes stay, so `scheme.bed`
  and the amplicon-sequencing tools keep working.

**One fix this uncovers.** `geometry` refuses a genome whose amplicons could not overlap their
neighbours, which for `n == 1` is a condition about a neighbour that does not exist. Over a whole
genome it never fires; over a short span it does, and it would refuse a 150 bp feature that one
amplicon covers perfectly well. With `n == 1` that condition no longer applies, so a
single-amplicon span is limited only by the primers themselves. Without the fix the smallest
usable span is about 250 bp at the default settings.

**Reused unchanged:** candidate enumeration, conserved-primer picking, scoring against every
sequence, degenerate bases, the dimer swap, and every output file.

### Unit 2: annotation (browser)

**What it does.** Turns an accession into a list of features.

- Fetches `efetch.fcgi?db=nuccore&id=<accession>&rettype=gbc&retmode=xml` and parses it with
  `DOMParser`: organism, record length, and for each feature its kind, location, gene name and
  product.
- Locations are parsed to a list of intervals plus a strand, so `join(...)` and `complement(...)`
  are handled; a spliced feature keeps its pieces.
- Kept kinds: `gene`, `CDS`, `mat_peptide`, `5'UTR`, `3'UTR`, `stem_loop`, `misc_feature`.
  `source` is dropped.
- Cached in `localStorage` under the accession, so the same organism works offline afterwards.
  Every read and write is wrapped in try/catch; the cache is a convenience, never required.
- Requests are serialised with a small delay and retried with backoff, to stay inside NCBI's
  limits for clients without an API key.
- The same parser accepts a GenBank flat file the user loads by hand, so the mode needs no
  network at all.

**Refuses:** a record longer than 500 kb, with a message saying to pick the gene's record
instead of a chromosome.

### Unit 3: the anchor map (browser)

**What it does.** Converts positions on the record to positions on the run's reference.

- Indexes the unique 20-mers of the reference, walks the record every 10 bases, keeps the hits,
  and takes the longest increasing chain of them.
- Between anchors it interpolates; beyond the first and last anchor it extrapolates by offset.
- Reports its own quality: the number of anchors and the span of the record they cover. Below a
  threshold (fewer than 10 anchors, or less than half the record covered) the UI does not pretend
  to map: it draws the diagram in record coordinates and says so.
- The part of the record outside the anchor span is where the user's sequences have nothing. The
  features there are drawn greyed and cannot be selected.
- Measured accuracy: 13 nt worst case on hMPV. Spans sent to the engine are therefore padded by
  50 nt on each side, which is a search window, not a primer position, so the padding costs
  nothing and absorbs the error.

**A pure module.** No DOM, no network, so it is unit-tested with `node --test`.

### Unit 4: blocks and the diagram (browser)

**The diagram.** One SVG, built with the same helpers and CSS tokens as the existing amplicon
map: a position ruler in reference coordinates, then lanes for genes, for proteins (CDS), and for
mature peptides and other regions. Features carry their name, product, coordinates and length in
a tooltip and in `aria-label`, are reachable with the keyboard, and highlight in step with the
list beside them, exactly as the amplicon map does today.

**Blocks.** A block is a span the user has chosen, with a name.

- Clicking a feature creates the block that covers it, padded on each side by 50 nt for the
  anchor map's error plus the end-loss allowance, so the feature itself always ends up inside the
  inserts rather than under a primer.
- Handles at both ends adjust it; dragging snaps to feature edges and to 10 bp.
- Two controls apply to every block: the target amplicon length and the minimum overlap. Changing
  either recomputes how many amplicons each block needs, shown live.
- Overlapping blocks merge, with a note saying so.
- A list beside the diagram shows name, interval, length and amplicon count, with a remove button.
- A "whole genome" shortcut produces today's tiling.

**Pure logic in `blocks.js`:** creating a block from a feature, merging, snapping, and the
amplicon count. Tested with `node --test`. Drawing and dragging stay in `app.js` and are checked
in a real browser.

**Submitting.** The mode sends `--tile-region` per block, in reference coordinates, plus
`--amplicon-len` and `--overlap`. Nothing else about the request changes, so progress, results,
the map, the table and every download work as they do now.

## Data flow

1. The user drops a FASTA. The app picks the reference (first sequence, or their choice) as today.
2. They pick an organism: catalogue entry, name search, or a GenBank file they load.
3. The browser fetches the annotation (the accession is all that leaves the machine) and caches it.
4. The browser builds the anchor map from the reference and the record, and converts every
   feature to reference coordinates.
5. The diagram is drawn, cropped to the mapped part.
6. The user clicks features and adjusts blocks.
7. Pressing design submits `--tile-region` spans; the engine runs the span tiling and writes the
   usual `primers.tsv`, `amplicons.tsv`, `scheme.bed`, `report.txt` and `alignment.fasta`.

Step 3 is the only one that needs the network, and only the first time for a given accession.

## Error handling

| Situation | What happens |
|---|---|
| No network, or NCBI is down | The fetch fails with a clear message; the catalogue, the cache and the GenBank upload still work, and the regions mode is unaffected. |
| Wrong organism picked | Few anchors, or anchors covering little of the record. The UI says the record does not look like these sequences, shows the numbers, and lets the user continue or change. |
| Record over 500 kb | Refused, with a message to pick the gene's record. |
| Feature too short for any amplicon | The block is padded symmetrically to the smallest span the settings allow, and the UI says it was widened and by how much. |
| A span has no acceptable primers | The engine's existing message, naming the primer-free stretch, surfaces as today. |
| NCBI rate limit | Serialised requests, backoff, and the cache; the message names the wait. |
| Annotation with no CDS or gene features | The diagram shows what there is; if there is nothing to click, the UI says so and offers manual handles over the whole record. |

## Testing

**Engine, test first:**

- `amplicons_for` on the boundaries: a span shorter than one amplicon, exactly one, one base over.
- `geometry` with `n == 1` accepts a span that the neighbour-overlap condition used to refuse,
  and still refuses a span too short to hold two primers.
- `tile` with a span: candidates outside it are never used; the end rules apply to the span's
  edges, not the genome's; amplicon coordinates stay absolute.
- `tile_spans`: numbering across spans, pools alternating, overlaps inside each span, dimers
  checked across spans, names (`F_LEFT`, `L_1_LEFT`).
- CLI: `--tile-region` parsing, mutual exclusions, bad coordinates, a span outside the reference.
- Both an end-to-end run on the synthetic genome family and one on a span that needs a single
  amplicon, so the single-amplicon path is covered.

**Annotation and anchors:**

- Saved INSDSeq XML fixtures (hMPV, one influenza segment, one parasite gene record) so no test
  touches the network.
- Location parsing: plain, `complement`, `join`, and a malformed location.
- The anchor map against a synthetic reference built from a record with known insertions and
  deletions: features land where expected, the quality report fires when the record is unrelated,
  and the unmapped part is reported.

**Server:** the new mode's settings pass through `prepare.py` and `params.py` with the existing
validation tests extended; no new network code to test.

**Browser:** `node --test` over `anchors.js` and `blocks.js`, which have no DOM. The diagram,
dragging and the whole flow are verified in a real browser against the user's hMPV data, the way
the amplicon map fix was.

## Order of work

1. Span tiling in the engine, with the command-line option and its tests.
2. Annotation fetch and parse, plus the anchor map, both as pure modules with tests.
3. The diagram and block picking in the web app.
4. The catalogue and the name search.

Each step leaves the tool working. Step 1 alone is already useful on the command line.

## Phase 2, not in this spec

Automatic organism detection: a button that sends about 1 kb of the consensus to a BLAST service,
with a notice, a per-user contact address kept in the browser, and the catalogue as the fallback
when the service is unavailable. It fills in the accession that phase 1 asks the user for, so
nothing in phase 1 depends on it.
