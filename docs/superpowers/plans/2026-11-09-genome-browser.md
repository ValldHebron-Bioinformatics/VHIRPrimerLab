# Genome browser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a user pick the organism, see its genes drawn over their own sequences, click a protein, and get the primers that cover it.

**Architecture:** The engine gains span tiling: the existing tiling search run inside given spans of the reference, merged into one scheme. Everything about RefSeq stays in the browser — three dependency-free ES modules fetch the annotation, map record coordinates onto the run's reference with shared k-mers, and turn clicked features into spans — so the server makes no network calls and the mode also works in the GitHub Pages build.

**Tech Stack:** Python 3.10+ (biopython, primer3-py), FastAPI for the server, plain ES modules and SVG in the browser, pytest and `node --test`.

**Spec:** `docs/superpowers/specs/2026-11-09-genome-browser-design.md`

## Global Constraints

- No new Python dependencies: biopython and primer3-py only.
- No new JavaScript runtime dependencies; browser code is plain ES modules; JS tests run with `node --test` and nothing else.
- The server makes no network calls. Every NCBI request happens in the browser.
- Annotation records over **500 kb** are refused.
- Anchor map: k = **20**, step = **10**; the map is trusted only with **at least 10 anchors** covering **at least 50 %** of the record.
- A block is padded on each side by **50 nt** plus the end-loss allowance.
- `--amplicon-len` default **400** bp, with the existing ±20/35/50 % tolerance.
- Primer names keep the `_LEFT` / `_RIGHT` suffixes, so `scheme.bed` stays ARTIC-compatible.
- Every `localStorage` read and write is wrapped in try/catch; the cache is never required.
- One molecule per run: no segmented-genome handling.

## Review Focus

1. **A span outside the reference, or with start ≥ end** — a clear error naming the span, never a crash or a silent empty result. Test in Task 5.
2. **Two blocks that touch exactly** (`a.end + 1 == b.start`) — merged into one, so no amplicon is designed twice and the pools keep alternating. Test in Task 9.
3. **Features on the complement strand** — a `complement(...)` or `join(complement(...))` location still yields intervals in increasing coordinates, with the strand kept separately. Test in Task 8.
4. **A record much shorter than the sequences** (a gene record against whole-genome uploads, and the reverse) — anchors cover only part; the mapped range is reported and the diagram crops instead of stretching. Test in Task 7.
5. **`localStorage` unavailable or holding a corrupt entry** — the annotation is fetched afresh and the mode works; nothing throws. Test in Task 8.

---

### Task 1: amplicon count, and geometry for a single amplicon

**Files:**
- Modify: `vhirprimerlab/tiling.py` (`geometry`, new `amplicons_for`)
- Test: `tests/test_tiling_geometry.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `amplicons_for(span_len: int, amplicon_len: int, overlap: int) -> int`; `geometry` unchanged in signature, but accepts spans it used to refuse when `n == 1`.

- [ ] **Step 1: Write the failing tests**

```python
def test_amplicons_for() -> None:
    assert amplicons_for(400, 400, 50) == 1
    assert amplicons_for(150, 400, 50) == 1          # shorter than one amplicon
    assert amplicons_for(750, 400, 50) == 2          # (750 - 50) / 350 = 2
    assert amplicons_for(751, 400, 50) == 3
    assert amplicons_for(6018, 400, 50) == 18        # the hMPV L gene


def test_geometry_allows_one_short_amplicon() -> None:
    # 150 bp holds two 25 bp primers and an insert; with n == 1 there is no neighbour,
    # so the overlap condition must not apply
    g = geometry(150, 1, 50, 0.20, 25)
    assert g is not None and g.len_min == 120
    assert geometry(40, 1, 50, 0.20, 25) is None     # too short for two primers


def test_geometry_still_refuses_short_amplicons_with_neighbours() -> None:
    assert geometry(300, 2, 50, 0.20, 25) is None
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_tiling_geometry.py -v`
Expected: FAIL, `amplicons_for` is not defined and `geometry(150, 1, …)` returns None.

- [ ] **Step 3: Implement both in `vhirprimerlab/tiling.py`**

`amplicons_for(span_len, amplicon_len, overlap) -> int` is `max(1, ceil((span_len - overlap) / (amplicon_len - overlap)))`, guarding `amplicon_len <= overlap` by raising `PipelineError("--amplicon-len must be larger than --overlap")`.

In `geometry`, the `max_outer < overlap + 2 * primer_len_max` rejection applies only when `n > 1`. With `n == 1` the requirement is instead that the amplicon can hold two primers and at least one base of insert: `len_min >= 2 * primer_len_max + 1`.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS, 280 tests plus the new ones.

- [ ] **Step 5: Commit**

```bash
git add vhirprimerlab/tiling.py tests/test_tiling_geometry.py
git commit -m "feat(tiling): amplicons_for, and allow a single short amplicon"
```

---

### Task 2: tile one span

**Files:**
- Modify: `vhirprimerlab/candidates.py` (`end_loss`, `enumerate_primers`, `score_primers`)
- Modify: `vhirprimerlab/tiling.py` (`choose`, `_no_fit`, `tile`)
- Test: `tests/test_tiling.py`, `tests/test_candidates.py`

**Interfaces:**
- Consumes: `amplicons_for` from Task 1.
- Produces:
  - `end_loss(c: PrimerCandidate, span: tuple[int, int]) -> int` — bases of the span left outside the insert at the end this primer reaches. Replaces today's `end_loss(c, genome_len)`; the whole genome is the span `(0, genome_len - 1)`.
  - `enumerate_primers(cons_seq: str, params: TilingParams, span: tuple[int, int] | None = None)` and `score_primers(..., span: tuple[int, int] | None = None)`, where a span restricts candidates to it and applies the end rules at its edges.
  - `tile(aln, cons, stats, cmap, params, progress=None, span: tuple[int, int] | None = None) -> TilingResult`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_tiling.py
def test_tile_within_a_span_stays_inside_it(family) -> None:
    aln, cons, stats, cmap = family
    result = tile(aln, cons, stats, cmap, _params(n_amplicons=2), span=(1000, 2999))
    for a in result.amplicons:
        assert a.fwd.start >= 1000 and a.rev.end <= 2999
    assert len(result.amplicons) == 2


def test_the_span_edges_carry_the_end_rules(family) -> None:
    aln, cons, stats, cmap = family
    result = tile(aln, cons, stats, cmap, _params(n_amplicons=2), span=(1000, 2999))
    assert result.amplicons[0].fwd.end + 1 - 1000 <= 50
    assert 2999 + 1 - result.amplicons[-1].rev.start <= 50


def test_a_span_with_no_primers_names_the_span_edge(family) -> None:
    aln, cons, stats, cmap = family
    with pytest.raises(PipelineError, match="of the span start"):
        tile(aln, cons, stats, cmap, _params(n_amplicons=1, max_start_loss=0), span=(1000, 2999))
```

```python
# tests/test_candidates.py
def test_end_rules_follow_the_span(family) -> None:
    aln, stats, cons = family
    fwds, revs = enumerate_primers(cons.seq, PARAMS, span=(1000, 2999))
    assert all(1000 <= c.start and c.end <= 2999 for c in fwds + revs)
    # inside the span's first 50 bases the end rules apply: longer primers, any GC
    near = [c for c in fwds if end_loss(c, (1000, 2999)) <= 50]
    assert near and any(c.length > 25 or c.gc < 40 for c in near)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_tiling.py tests/test_candidates.py -q`
Expected: FAIL, `tile()` takes no `span` argument.

- [ ] **Step 3: Implement the span in `candidates.py` and `tiling.py`**

`end_loss` takes the span and measures from its edges. `enumerate_primers` walks only the span's anchor positions and compares each against the span's edges when deciding whether the end rules apply. `score_primers` passes the span on for the GC relaxation.

In `tiling.py`: `tile` defaults `span` to `(0, len(cons.seq) - 1)`, derives `genome_len` from the span, passes the span to `enumerate_primers`, `score_primers` and `choose`, and computes `uncovered_ends` from the span's edges. `choose` takes `span` instead of `genome_len`. `_no_fit`'s message says "of the span start/end" when the span is not the whole consensus, and keeps today's wording otherwise.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vhirprimerlab/candidates.py vhirprimerlab/tiling.py tests/test_tiling.py tests/test_candidates.py
git commit -m "feat(tiling): tile inside a span of the reference"
```

---

### Task 3: amplicon names come from the model

**Files:**
- Modify: `vhirprimerlab/models.py` (`Amplicon`), `vhirprimerlab/tiling.py` (`_primer_names`, `_dimer_pairs`, `tile`), `vhirprimerlab/io.py` (`_names` and its four call sites)
- Test: `tests/test_tiling_io.py`

**Interfaces:**
- Consumes: Task 2's `tile`.
- Produces: `Amplicon.name: str` — the pair name without a side suffix, `VHIRPrimerLab_3` for whole-genome tiling. `io._names(name: str) -> tuple[str, str]` and `tiling._primer_names(name: str) -> tuple[str, str]` both return `(f"{name}_LEFT", f"{name}_RIGHT")`.

- [ ] **Step 1: Write the failing tests**

```python
def test_whole_genome_names_are_unchanged(result) -> None:
    assert [a.name for a in result.amplicons[:2]] == ["VHIRPrimerLab_1", "VHIRPrimerLab_2"]
    assert tiling_order(result)[:2] == [
        ("VHIRPrimerLab_1_LEFT", result.amplicons[0].fwd.seq),
        ("VHIRPrimerLab_1_RIGHT", result.amplicons[0].rev.seq),
    ]


def test_the_order_sheet_follows_the_amplicon_name(result) -> None:
    result.amplicons[0].name = "F"
    assert tiling_order(result)[0][0] == "F_LEFT"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_tiling_io.py -q`
Expected: FAIL, `Amplicon` has no `name`.

- [ ] **Step 3: Add the field and route every name through it**

`Amplicon.name` is set in `tile` to `f"VHIRPrimerLab_{index}"`. `io._names` and `tiling._primer_names` take the name. The four `io.py` call sites (`tiling_order`, `write_tiling_tsv`, `write_scheme_bed`, the report section) and `_dimer_pairs` use `a.name`.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS, with every existing output file byte-identical.

- [ ] **Step 5: Commit**

```bash
git add vhirprimerlab/models.py vhirprimerlab/tiling.py vhirprimerlab/io.py tests/test_tiling_io.py
git commit -m "refactor(tiling): primer names come from the amplicon"
```

---

### Task 4: tile several spans into one scheme

**Files:**
- Modify: `vhirprimerlab/tiling.py` (new `tile_spans`), `vhirprimerlab/models.py` (`TilingResult.spanned`, `TilingParams.amplicon_len`), `vhirprimerlab/io.py` (the report's uncovered-ends line)
- Test: `tests/test_tiling.py`

**Interfaces:**
- Consumes: `tile` (Task 2), `Amplicon.name` (Task 3), `amplicons_for` (Task 1).
- Produces: `TilingParams.amplicon_len: int = 400`, and
  `tile_spans(aln, cons, stats, cmap, params, spans: list[Region], progress=None) -> TilingResult`
  with `TilingResult.spanned = True`. `Region` is the existing `vhirprimerlab.models.Region`
  (name, start, end in reference coordinates).

- [ ] **Step 1: Write the failing tests**

```python
def test_tile_spans_numbers_and_names_across_spans(family) -> None:
    aln, cons, stats, cmap = family
    spans = [Region("F", 1001, 3000), Region("SH", 4001, 4600)]
    res = tile_spans(aln, cons, stats, cmap, _params(amplicon_len=700), spans)
    assert [a.name for a in res.amplicons] == ["F_1", "F_2", "F_3", "SH"]
    assert [a.index for a in res.amplicons] == [1, 2, 3, 4]
    assert [a.pool for a in res.amplicons] == [1, 2, 1, 2]
    assert res.spanned is True


def test_tile_spans_overlaps_inside_a_span_only(family) -> None:
    aln, cons, stats, cmap = family
    spans = [Region("F", 1001, 3000), Region("SH", 4001, 4600)]
    res = tile_spans(aln, cons, stats, cmap, _params(amplicon_len=700), spans)
    assert res.amplicons[0].overlap_next >= 50          # inside F
    assert res.amplicons[2].overlap_next is None        # last of F: no neighbour
    assert res.amplicons[-1].overlap_next is None


def test_tile_spans_checks_dimers_across_spans(family, monkeypatch) -> None:
    aln, cons, stats, cmap = family
    monkeypatch.setattr(tiling, "dimer_dg", lambda a, b: -20_000)
    spans = [Region("F", 1001, 3000), Region("SH", 4001, 4600)]
    res = tile_spans(aln, cons, stats, cmap, _params(amplicon_len=700), spans)
    assert any("F_1_LEFT" in w and "SH" in w for w in res.dimer_warnings)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_tiling.py -k spans -q`
Expected: FAIL, `tile_spans` is not defined.

- [ ] **Step 3: Implement `tile_spans`**

For each span in reference order: convert it with `cmap.region_to_included` to a consensus interval, take `n = amplicons_for(span_len, params.amplicon_len, params.overlap)`, call `tile` with that span and `n`, and name its amplicons after the span — the span's own name when it has one amplicon, `name_1`, `name_2`, … when it has several. Then renumber `index` across the scheme, set `pool` by `1 if index % 2 else 2`, set `overlap_next` to `None` on each span's last amplicon, run the existing dimer check once over the merged list, and set `uncovered_ends` to the largest loss at any span start and at any span end. A `PipelineError` from one span is re-raised with `f"Region '{span.name}': "` in front.

`TilingParams.amplicon_len: int = 400` and `TilingResult.spanned: bool = False`; in `io.py` the report line reads `Not covered at the span edges: at most N bp (start), M bp (end)` when `spanned`, and keeps today's wording otherwise.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vhirprimerlab/tiling.py vhirprimerlab/models.py vhirprimerlab/io.py tests/test_tiling.py
git commit -m "feat(tiling): tile several spans into one scheme"
```

---

### Task 5: the command line

**Files:**
- Modify: `vhirprimerlab/cli.py` (`build_parser`, `TILING_CONFLICTS`, `check_tiling_params`, `parse_args`, `validate_args`, `run_tiling`)
- Modify: `README.md` (the tiling options table and a span-tiling example)
- Test: `tests/test_tiling_cli.py`, `tests/test_cli.py`

**Interfaces:**
- Consumes: `tile_spans` (Task 4).
- Produces: `--tile-region NAME:START-END` (repeatable), `--amplicon-len` (default 400), and `cli.check_span(region: Region, ref_length: int) -> None`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_cli.py
def test_tile_region_routes_to_tiling() -> None:
    args, runner, _ = parse_args(["-i", "in.fa", "--tile-region", "F:101-900", "-o", "out"])
    assert runner is run_tiling and args.tile_region == ["F:101-900"]
    assert args.amplicon_len == 400


def test_tile_region_cannot_be_combined() -> None:
    for clash in (["-n", "4"], ["--region", "x:1-100"], ["--regions-file", "r.tsv"]):
        with pytest.raises(SystemExit):
            parse_args(["-i", "in.fa", "--tile-region", "F:101-900", *clash, "-o", "out"])


def test_a_span_outside_the_reference_is_refused() -> None:
    with pytest.raises(PipelineError, match="outside the reference"):
        check_span(Region("F", 9000, 9500), ref_length=1000)


def test_a_backwards_span_is_refused() -> None:
    with pytest.raises(PipelineError, match="start must be smaller than end"):
        check_span(Region("F", 900, 101), ref_length=1000)
```

```python
# tests/test_tiling_cli.py
def test_span_tiling_end_to_end(tmp_path) -> None:
    out = tmp_path / "out"
    code = main(["-i", str(FASTA), "--tile-region", "F:1001-3000", "--tile-region", "SH:4001-4600",
                 "--amplicon-len", "700", "-o", str(out)])
    assert code == 0
    names = [r["forward_name"] for r in read_rows(out / "amplicons.tsv")]
    assert names[:2] == ["F_1_LEFT", "F_2_LEFT"] and names[-1] == "SH_LEFT"
    assert (out / "scheme.bed").read_text().count("SH_RIGHT") == 1
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_cli.py tests/test_tiling_cli.py -q`
Expected: FAIL, unrecognised argument `--tile-region`.

- [ ] **Step 3: Add the arguments and the route**

`--tile-region` (`action="append"`, metavar `[NAME:]START-END`) and `--amplicon-len` (int, default 400) join `build_parser`. `TILING_CONFLICTS` gains `--tile-region`, and `--tile-region` clashes with `-n`, `--region` and `--regions-file` through `parser.error`. In `parse_args`, `args.tile_region` routes to `run_tiling` the way `-n` does, and `args.max_degeneracy` defaults to 4 there too. `validate_args` sends it to `check_tiling_params`, which requires `--amplicon-len > --overlap` and accepts a missing `-n`.

`run_tiling` resolves the spans with the existing `resolve_regions(a.tile_region, None, inp.ref_len)`, validates each with `check_span`, and calls `tile_spans` instead of `tile` when spans were given.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vhirprimerlab/cli.py README.md tests/test_cli.py tests/test_tiling_cli.py
git commit -m "feat(cli): --tile-region to tile chosen spans of the reference"
```

---

### Task 6: the web request

**Files:**
- Modify: `vhirprimerlab/web/params.py` (`MODES`, `SETTINGS`), `vhirprimerlab/web/prepare.py` (`prepare_job`), `vhirprimerlab/web/result.py` (`build_result`)
- Test: `tests/test_web_params.py`, `tests/test_web_api.py`

**Interfaces:**
- Consumes: Task 5's flags.
- Produces: the mode `"genes"`, whose settings are `tile_regions` (a list of `NAME:START-END`, flag `--tile-region`, kind `"regions"`), `amplicon_len`, `overlap`, plus the shared design settings.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_web_params.py
def test_genes_mode_builds_tile_region_flags() -> None:
    argv = build_argv("genes", {"tile_regions": ["F:3012-4631", "L:7074-13090"],
                                "amplicon_len": 400})
    assert argv == ["--tile-region", "F:3012-4631", "--tile-region", "L:7074-13090",
                    "--amplicon-len", "400"]


def test_genes_mode_rejects_tiling_only_settings() -> None:
    with pytest.raises(ValueError, match="unknown setting"):
        build_argv("genes", {"amplicons": 4})


def test_genes_mode_needs_at_least_one_region(limits) -> None:
    with pytest.raises(ValueError, match="Pick at least one region"):
        prepare_job("genes", {"tile_regions": []}, FASTA_BYTES, "", limits, ["mafft"])
```

```python
# tests/test_web_api.py
def test_a_genes_job_reads_the_amplicon_table(tmp_path) -> None:
    res = build_result(_finished_job(tmp_path), ["-i", "input.fasta", "--tile-region", "F:1-900"], 0)
    assert res["mode"] == "tile"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest tests/test_web_params.py tests/test_web_api.py -q`
Expected: FAIL, unknown mode `genes`.

- [ ] **Step 3: Add the mode**

`MODES` gains `"genes"`; `SETTINGS["genes"]` holds `tile_regions`, `amplicon_len`, `overlap` and `**_DESIGN`. In `prepare_job`, the `"genes"` mode refuses an empty `tile_regions` with `"Pick at least one region on the genome first."` and refuses a span that ends past the longest uploaded sequence. In `build_result`, the mode is `"tile"` when `-n` **or** `--tile-region` is in `argv`.

`engine.js` and `engine-worker.js` need no change: they pass the mode through, and their `FILES` map is keyed by the *result's* mode, which is `"tile"` for a genes run.

- [ ] **Step 4: Run the whole suite**

Run: `python -m pytest -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add vhirprimerlab/web tests/test_web_params.py tests/test_web_api.py
git commit -m "feat(web): accept the genes mode and its spans"
```

---

### Task 7: the anchor map

**Files:**
- Create: `vhirprimerlab/web/static/anchors.js`
- Create: `tests/browser/anchors.test.mjs`, `tests/test_static_modules.py`
- Modify: `tests/test_build_site.py` (the list of files the built site must contain)

**Interfaces:**
- Consumes: nothing.
- Produces:
  ```js
  buildAnchorMap(record, reference, { k = 20, step = 10 } = {}) -> {
    anchors: [[recordPos, referencePos], …],   // 1-based, strictly increasing in both
    count, covered,                            // anchors, and the fraction of the record they span
    trusted,                                   // count >= 10 && covered >= 0.5
    mappedRange: [recordStart, recordEnd],     // record positions the anchors cover
    convert(pos),                              // record position -> reference position
  }
  ```

- [ ] **Step 1: Write the failing tests**

```js
// tests/browser/anchors.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { buildAnchorMap } from '../../vhirprimerlab/web/static/anchors.js';

// deterministic pseudo-random ACGT, so no 20-mer repeats
const seq = (n, seed = 1) => { /* implement in the test file */ };

test('an identical record maps one to one', () => {
  const s = seq(2000);
  const map = buildAnchorMap(s, s);
  assert.ok(map.trusted);
  assert.equal(map.convert(1), 1);
  assert.equal(map.convert(1500), 1500);
});

test('an insertion in the reference shifts later positions', () => {
  const s = seq(2000);
  const reference = s.slice(0, 1000) + seq(30, 9) + s.slice(1000);
  const map = buildAnchorMap(s, reference);
  assert.equal(map.convert(500), 500);
  assert.ok(Math.abs(map.convert(1500) - 1530) <= 15);
});

test('a record longer than the reference reports the mapped range only', () => {
  const s = seq(3000);
  const map = buildAnchorMap(s, s.slice(800, 1800));   // the reference is an inner gene
  assert.ok(map.mappedRange[0] >= 790 && map.mappedRange[1] <= 1810);
});

test('an unrelated record is not trusted', () => {
  assert.equal(buildAnchorMap(seq(2000, 1), seq(2000, 7)).trusted, false);
});
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/browser && node --test anchors.test.mjs`
Expected: FAIL, the module does not exist.

- [ ] **Step 3: Implement `anchors.js`**

Index the reference's k-mers, keeping only those that occur once. Walk the record every `step` bases, keep the hits, and take the longest strictly increasing chain of reference positions by patience sorting with back-pointers. `convert` interpolates linearly between the surrounding anchors and extrapolates by a constant offset outside them. Export `buildAnchorMap` as a named export; the module touches neither the DOM nor the network.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/browser && node --test anchors.test.mjs`
Expected: PASS.

- [ ] **Step 5: Run these tests from pytest too**

Create `tests/test_static_modules.py`: it runs `node --test anchors.test.mjs` in `tests/browser`
and asserts a zero exit code, skipping only when `shutil.which("node")` is None. Unlike
`tests/test_browser_engine.py` it must not require Pyodide or the built engine, because these
modules are plain JavaScript. Later tasks add their own `.mjs` files to its argument list.

Add `anchors.js` to the file list in `tests/test_build_site.py`, so the module is checked into
the published site build.

- [ ] **Step 6: Commit**

```bash
git add vhirprimerlab/web/static/anchors.js tests/browser/anchors.test.mjs \
        tests/test_static_modules.py tests/test_build_site.py
git commit -m "feat(web): map record coordinates onto the reference with shared k-mers"
```

---

### Task 8: the annotation

**Files:**
- Create: `vhirprimerlab/web/static/annotation.js`
- Create: `tests/browser/fixtures/NC_039199.xml` (hMPV, trimmed to its features), `tests/browser/fixtures/influenza-segment.xml`, `tests/browser/fixtures/parasite-gene.xml`, `tests/browser/fixtures/hmpv.gb` (a short GenBank flat file)
- Test: `tests/browser/annotation.test.mjs`

**Interfaces:**
- Consumes: nothing.
- Produces:
  ```js
  parseLocation(text) -> { strand: '+' | '-', intervals: [[start, end], …] }  // 1-based, increasing
  parseInsdseq(xmlText) -> { accession, organism, length, sequence, features, skipped }
  parseGenbank(text)    -> the same shape
  fetchAnnotation(accession, { fetch, storage } = {}) -> Promise<annotation>
  MAX_RECORD_LENGTH = 500_000
  ```
  A feature is `{ kind, name, product, strand, intervals, start, end }`, where `start`/`end` bound its intervals.

- [ ] **Step 1: Write the failing tests**

```js
test('a plain location', () => {
  assert.deepEqual(parseLocation('40..1224'), { strand: '+', intervals: [[40, 1224]] });
});

test('a complement location keeps increasing intervals', () => {
  assert.deepEqual(parseLocation('complement(40..1224)'),
    { strand: '-', intervals: [[40, 1224]] });
});

test('a joined complement location', () => {
  const loc = parseLocation('join(complement(20..30),complement(1..10))');
  assert.deepEqual(loc, { strand: '-', intervals: [[1, 10], [20, 30]] });
});

test('a malformed location throws', () => {
  assert.throws(() => parseLocation('nonsense'), /location/);
});

test('the hMPV record yields its genes and proteins', () => {
  const a = parseInsdseq(readFileSync('tests/browser/fixtures/NC_039199.xml', 'utf8'));
  assert.equal(a.organism, 'human metapneumovirus');
  assert.equal(a.length, 13350);
  const f = a.features.find((x) => x.kind === 'CDS' && x.name === 'F');
  assert.deepEqual([f.start, f.end], [3052, 4671]);
  assert.equal(f.product, 'fusion protein');
  assert.ok(!a.features.some((x) => x.kind === 'source'));
});

test('a record over the cap is refused', async () => {
  const fetch = async () => ({ ok: true, text: async () => bigRecordXml });
  await assert.rejects(fetchAnnotation('NC_000000.1', { fetch }), /500/);
});

test('a corrupt cache entry is ignored', async () => {
  const storage = { getItem: () => '{not json', setItem() {} };
  const fetch = async () => ({ ok: true, text: async () => hmpvXml });
  assert.equal((await fetchAnnotation('NC_039199.1', { fetch, storage })).organism,
               'human metapneumovirus');
});

test('storage that throws does not break the fetch', async () => {
  const storage = { getItem() { throw new Error('private'); },
                    setItem() { throw new Error('private'); } };
  const fetch = async () => ({ ok: true, text: async () => hmpvXml });
  assert.ok(await fetchAnnotation('NC_039199.1', { fetch, storage }));
});

test('a GenBank file parses to the same shape', () => {
  const a = parseGenbank(readFileSync('tests/browser/fixtures/hmpv.gb', 'utf8'));
  assert.equal(a.features.find((x) => x.name === 'F').start, 3052);
});
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/browser && node --test annotation.test.mjs`
Expected: FAIL, the module does not exist.

- [ ] **Step 3: Implement `annotation.js`**

`parseInsdseq` uses `DOMParser` when one is available and falls back to a small regex reader under Node, so the fixtures parse in tests without a DOM library. Kept kinds: `gene`, `CDS`, `mat_peptide`, `5'UTR`, `3'UTR`, `stem_loop`, `misc_feature`. `fetchAnnotation` takes its `fetch` and `storage` by injection, defaulting to `globalThis.fetch` and `localStorage`; it reads the cache inside try/catch, refuses a record over `MAX_RECORD_LENGTH` with a message naming the cap and advising a gene record, retries twice with 1 s and 3 s backoff on a 429 or a 5xx, and writes the cache inside try/catch.

The URL is `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbc&retmode=xml&id=<accession>`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/browser && node --test annotation.test.mjs`
Expected: PASS.

- [ ] **Step 5: Add the module to both test lists and commit**

Add `annotation.test.mjs` to `tests/test_static_modules.py` and `annotation.js` to
`tests/test_build_site.py`.

```bash
git add vhirprimerlab/web/static/annotation.js tests/browser/annotation.test.mjs \
        tests/browser/fixtures tests/test_static_modules.py tests/test_build_site.py
git commit -m "feat(web): read a RefSeq or GenBank annotation"
```

---

### Task 9: blocks

**Files:**
- Create: `vhirprimerlab/web/static/blocks.js`
- Test: `tests/browser/blocks.test.mjs`

**Interfaces:**
- Consumes: a feature from Task 8, a map from Task 7.
- Produces:
  ```js
  PAD = 50
  blockFromFeature(feature, { convert, pad = PAD, endLoss = 50, limit }) -> { name, start, end }
  mergeBlocks(blocks) -> blocks                 // sorted, touching or overlapping ones merged
  ampliconCount(block, ampliconLen, overlap) -> number
  snapPosition(pos, edges, grid = 10) -> number
  smallestSpan(ampliconLen, overlap, primerLenMax, tolerance) -> number
  ```

- [ ] **Step 1: Write the failing tests**

```js
test('a block covers its feature with padding on both sides', () => {
  const b = blockFromFeature({ name: 'F', start: 3052, end: 4671 },
                             { convert: (p) => p - 40, limit: 13232 });
  assert.equal(b.name, 'F');
  assert.ok(b.start <= 3012 - 99 && b.end >= 4631 + 99);
});

test('padding never leaves the reference', () => {
  const b = blockFromFeature({ name: 'N', start: 1, end: 200 },
                             { convert: (p) => p, limit: 500 });
  assert.equal(b.start, 1);
  assert.equal(b.end, 500);
});

test('blocks that touch exactly are merged', () => {
  const merged = mergeBlocks([{ name: 'A', start: 1, end: 100 },
                              { name: 'B', start: 101, end: 200 }]);
  assert.equal(merged.length, 1);
  assert.deepEqual([merged[0].start, merged[0].end], [1, 200]);
  assert.equal(merged[0].name, 'A+B');
});

test('separate blocks stay separate and sorted', () => {
  const merged = mergeBlocks([{ name: 'B', start: 500, end: 600 },
                              { name: 'A', start: 1, end: 100 }]);
  assert.deepEqual(merged.map((b) => b.name), ['A', 'B']);
});

test('the amplicon count matches the engine', () => {
  assert.equal(ampliconCount({ start: 1, end: 400 }, 400, 50), 1);
  assert.equal(ampliconCount({ start: 1, end: 750 }, 400, 50), 2);
  assert.equal(ampliconCount({ start: 1, end: 751 }, 400, 50), 3);
});

test('dragging snaps to feature edges first, then to the grid', () => {
  assert.equal(snapPosition(3055, [3052, 4671]), 3052);
  assert.equal(snapPosition(3104, [3052, 4671]), 3100);
});
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd tests/browser && node --test blocks.test.mjs`
Expected: FAIL, the module does not exist.

- [ ] **Step 3: Implement `blocks.js`**

`blockFromFeature` converts the feature's bounds with `convert`, pads by `pad + endLoss` on each side, clamps to `[1, limit]`, and widens symmetrically to `smallestSpan` when shorter. `mergeBlocks` sorts by start and merges while `next.start <= current.end + 1`, joining names with `+`. `ampliconCount` is the same formula as the engine's `amplicons_for` — deliberately duplicated, with a comment naming `vhirprimerlab/tiling.py` as the other copy and the test above as the contract. `snapPosition` snaps to any edge within 10 bp, else rounds to the grid.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd tests/browser && node --test blocks.test.mjs`
Expected: PASS.

- [ ] **Step 5: Add the module to both test lists and commit**

Add `blocks.test.mjs` to `tests/test_static_modules.py` and `blocks.js` to
`tests/test_build_site.py`.

```bash
git add vhirprimerlab/web/static/blocks.js tests/browser/blocks.test.mjs \
        tests/test_static_modules.py tests/test_build_site.py
git commit -m "feat(web): turn clicked features into amplicon blocks"
```

---

### Task 10: the genome page

**Files:**
- Modify: `vhirprimerlab/web/static/index.html` (a nav entry, the organism fieldset, the diagram and the block list)
- Modify: `vhirprimerlab/web/static/app.js` (`mode`, `applyMode`, `settings`, `setPage`, and the new drawing and dragging)
- Modify: `vhirprimerlab/web/static/app.css`
- Test: verified in a real browser; the pure logic is already covered by Tasks 7–9

**Interfaces:**
- Consumes: `buildAnchorMap`, `fetchAnnotation`, `parseGenbank`, `blockFromFeature`, `mergeBlocks`, `ampliconCount`, `snapPosition`.
- Produces: `settings('genes')` returning `{ tile_regions: ['F:3012-4631', …], amplicon_len, overlap, … }`.

- [ ] **Step 1: Add the page shell**

A third `.modes` link, `data-mode="genes"`, labelled "Pick by gene". `mode()` returns `'genes'` when `state.page === 'genes'`. The page keeps the existing sequences fieldset and adds: an organism fieldset (a `<select>` for the catalogue, a search box, a "Load a GenBank file" button and an error paragraph), an `<svg id="genome">` inside a figure, the two block controls (target amplicon length, minimum overlap) and a `<ul id="blocks">`.

- [ ] **Step 2: Draw the genome**

`drawGenome(annotation, map, blocks)` builds the SVG with the same `el()` helper, ruler and CSS tokens as `drawMap`: a position ruler in reference coordinates, then one lane for `gene`, one for `CDS`, one for everything else. Each feature is a `<g tabindex="0">` with `<title>` and `aria-label` giving name, product, span and length, and features outside `map.mappedRange` get the class `off` and no click handler. When `map.trusted` is false the ruler is labelled with record coordinates and a note says so.

- [ ] **Step 3: Pick and adjust blocks**

Clicking a feature calls `blockFromFeature`, then `mergeBlocks`, then redraws. Each block is drawn as a bracket above the lanes with a handle at each end; dragging a handle moves it through `snapPosition` and redraws live. The list shows name, span, length and `ampliconCount`, with a remove button, and hovering a row highlights its block, reusing today's `select()` pattern. Changing either control recomputes every count.

- [ ] **Step 4: Submit**

`settings('genes')` returns the merged blocks as `tile_regions` strings, plus `amplicon_len` and `overlap`. `submit` sends mode `genes`. Nothing else in the run path changes.

- [ ] **Step 5: Verify in a browser**

Serve `python -m vhirprimerlab web --port 8791 --no-browser`, load `hmpv_aligned.fasta`, pick human metapneumovirus, click the F gene and then the L gene, check the diagram, the counts and the merge note, run the design, and confirm `primers.tsv` names start with `F_` and `L_`.

- [ ] **Step 6: Commit**

```bash
git add vhirprimerlab/web/static
git commit -m "feat(web): pick amplicons by clicking genes on the genome"
```

---

### Task 11: the catalogue and the name search

**Files:**
- Create: `vhirprimerlab/web/static/organisms.json`
- Modify: `vhirprimerlab/web/static/annotation.js` (`searchOrganisms`), `vhirprimerlab/web/static/app.js` (wiring)
- Modify: `README.md` (a "Pick by gene" section)
- Test: `tests/browser/annotation.test.mjs`, plus a browser check

**Interfaces:**
- Consumes: `fetchAnnotation`.
- Produces: `searchOrganisms(query, { fetch } = {}) -> Promise<[{ accession, title, length }]>`; `organisms.json` as `[{ "label": "Human metapneumovirus", "accession": "NC_039199.1", "note": "13.4 kb, 8 genes" }]`.

- [ ] **Step 1: Write the failing test**

```js
test('a name search returns accessions with their titles', async () => {
  const fetch = async (url) => ({ ok: true, text: async () =>
    url.includes('esearch') ? esearchXml : esummaryXml });
  const hits = await searchOrganisms('human metapneumovirus', { fetch });
  assert.equal(hits[0].accession, 'NC_039199.1');
  assert.match(hits[0].title, /metapneumovirus/i);
});

test('a search with no hits returns an empty list', async () => {
  const fetch = async () => ({ ok: true, text: async () => emptyEsearchXml });
  assert.deepEqual(await searchOrganisms('zzzz', { fetch }), []);
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cd tests/browser && node --test annotation.test.mjs`
Expected: FAIL, `searchOrganisms` is not exported.

- [ ] **Step 3: Implement the search and the catalogue**

`searchOrganisms` calls `esearch.fcgi?db=nuccore&term=<query> AND srcdb_refseq[prop]&retmax=10`, then `esummary.fcgi` for the ids, and returns accession, title and length, leaving out records over the cap. `organisms.json` lists the usual targets: human metapneumovirus, RSV-A, RSV-B, SARS-CoV-2 and influenza A segments 4 and 6. The select fills from it; picking an entry calls `fetchAnnotation`.

- [ ] **Step 4: Run every test**

Run: `python -m pytest -q`, which now includes the Node module tests through
`tests/test_static_modules.py`. Add `organisms.json` to `tests/test_build_site.py`.
Expected: PASS.

- [ ] **Step 5: Verify in a browser and commit**

Search for an organism that is not in the catalogue, pick it, and confirm the diagram draws.

```bash
git add vhirprimerlab/web/static README.md tests/browser/annotation.test.mjs
git commit -m "feat(web): an organism catalogue and a name search"
```
