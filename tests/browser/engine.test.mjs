// The browser engine, end to end, in Node (same Pyodide and WebAssembly as browsers).
//   cd tests/browser && npm install && PYTHON=<python with vhirprimerlab> node --test
import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { createRequire } from "node:module";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "../..");
const PYTHON = process.env.PYTHON || "python3";
const require = createRequire(import.meta.url);

const site = fs.mkdtempSync(path.join(os.tmpdir(), "vhirprimerlab-site-"));
const built = spawnSync(PYTHON, [path.join(ROOT, "tools/build_site.py"), path.join(site, "_site")], { encoding: "utf8" });
assert.equal(built.status, 0, built.stderr);
const APP = path.join(site, "_site", "app");

const { loadPyodide } = require(path.join(HERE, "node_modules/pyodide/pyodide.js"));
const { createEngineCore } = require(path.join(APP, "engine-worker.js"));
const mafft = require(path.join(APP, "mafft-runner.js"));
const mafftFactories = Object.fromEntries(mafft.MAFFT_PROGRAMS.map(
  (name) => [name, require(path.join(APP, `engine/mafft-${name}.js`))]));
const core = await createEngineCore({
  loadPyodide: () => loadPyodide(),
  fetchBytes: async (rel) => new Uint8Array(fs.readFileSync(path.join(APP, rel))),
  mafftFactories,
  mafft,
});
const genomes = fs.readFileSync(path.join(ROOT, "examples/genomes.fasta"), "utf8");

test("tiles the example genomes like the command line", () => {
  const stages = [];
  const out = core.run({ mode: "tile", settings: { amplicons: 8 }, fasta: genomes, primers: "" },
    (stage) => stages.push(stage));
  assert.equal(out.error, undefined, out.error);
  const { result, files, zip } = out;
  assert.equal(result.exit_code, 0);
  assert.equal(result.mode, "tile");
  assert.equal(result.rows.length, 8);
  assert.equal(result.reference_length, 9000);
  assert.ok(stages.includes("aligning") && stages.includes("done"));
  assert.ok(files["primers.tsv"] && files["scheme.bed"] && files["report.txt"] && files["alignment.fasta"]);
  assert.equal(String.fromCharCode(zip[0], zip[1]), "PK");

  // the command line on the alignment the browser produced gives the same primers
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "vhirprimerlab-native-"));
  fs.writeFileSync(path.join(dir, "aln.fasta"), files["alignment.fasta"]);
  const native = spawnSync(PYTHON, ["-m", "vhirprimerlab", "--aligned", path.join(dir, "aln.fasta"), "-n", "8",
    "-o", path.join(dir, "out")], { encoding: "utf8", cwd: ROOT });
  assert.equal(native.status, 0, native.stderr);
  assert.equal(files["primers.tsv"], fs.readFileSync(path.join(dir, "out/primers.tsv"), "utf8"));
});

test("checks primer pairs", () => {
  const tiled = core.run({ mode: "tile", settings: { amplicons: 2 }, fasta: genomes, primers: "" });
  const pairs = tiled.result.rows.map((r) => `p${r.amplicon}\t${r.forward_seq}\t${r.reverse_seq}`).join("\n");
  const out = core.run({ mode: "check", settings: {}, fasta: genomes, primers: `${pairs}\n` });
  assert.equal(out.result.mode, "check");
  assert.deepEqual(out.result.rows.map((r) => r.name), ["p1", "p2"]);
});

test("reports the same validation messages as the server", () => {
  assert.equal(core.run({ mode: "tile", settings: { amplicons: 500 }, fasta: genomes, primers: "" }).error,
    "The web app tiles up to 200 amplicons.");
  assert.match(core.run({ mode: "tile", settings: { amplicons: 4, aligner: "muscle" }, fasta: genomes, primers: "" }).error,
    /MUSCLE is not installed/);
  assert.match(core.run({ mode: "tile", settings: { amplicons: 4 }, fasta: "hello", primers: "" }).error, /./);
});

test("aligns like the command line, so it no longer warns about the method", () => {
  const out = core.run({ mode: "tile", settings: { amplicons: 4 }, fasta: genomes, primers: "" });
  assert.equal(out.result.exit_code, 0);
  assert.ok(!out.result.warnings.some((w) => /FFT-NS|command line/.test(w)), out.result.warnings.join(" | "));
});

test("a MAFFT failure is reported clearly, out of memory included", async () => {
  const failing = await createEngineCore({
    loadPyodide: () => loadPyodide(),
    fetchBytes: async (rel) => new Uint8Array(fs.readFileSync(path.join(APP, rel))),
    mafftFactories,
    mafft: { MAFFT_PROGRAMS: mafft.MAFFT_PROGRAMS, makeMafftRunner: () => () => { throw new Error("MAFFT failed: Cannot allocate memory"); } },
  });
  const out = failing.run({ mode: "tile", settings: { amplicons: 4 }, fasta: genomes, primers: "" });
  assert.equal(out.result.exit_code, 1);
  assert.match(out.result.error, /ran out of memory/);
});
