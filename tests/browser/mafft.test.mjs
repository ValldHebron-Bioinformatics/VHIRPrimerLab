// node --test tests/browser/mafft.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { createRequire } from "node:module";
import { synthGenomes } from "./synth.mjs";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "../..");
const ENGINE = path.join(ROOT, "site/app/engine");
const require = createRequire(import.meta.url);
const { MAFFT_PROGRAMS, MAFFT_VERSION, autoMethod, makeMafftRunner, parseFasta } =
  require(path.join(ROOT, "vhirprimerlab/web/static/mafft-runner.js"));

function loadPrograms() {
  const programs = {};
  for (const name of MAFFT_PROGRAMS) {
    programs[name] = {
      factory: require(path.join(ENGINE, `mafft-${name}.js`)),
      wasm: fs.readFileSync(path.join(ENGINE, `mafft-${name}.wasm`)),
    };
  }
  return programs;
}

const input = fs.readFileSync(path.join(ROOT, "examples/genomes.fasta"), "utf8");

test("picks the method mafft --auto picks", () => {
  const name = (n, len) => autoMethod(n, len).name;
  assert.equal(name(16, 800), "L-INS-i");
  assert.equal(name(99, 2999), "L-INS-i");
  assert.equal(name(100, 300), "L-INS-i");
  assert.equal(name(199, 999), "L-INS-i");
  assert.equal(name(100, 2999), "FFT-NS-i");
  assert.equal(name(200, 999), "FFT-NS-i");
  assert.equal(name(499, 9999), "FFT-NS-i");
  assert.equal(name(500, 100), "FFT-NS-2");
  assert.equal(name(10, 10000), "FFT-NS-2");

  // the exact programs and arguments mafft 7.526 runs (traced with bash -x)
  assert.deepEqual(autoMethod(16, 800).steps, [
    { program: "tbfast", args: "_ -u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1 -L _ -+ 16 -W 0.00001 -V -1.53 -s 0.0 -O -C 0 -b 62 -f -1.53 -Q 100.0 -h 0 -F -l 2.7 -X 0.1 -i infile".split(" ") },
    { program: "dvtditr", args: "-W 0.00001 -E 0.0 -s 0.0 -C 0 -t 0 -F -l 2.7 -z 50 -b 62 -f -1.53 -Q 100.0 -h 0 -I 16 -X 0.1 -p BAATARI2 -K 0 -i pre".split(" ") },
  ]);
  assert.ok(autoMethod(150, 900).steps[0].args.join(" ").includes("-+ 2 "));
  assert.deepEqual(autoMethod(120, 8000).steps, [
    { program: "disttbfast", args: "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 -i infile".split(" "), stdout: "pre" },
    { program: "dndpre", args: "-y hat2 -b 62 -M 2 -C 0 -i pre".split(" ") },
    { program: "dvtditr", args: "-W 0.00001 -E 0.0 -s 0.0 -C 0 -t 0 -F -z 50 -b 62 -f -1.53 -Q 100.0 -h 0 -I 2 -X 0.1 -p BAATARI2 -K 0 -i pre".split(" ") },
  ]);
  assert.deepEqual(autoMethod(160, 13444).steps, [
    { program: "disttbfast", args: "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 -i infile".split(" "), stdout: "pre" },
  ]);
  // two sequences: one guide-tree cycle
  assert.ok(autoMethod(2, 4000).steps[0].args.join(" ").includes("-E 1 "));
  assert.ok(autoMethod(2, 20000).steps[0].args.join(" ").includes("-E 1 "));
});

test("aligns the example genomes synchronously, twice in one process", () => {
  const runMafft = makeMafftRunner(loadPrograms());
  for (let round = 0; round < 2; round += 1) {
    const aligned = runMafft(input);          // a plain value, not a Promise
    assert.equal(typeof aligned, "string");
    const rows = parseFasta(aligned);
    const source = parseFasta(input);
    assert.equal(rows.length, 24);
    assert.deepEqual(rows.map((r) => r.id), source.map((r) => r.id));
    assert.deepEqual(new Set(rows.map((r) => r.seq.length)), new Set([9000]));
    rows.forEach((r, i) => assert.equal(r.seq.replace(/-/g, "").toUpperCase(), source[i].seq.toUpperCase()));
  }
});

test("reports MAFFT failures as errors", () => {
  const runMafft = makeMafftRunner(loadPrograms());
  assert.throws(() => runMafft("no fasta here\n"), /MAFFT failed/);
});

// The command-line MAFFT, when it is the version the browser runs
function cliMafft() {
  const dirs = [process.env.PYTHON && path.dirname(process.env.PYTHON), ...(process.env.PATH || "").split(path.delimiter)];
  for (const dir of dirs.filter(Boolean)) {
    const bin = path.join(dir, "mafft");
    if (!fs.existsSync(bin)) continue;
    const v = spawnSync(bin, ["--version"], { encoding: "utf8" });
    return { bin, version: `${v.stdout}${v.stderr}`.trim() };
  }
  return null;
}
const cli = cliMafft();
const sameVersion = cli && cli.version.startsWith(`v${MAFFT_VERSION} `);

for (const [label, nseq, length] of [
  ["L-INS-i, 16 iterations", 16, 800],
  ["L-INS-i, 2 iterations", 100, 300],
  ["FFT-NS-i", 24, 3200],
  ["FFT-NS-2", 4, 10200],
]) {
  test(`gives the same alignment as mafft --auto: ${label}`,
    { skip: !sameVersion && `needs mafft v${MAFFT_VERSION} on the PATH (found ${cli ? cli.version : "none"})` },
    () => {
      const fasta = synthGenomes(nseq, length, 7);
      const dir = fs.mkdtempSync(path.join(os.tmpdir(), "vhirprimerlab-mafft-"));
      fs.writeFileSync(path.join(dir, "in.fa"), fasta);
      const native = spawnSync(cli.bin, ["--auto", "--quiet", path.join(dir, "in.fa")], { encoding: "utf8", maxBuffer: 1 << 28 });
      assert.equal(native.status, 0, native.stderr);
      const want = parseFasta(native.stdout);
      const got = parseFasta(makeMafftRunner(loadPrograms())(fasta));
      assert.deepEqual(got.map((r) => r.id), want.map((r) => r.id));
      assert.deepEqual(got.map((r) => r.seq.toUpperCase()), want.map((r) => r.seq.toUpperCase()));
    });
}
