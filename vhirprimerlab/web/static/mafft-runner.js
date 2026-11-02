// MAFFT compiled to WebAssembly, run synchronously so Python code can call it like a function.
// It runs the same programs, with the same arguments, as `mafft --auto` on the command line, so
// the browser gives the same alignment.
// Classic script: importScripts() in a worker defines self.VhirMafft; require() in Node.
(function (root) {
'use strict';

const MAFFT_VERSION = '7.526';
const MAFFT_PROGRAMS = ['disttbfast', 'tbfast', 'dndpre', 'dvtditr'];

// What `mafft --auto` (MAFFT 7.526) runs for nseq nucleotide sequences, the longest maxLen bases
// long, up to 500 sequences: its size rules and the arguments it passes (traced with bash -x).
function autoMethod(nseq, maxLen) {
  const words = (s) => s.split(' ');
  const cycles = nseq === 2 ? 1 : 2; // two sequences: one guide-tree cycle
  const disttbfast = { program: 'disttbfast', stdout: 'pre',
    args: words(`-q 0 -E ${cycles} -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 -i infile`) };
  const refine = (iterations, local) => ({ program: 'dvtditr',
    args: words(`-W 0.00001 -E 0.0 -s 0.0 -C 0 -t 0 -F ${local}-z 50 -b 62 -f -1.53 -Q 100.0 -h 0 -I ${iterations} -X 0.1 -p BAATARI2 -K 0 -i pre`) });
  const linsi = (iterations) => ({ name: 'L-INS-i', steps: [
    { program: 'tbfast',
      args: words(`_ -u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1 -L _ -+ ${iterations} -W 0.00001 -V -1.53 -s 0.0 -O -C 0 -b 62 -f -1.53 -Q 100.0 -h 0 -F -l 2.7 -X 0.1 -i infile`) },
    refine(iterations, '-l 2.7 '),
  ] });

  // 1000 iterations, capped at 16 by the BAATARI2 strategy
  if (maxLen < 3000 && nseq < 100) return linsi(16);
  if (maxLen < 1000 && nseq < 200) return linsi(2);
  if (maxLen < 10000 && nseq < 500) {
    return { name: 'FFT-NS-i', steps: [
      disttbfast,
      { program: 'dndpre', args: words('-y hat2 -b 62 -M 2 -C 0 -i pre') },
      refine(2, ''),
    ] };
  }
  return { name: 'FFT-NS-2', steps: [disttbfast] };
}

function parseFasta(text) {
  const records = [];
  let current = null;
  for (const raw of text.split(/\r?\n/)) {
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith('>') || line.startsWith('=')) {
      current = { id: line.slice(1).split(/\s+/)[0], seq: '' };
      records.push(current);
    } else if (current) {
      current.seq += line;
    }
  }
  return records;
}

// Runs one MAFFT program in a fresh module, in a working folder holding `files`
// ({name: Uint8Array|string}); returns the folder's files afterwards.
function runProgram(program, step, files) {
  let mod = null;
  const out = [];
  const err = [];
  // Built with -sWASM_ASYNC_COMPILATION=0: with the bytes given, the instance starts
  // synchronously and the pre-run hook hands us the module before the factory returns.
  const started = program.factory({
    wasmBinary: program.wasm,
    noInitialRun: true,
    print: (line) => out.push(line),
    printErr: (line) => err.push(line),
    preRun: [(m) => { mod = m; }],
  });
  if (started && typeof started.catch === 'function') started.catch(() => {});
  if (!mod) throw new Error('MAFFT could not start');
  mod.FS.mkdir('/work');
  mod.FS.chdir('/work');
  for (const [name, data] of Object.entries(files)) mod.FS.writeFile(name, data);
  let code = 0;
  try {
    // -i <file>: MAFFT rewinds its input, which a stream cannot do
    code = mod.callMain(step.args.slice());
  } catch (e) {
    if (e && e.name === 'ExitStatus') code = e.status;
    else throw new Error(`MAFFT failed: ${e && e.message ? e.message : e}`);
  }
  if (code !== 0) {
    throw new Error(`MAFFT failed: ${err.filter(Boolean).slice(-3).join(' ') || `exit code ${code}`}`);
  }
  const result = {};
  for (const name of mod.FS.readdir('/work')) {
    if (name === '.' || name === '..') continue;
    if (mod.FS.isFile(mod.FS.stat(name).mode)) result[name] = mod.FS.readFile(name);
  }
  if (step.stdout) result[step.stdout] = out.length ? `${out.join('\n')}\n` : '';
  return result;
}

// programs: {disttbfast|tbfast|dndpre|dvtditr: {factory, wasm}} — the module factories
// (mafft-<name>.js) and their .wasm bytes.
// Returns runMafft(fastaText) -> aligned FASTA text, in the input order.
function makeMafftRunner(programs) {
  return function runMafft(fastaText) {
    const records = parseFasta(fastaText);
    const maxLen = records.reduce((m, r) => Math.max(m, r.seq.length), 0);
    const method = autoMethod(records.length, maxLen);
    // the files mafft's script starts with
    let files = { infile: fastaText, pre: '', 'hat3.seed': '' };
    method.steps.forEach((step, i) => {
      files = runProgram(programs[step.program], step, files);
      // after disttbfast, mafft's script runs `mv hat3.seed hat3`
      if (i === 0 && step.program === 'disttbfast' && 'hat3.seed' in files) {
        files.hat3 = files['hat3.seed'];
        delete files['hat3.seed'];
      }
    });
    const text = typeof files.pre === 'string' ? files.pre : new TextDecoder().decode(files.pre || new Uint8Array());
    if (!text.trim()) throw new Error('MAFFT failed: no alignment was produced');
    return text.replace(/^=/gm, '>');
  };
}

const api = { MAFFT_VERSION, MAFFT_PROGRAMS, autoMethod, parseFasta, makeMafftRunner };
if (typeof module === 'object' && module.exports) module.exports = api;
else root.VhirMafft = api;
})(typeof self !== 'undefined' ? self : globalThis);
