// vhirprimerlab's browser engine: runs the vhirprimerlab command line inside Pyodide, with primer3 and
// MAFFT compiled to WebAssembly. Loaded as a classic Web Worker by the page, or with require()
// in Node for tests (createEngineCore).
(function (root) {
  'use strict';

  const PYODIDE_VERSION = '0.27.7';
  // served from the site itself (tools/build_site.py downloads and verifies it), relative to this worker
  const PYODIDE_URL = 'engine/pyodide/';

  // Python side: MAFFT as vhirprimerlab's "mafft" aligner, progress to the page, and one call that
  // checks a request exactly like the server (prepare_job), runs vhirprimerlab and builds the result.
  const SHIM = `
import contextlib, io, json, os, re, shutil, sys, traceback, zipfile
from pathlib import Path

import vhirprimerlab_js
from vhirprimerlab.align import register_aligner
from vhirprimerlab.cli import main, set_progress_reporter
from vhirprimerlab.io import alignment_from_records, load_fasta
from vhirprimerlab.models import PipelineError
from vhirprimerlab.web.params import Limits
from vhirprimerlab.web.prepare import prepare_job
from vhirprimerlab.web.result import OUTPUT_FILES, build_result

LIMITS = Limits()
OOM = "Your browser ran out of memory. Try fewer sequences or the command line."


def _mafft(records):
    # the same programs and arguments as \`mafft --auto\` on the command line (mafft-runner.js)
    text = "".join(f">{r.id}\\n{r.seq}\\n" for r in records)
    try:
        aligned = vhirprimerlab_js.run_mafft(text)
    except Exception as exc:  # a JavaScript error from the MAFFT module
        message = str(exc)
        if re.search(r"memory|alloc|out of bounds", message, re.IGNORECASE):
            raise PipelineError(OOM) from None
        raise PipelineError("MAFFT failed: " + message.split("MAFFT failed:")[-1].strip()) from None
    out = Path("/tmp/mafft_out.fasta")
    out.write_text(aligned)
    return alignment_from_records(load_fasta(out, allow_gaps=True))


register_aligner("mafft", _mafft)
set_progress_reporter(lambda stage, fraction: vhirprimerlab_js.progress(stage, fraction))


def run_job(mode, settings_json, fasta, primers, background=""):
    try:
        argv, files = prepare_job(mode, json.loads(settings_json), fasta.encode(), primers,
                                  LIMITS, ["mafft"], background)
    except ValueError as exc:
        return json.dumps({"error": str(exc)}), b""
    job = Path("/job")
    os.chdir("/")  # the previous run left us inside /job
    shutil.rmtree(job, ignore_errors=True)
    job.mkdir()
    for name, data in files.items():
        (job / name).write_bytes(data)
    os.chdir(job)
    log = io.StringIO()
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            code = main([*argv, "-o", "out"])
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
        except MemoryError:
            print(f"ERROR: {OOM}", file=sys.stderr)
            code = 1
        except Exception:  # reported to the user as a failed run
            traceback.print_exc()
            code = 1
    (job / "log.txt").write_text(log.getvalue())
    result = build_result(job, argv, code)
    out = job / "out"
    texts = {n: (out / n).read_text() for n in OUTPUT_FILES if (out / n).is_file()}
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, text in texts.items():
            zf.writestr(name, text)
    return json.dumps({"result": result, "files": texts}), buf.getvalue()
`;

  // loadPyodide(): a Pyodide instance; fetchBytes(relativePath): Uint8Array of a site file;
  // mafftFactories: {disttbfast|tbfast|dndpre|dvtditr: module factory}; mafft: mafft-runner.js;
  // report(stage, fraction).
  async function createEngineCore({ loadPyodide, fetchBytes, mafftFactories, mafft, report }) {
    const say = report || (() => {});
    say('Loading the analysis tools', 0.02);
    const manifest = JSON.parse(new TextDecoder().decode(await fetchBytes('engine/manifest.json')));
    const py = await loadPyodide();
    say('Loading Biopython', 0.3);
    await py.loadPackage(['biopython', 'micropip'], { messageCallback: () => {} });
    say('Loading primer3 and VHIRPrimerLab', 0.6);
    const micropip = py.pyimport('micropip');
    for (const name of [manifest.primer3, manifest.vhirprimerlab]) {
      py.FS.writeFile(`/tmp/${name}`, await fetchBytes(`engine/${name}`));
      await micropip.install.callKwargs(`emfs:/tmp/${name}`, { deps: false });
    }
    micropip.destroy();
    say('Loading MAFFT', 0.9);
    const programs = {};
    for (const name of mafft.MAFFT_PROGRAMS) {
      programs[name] = { factory: mafftFactories[name], wasm: await fetchBytes(`engine/mafft-${name}.wasm`) };
    }
    const runMafft = mafft.makeMafftRunner(programs);
    let progress = () => {};
    py.registerJsModule('vhirprimerlab_js', {
      run_mafft: (text) => runMafft(text),
      progress: (stage, fraction) => progress(stage, fraction),
    });
    py.runPython(SHIM);
    const runJob = py.globals.get('run_job');
    say('Ready', 1);

    return {
      // job: {mode, settings, fasta, primers, background}; returns {error} or {result, files, zip}
      run(job, onProgress) {
        progress = onProgress || (() => {});
        const out = runJob(job.mode, JSON.stringify(job.settings || {}), job.fasta || '', job.primers || '',
          job.background || '');
        try {
          const [text, zip] = out.toJs();
          return { ...JSON.parse(text), zip };
        } finally {
          out.destroy();
          progress = () => {};
        }
      },
    };
  }

  function friendly(err) {
    const text = String((err && err.message) || err);
    if (/memory|out of bounds|RangeError|Maximum call stack/i.test(text)) {
      return 'Your browser ran out of memory. Try fewer sequences or the command line.';
    }
    if (/fetch|NetworkError|Failed to load|importScripts|load script/i.test(text)) {
      return 'Could not load the analysis tools. Check your connection and reload the page.';
    }
    return `The analysis stopped: ${text.split('\n').filter(Boolean).slice(-1)[0] || 'unknown error'}`;
  }

  const api = { createEngineCore, friendly, PYODIDE_URL, PYODIDE_VERSION };
  if (typeof module === 'object' && module.exports) {
    module.exports = api;
    return;
  }

  // ---- running as a Web Worker ----
  let core = null;
  const post = (message, transfer) => root.postMessage(message, transfer || []);
  function ensureCore() {
    if (!core) {
      core = (async () => {
        root.importScripts(`${PYODIDE_URL}pyodide.js`, 'mafft-runner.js');
        const mafftFactories = {};
        for (const name of root.VhirMafft.MAFFT_PROGRAMS) {
          root.importScripts(`engine/mafft-${name}.js`); // defines createMafft_<name>
          mafftFactories[name] = root[`createMafft_${name}`];
        }
        return createEngineCore({
          loadPyodide: () => root.loadPyodide({ indexURL: new URL(PYODIDE_URL, root.location.href).href }),
          fetchBytes: async (rel) => {
            const r = await fetch(rel);
            if (!r.ok) throw new Error(`Failed to load ${rel} (${r.status})`);
            return new Uint8Array(await r.arrayBuffer());
          },
          mafftFactories,
          mafft: root.VhirMafft,
          report: (stage, fraction) => post({ type: 'loading', stage, progress: fraction }),
        });
      })();
      core.catch(() => { core = null; });
    }
    return core;
  }

  root.onmessage = async (event) => {
    const msg = event.data || {};
    try {
      const engine = await ensureCore();
      if (msg.type === 'warmup') {
        post({ type: 'ready' });
        return;
      }
      const out = engine.run(msg, (stage, fraction) => post({ type: 'progress', stage, progress: fraction }));
      if (out.error) {
        post({ type: 'invalid', message: out.error });
      } else {
        post({ type: 'done', result: out.result, files: out.files, zip: out.zip }, [out.zip.buffer]);
      }
    } catch (err) {
      post({ type: 'error', message: friendly(err) });
    }
  };
})(typeof self !== 'undefined' ? self : globalThis);
