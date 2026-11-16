// The two ways the page can run an analysis, behind one interface:
//   engine.name            'server' | 'browser'
//   engine.limits()        Promise of the upload limits and installed aligners
//   engine.run(job, onEvent) -> {promise, cancel}
//        job: {mode, settings, fasta (text), primers (text), background (text, optional FASTA)}
//        onEvent({phase: 'loading'|'queued'|'running', stage, progress, queuePosition, jobId})
//        promise resolves {status: 'done'|'failed', result, error, downloads: [{name, href}], zip}
//        and rejects with an Error (err.cancelled when cancelled, err.invalid for bad input)
// Classic script: defines window.VhirEngines.
(function (root) {
  'use strict';

  const FILES = {
    tile: ['primers.tsv', 'amplicons.tsv', 'scheme.bed', 'report.txt', 'alignment.fasta'],
    regions: ['primers.tsv', 'pairs.tsv', 'report.txt', 'alignment.fasta'],
    check: ['check.tsv', 'report.txt', 'alignment.fasta'],
  };
  const BROWSER_LIMITS = {
    max_upload_mb: 20, max_seqs: 500, max_seq_len: 35000, max_amplicons: 200,
    retention_hours: null, timeout_minutes: null, aligners: ['mafft'],
  };

  function cancelledError() {
    const err = new Error('Cancelled');
    err.cancelled = true;
    return err;
  }

  // ---------- server: POST the job, then poll it ----------
  function createServerEngine(limits) {
    function follow(jobId, onEvent) {
      let stopped = false;
      let timer = null;
      let rejectFn = null;
      const promise = new Promise((resolve, reject) => {
        rejectFn = reject;
        async function poll() {
          if (stopped) return;
          let st;
          try {
            const r = await fetch(`/api/jobs/${jobId}`);
            if (r.status === 404) {
              reject(new Error(`That job is no longer available: results are kept for ${limits.retention_hours || 24} hours.`));
              return;
            }
            st = await r.json();
          } catch (err) {
            timer = setTimeout(poll, 3000);
            return;
          }
          if (stopped) return;
          if (st.status === 'done' || st.status === 'failed') {
            const m = (st.result && st.result.mode) || 'tile';
            resolve({
              status: st.status, result: st.result, error: st.error, jobId,
              downloads: FILES[m].map((name) => ({ name, href: `/api/jobs/${jobId}/files/${name}` })),
              zip: { name: `vhirprimerlab-${jobId.slice(0, 8)}.zip`, href: `/api/jobs/${jobId}/zip` },
            });
            return;
          }
          onEvent({ phase: st.status === 'queued' ? 'queued' : 'running', stage: st.stage,
            progress: st.progress, queuePosition: st.queue_position, jobId });
          timer = setTimeout(poll, 1000);
        }
        poll();
      });
      return {
        promise,
        cancel() {
          stopped = true;
          clearTimeout(timer);
          // stop the job on the server too, so it frees its worker slot
          fetch(`/api/jobs/${jobId}`, { method: 'DELETE' }).catch(() => {});
          if (rejectFn) rejectFn(cancelledError());
        },
      };
    }

    return {
      name: 'server',
      limits: async () => limits,
      resume: follow,
      run(job, onEvent) {
        let inner = null;
        let cancelled = false;
        const promise = (async () => {
          const body = new FormData();
          body.append('mode', job.mode);
          body.append('params', JSON.stringify(job.settings || {}));
          body.append('fasta', new Blob([job.fasta], { type: 'text/plain' }), 'input.fasta');
          if (job.mode === 'check') body.append('primers_text', job.primers || '');
          if (job.background) body.append('background_text', job.background);
          let response;
          try {
            response = await fetch('/api/jobs', { method: 'POST', body });
          } catch (err) {
            throw new Error('Could not reach the server. Check your connection and try again.');
          }
          const data = await response.json().catch(() => ({}));
          if (!response.ok) {
            const err = new Error(data.detail || 'The server could not start the job.');
            err.invalid = response.status === 400;
            throw err;
          }
          if (cancelled) throw cancelledError();
          onEvent({ phase: 'queued', jobId: data.id, progress: 0 });
          inner = follow(data.id, onEvent);
          return inner.promise;
        })();
        return {
          promise,
          cancel() {
            cancelled = true;
            if (inner) inner.cancel();
          },
        };
      },
    };
  }

  // ---------- browser: a worker running vhirprimerlab in Pyodide ----------
  function createBrowserEngine(workerUrl) {
    let worker = null;
    let urls = [];

    function startWorker() {
      if (!worker) worker = new Worker(workerUrl || 'engine-worker.js');
      return worker;
    }
    function stopWorker() {
      if (worker) worker.terminate();
      worker = null;
    }

    return {
      name: 'browser',
      limits: async () => BROWSER_LIMITS,
      warmup() {
        startWorker().postMessage({ type: 'warmup' });
      },
      run(job, onEvent) {
        let rejectFn = null;
        const promise = new Promise((resolve, reject) => {
          rejectFn = reject;
          const w = startWorker();
          w.onmessage = (event) => {
            const m = event.data || {};
            if (m.type === 'loading') onEvent({ phase: 'loading', stage: m.stage, progress: m.progress });
            else if (m.type === 'progress') onEvent({ phase: 'running', stage: m.stage, progress: m.progress });
            else if (m.type === 'invalid') {
              const err = new Error(m.message);
              err.invalid = true;
              reject(err);
            } else if (m.type === 'error') {
              stopWorker(); // a crashed Pyodide cannot be reused; the next run starts a fresh one
              reject(new Error(m.message));
            } else if (m.type === 'done') {
              urls.forEach((u) => URL.revokeObjectURL(u));
              urls = [];
              const link = (name, data, type) => {
                const href = URL.createObjectURL(new Blob([data], { type }));
                urls.push(href);
                return { name, href };
              };
              const mode = m.result.mode;
              resolve({
                status: m.result.exit_code === 1 ? 'failed' : 'done',
                result: m.result,
                error: m.result.error,
                downloads: FILES[mode].filter((n) => n in m.files).map((n) => link(n, m.files[n], 'text/plain')),
                zip: link('vhirprimerlab-results.zip', m.zip, 'application/zip'),
              });
            }
          };
          w.onerror = (event) => {
            event.preventDefault();
            stopWorker();
            reject(new Error('Could not load the analysis tools. Check your connection and reload the page.'));
          };
          w.postMessage({ type: 'run', mode: job.mode, settings: job.settings || {},
            fasta: job.fasta, primers: job.primers || '', background: job.background || '' });
        });
        return {
          promise,
          cancel() {
            stopWorker();
            if (rejectFn) rejectFn(cancelledError());
          },
        };
      },
    };
  }

  // The server engine when this page is served by vhirprimerlab web / Docker (the server marks the
  // page with <meta name="vhirprimerlab-server">), else the browser.
  async function pickEngine() {
    const marked = typeof document !== 'undefined'
      && document.querySelector('meta[name="vhirprimerlab-server"]');
    if (!marked) return createBrowserEngine(); // GitHub Pages or any static copy of the page
    const controller = typeof AbortController === 'function' ? new AbortController() : null;
    const timer = setTimeout(() => controller && controller.abort(), 2000);
    try {
      const r = await fetch('/api/limits', controller ? { signal: controller.signal } : undefined);
      if (r.ok) {
        const limits = await r.json();
        if (limits && typeof limits.max_seqs === 'number') return createServerEngine(limits);
      }
    } catch (err) {
      // no server: GitHub Pages or a plain file server
    } finally {
      clearTimeout(timer);
    }
    return createBrowserEngine();
  }

  root.VhirEngines = { pickEngine, createServerEngine, createBrowserEngine, FILES, BROWSER_LIMITS };
})(typeof window !== 'undefined' ? window : globalThis);
