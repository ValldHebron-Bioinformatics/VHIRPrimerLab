/* VHIRPrimerLab landing page behaviour. Data comes from data.js (a real run of the README example). */
(() => {
  'use strict';
  const VP = window.VP;
  const $ = (s, el = document) => el.querySelector(s);
  const $$ = (s, el = document) => [...el.querySelectorAll(s)];
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

  /* Swap content with a short crossfade */
  function swap(el, fn) {
    if (reduced) { fn(); return; }
    el.classList.add('is-swapping');
    setTimeout(() => { fn(); el.classList.remove('is-swapping'); }, 160);
  }

  /* ---------- Entrances ---------- */
  $$('[data-reveal]').forEach((sec) => $$('.r', sec).forEach((el, i) => el.style.setProperty('--i', i)));
  if ('IntersectionObserver' in window) {
    const io = new IntersectionObserver((entries) => {
      entries.forEach((e) => { if (e.isIntersecting) { e.target.classList.add('in'); io.unobserve(e.target); } });
    }, { rootMargin: '0px 0px -8% 0px', threshold: 0.05 });
    $$('[data-reveal]').forEach((el) => io.observe(el));
  } else {
    $$('[data-reveal]').forEach((el) => el.classList.add('in'));
  }

  /* ---------- Hero: live coverage check (same rules as vhirprimerlab/scoring.py) ---------- */
  const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };

  function match(site, seq, three, maxMm, window) {
    const n = site.length, w = Math.min(window, n);
    const inWin = (i) => (three === 'right' ? i >= n - w : i < w);
    const marks = [];
    let mm = 0, inW = 0;
    for (let i = 0; i < n; i++) {
      const s = seq.bases[i];
      const bad = !IUPAC[site[i]].includes(s) || !'ACGT'.includes(s);
      if (bad) { mm++; if (inWin(i)) inW++; }
      marks.push(bad);
    }
    for (const i of seq.ins) {
      mm++;
      if (three === 'right' ? i - 1 >= n - w : i <= w - 1) inW++;
    }
    const reason = inW ? '3′ end' : mm > maxMm ? `${mm} mismatches` : null;
    return { ok: !reason, reason, marks, inWin };
  }

  function siteHtml(site, seq, m) {
    let out = '';
    for (let i = 0; i < site.length; i++) {
      if (seq && seq.ins.includes(i)) out += '<b class="ins" title="extra base"></b>';
      const cls = [seq && m.marks[i] ? 'x' : '', m.inWin(i) ? 'w' : ''].join(' ').trim();
      const ch = !seq ? site[i] : m.marks[i] ? seq.bases[i] : '·';
      out += `<b${cls ? ` class="${cls}"` : ''}>${ch}</b>`;
    }
    return `<code translate="no" class="site">${out}</code>`;
  }

  const cov = $('#coverage');
  if (cov && VP) {
    const state = { region: 'cons', mm: 2, win: 5 };
    let last = {};
    const body = $('#cov-body');

    function render() {
      const p = VP.pairs.find((x) => x.region === state.region);
      const fwdWin = { inWin: (i) => i >= p.fwd.length - Math.min(state.win, p.fwd.length), marks: [] };
      const revWin = { inWin: (i) => i < Math.min(state.win, p.revSite.length), marks: [] };
      let rows = `<tr class="primer"><td class="id"><code translate="no">primer</code></td><td class="c-fwd">${siteHtml(p.fwd, null, fwdWin)}</td><td class="c-rev">${siteHtml(p.revSite, null, revWin)}</td><td><span class="amp-t">${p.ampLen} bp amplicon</span></td></tr>`;
      let ok = 0;
      const now = {};
      p.fwdSeqs.forEach((fs, k) => {
        const rs = p.revSeqs[k];
        const f = match(p.fwd, fs, 'right', state.mm, state.win);
        const r = match(p.revSite, rs, 'left', state.mm, state.win);
        const pass = f.ok && r.ok;
        if (pass) ok++;
        const why = [f.reason && `fwd ${f.reason}`, r.reason && `rev ${r.reason}`].filter(Boolean).join(' · ');
        now[fs.id] = pass;
        const changed = last.region === state.region && last[fs.id] !== undefined && last[fs.id] !== pass;
        rows += `<tr${changed ? ' class="flash"' : ''}><td class="id"><code translate="no">${esc(fs.id)}</code></td><td class="c-fwd">${siteHtml(p.fwd, fs, f)}</td><td class="c-rev">${siteHtml(p.revSite, rs, r)}</td><td class="res ${pass ? 'ok' : 'no'}"><span class="res-t">${pass ? 'amplified' : esc(why)}</span></td></tr>`;
      });
      body.innerHTML = rows;
      last = { region: state.region, ...now };
      const n = p.fwdSeqs.length;
      $('#cov-pct').textContent = `${Math.round((100 * ok) / n)}%`;
      $('#cov-n').textContent = `${ok} of ${n} amplified`;
    }

    $$('[data-region]', cov).forEach((b) => b.addEventListener('click', () => {
      state.region = b.dataset.region;
      $$('[data-region]', cov).forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
      render();
    }));
    $$('button[data-strand]', cov).forEach((b) => b.addEventListener('click', () => {
      cov.dataset.strand = b.dataset.strand;
      $$('button[data-strand]', cov).forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
    }));
    ['mm', 'win'].forEach((key) => {
      const input = $('#' + key), out = $('#' + key + '-out');
      input.addEventListener('input', () => { state[key] = +input.value; out.textContent = input.value; render(); });
    });
    render();
  }

  /* ---------- Why: slow crossfade of micrographs ---------- */
  const stack = $('.fade-stack');
  if (stack) {
    const imgs = $$('img', stack);
    const caps = ['Influenza A H1N1, transmission electron micrograph', 'Avian influenza A H5N1 in cultured cells', 'SARS-CoV-2 particles inside an infected cell', 'MERS coronavirus on the surface of a cell'];
    let i = 0, timer = null;
    const step = () => {
      imgs[i].classList.remove('on');
      i = (i + 1) % imgs.length;
      imgs[i].classList.add('on');
      $('#fade-cap').textContent = caps[i];
    };
    let paused = false;
    const start = () => { if (!timer && !reduced && !paused) timer = setInterval(step, 5600); };
    const stop = () => { clearInterval(timer); timer = null; };
    document.addEventListener('visibilitychange', () => (document.hidden ? stop() : start()));
    const btn = $('[data-motion="fade"]');
    if (reduced) btn.hidden = true;
    btn.addEventListener('click', () => {
      paused = !paused;
      btn.setAttribute('aria-pressed', String(paused));
      btn.textContent = paused ? 'Play' : 'Pause';
      stack.classList.toggle('paused', paused);
      paused ? stop() : start();
    });
    start();
  }

  /* ---------- Bench: play the clip only while it is on screen ---------- */
  const clip = $('.bench video');
  if (clip) {
    const btn = $('[data-motion="clip"]');
    let paused = reduced, visible = false;
    const sync = () => {
      btn.setAttribute('aria-pressed', String(paused));
      btn.textContent = paused ? 'Play' : 'Pause';
      if (visible && !paused) { clip.preload = 'auto'; clip.play().catch(() => {}); } else clip.pause();
    };
    btn.addEventListener('click', () => { paused = !paused; sync(); });
    if ('IntersectionObserver' in window) {
      new IntersectionObserver(([e]) => { visible = e.isIntersecting; sync(); }, { threshold: 0.25 }).observe(clip);
    }
    sync();
  }

  /* ---------- How it works: one slider, four stages ---------- */
  const stages = [
    ['First, line up every sequence', 'MAFFT or MUSCLE 5 aligns the family. Sequences far from the median length are set aside first, so a truncated genome cannot skew the alignment.', 'Align'],
    ['Then, take the majority base', 'Each column votes for its most common base. Columns that are gaps in at least half the sequences are left out of the consensus.', 'Consensus'],
    ['Next, ask primer3 for pairs', 'primer3 designs up to 200 candidate pairs on the consensus, keeping the whole amplicon inside each region and inside your Tm, GC and length limits.', 'Design'],
    ['Finally, test each pair on each strain', 'Every candidate is scored against every sequence, degenerate bases are added where they rescue strains, and pairs are ranked by coverage, then penalty.', 'Re-check'],
  ];

  function alnRows(withCons) {
    const w = VP.win;
    const varSet = new Set(w.var), siteSet = new Set(w.site);
    const line = (s, mark) => [...s].map((c, i) => {
      const cls = c === '-' ? 'g' : mark && siteSet.has(i) ? 's' : withCons && varSet.has(i) ? 'v' : '';
      return `<b${cls ? ` class="${cls}"` : ''}>${c}</b>`;
    }).join('');
    let h = w.ids.map((id, k) => `<div><code translate="no" class="aid">${esc(id)}</code><code translate="no">${line(w.rows[k], false)}</code></div>`).join('');
    if (withCons) h += `<div class="cons"><code translate="no" class="aid">cons</code><code translate="no">${line(w.cons, true)}</code></div>`;
    return `<div class="aln">${h}</div>`;
  }

  function stageHtml(k) {
    if (k === 0) return `<p class="frag-title">alignment.fasta <span>36 of 1,202 columns</span></p>${alnRows(false)}`;
    if (k === 1) return `<p class="frag-title">Consensus <span>grey: strains differ · magenta: primer site</span></p>${alnRows(true)}`;
    if (k === 2) {
      const r = VP.tsv.filter((x) => x.region === 'cons');
      return `<p class="frag-title">primer3 · region cons <span>200 candidate pairs</span></p>
        <div class="cov-scroll"><table class="ftable"><thead><tr><th>Forward</th><th>Reverse</th><th>bp</th><th>Penalty</th></tr></thead><tbody>
        ${r.map((x) => `<tr><td><code translate="no">${x.forward_seq}</code></td><td><code translate="no">${x.reverse_seq}</code></td><td>${x.amplicon_len}</td><td>${x.penalty}</td></tr>`).join('')}
        </tbody></table></div>`;
    }
    return `<p class="frag-title">Re-check <span>best pair per region</span></p>
      <div class="bars">
        <div class="bar-row"><p><span>cons 61–540 · 200 of 200 pairs pass</span><strong>100%</strong></p><div class="meter"><i style="width:100%"></i></div></div>
        <div class="bar-row"><p><span>hv 701–1100 · 0 of 200 pairs pass</span><strong>30%</strong></p><div class="meter"><i class="low" style="width:30%"></i></div></div>
        <div><p class="frag-title">Why hv’s best pair misses <span>sequences</span></p>
          <ul class="why-miss" role="list"><li>rev 3′ mismatch <span>5</span></li><li>fwd 3′ mismatch <span>4</span></li><li>fwd too many mismatches <span>3</span></li><li>rev too many mismatches <span>2</span></li></ul></div>
      </div>`;
  }

  const stageInput = $('#stage');
  if (stageInput && VP) {
    const frag = $('#stage-frag'), clause = $('#clause'), note = $('#clause-note'), fill = $('#rail-fill');
    let cur = -1;
    const set = (k) => {
      if (k === cur) return;
      const first = cur === -1;
      cur = k;
      stageInput.value = k;
      stageInput.setAttribute('aria-valuetext', stages[k][2]);
      fill.style.setProperty('--p', `${(k / 3) * 100}%`);
      $$('.stops button').forEach((b) => (+b.dataset.stop === k ? b.setAttribute('aria-current', 'step') : b.removeAttribute('aria-current')));
      const apply = () => { frag.innerHTML = stageHtml(k); clause.textContent = stages[k][0]; note.textContent = stages[k][1]; };
      if (first) { apply(); return; }
      swap(clause, () => {}); swap(note, () => {}); swap(frag, apply);
    };
    stageInput.addEventListener('input', () => set(+stageInput.value));
    $$('.stops button').forEach((b) => b.addEventListener('click', () => set(+b.dataset.stop)));
    set(0);
  }

  /* ---------- Outputs: the files a run writes ---------- */
  const notes = {
    tsv: 'One row per reported pair: both sequences, Tm, GC, amplicon coordinates on your reference, coverage, and which strains it misses and why. primers.tsv holds just the names and sequences, ready to order.',
    report: 'Parameters, input and alignment statistics, one section per region and, when a region falls short, a count of why its best pair misses.',
    fasta: 'The alignment that was actually used, after outliers were dropped, so you can open it in any viewer and check the primer sites yourself.',
    exit: '0 when every region has a passing pair, 2 when one falls short but every file is written, 1 on bad input. Easy to branch on in a pipeline.',
  };
  function sheetHtml(k) {
    if (k === 'tsv') {
      const cols = ['region', 'rank', 'forward_seq', 'reverse_seq', 'amplicon_len', 'coverage_pct', 'degenerate', 'below_threshold'];
      return `<div class="sheet-head"><span>pairs.tsv</span><span>10 rows · 22 columns</span></div>
        <table class="tsv"><thead><tr>${cols.map((c) => `<th><code translate="no">${c}</code></th>`).join('')}</tr></thead><tbody>
        ${VP.tsv.map((r) => `<tr${r.below_threshold === 'True' ? ' class="below"' : ''}>${cols.map((c) => `<td><code translate="no">${esc(r[c])}</code></td>`).join('')}</tr>`).join('')}
        </tbody></table>`;
    }
    if (k === 'report') {
      return `<div class="sheet-head"><span>report.txt</span><span>excerpt</span></div><pre>Alignment
---------
Alignment columns: 1202
Consensus length: 1200
Fully conserved columns: 802
Mean majority-base frequency: 0.921

Region hv (701-1100)
--------------------
primer3 candidate pairs: 200
Pairs meeting --min-coverage (95%): 0
Best coverage: 30.0%
WARNING: no pair reached 95% coverage; the best
available pairs are reported with below_threshold=True.
Why the best pair misses sequences:
  rev 3' mismatch: 5
  fwd 3' mismatch: 4
  fwd too many mismatches: 3
  rev too many mismatches: 2

Summary
-------
Regions without a pair meeting --min-coverage (95%): hv</pre>`;
    }
    if (k === 'fasta') {
      return `<div class="sheet-head"><span>alignment.fasta</span><span>10 records · 1,202 columns</span></div><pre>${VP.fasta.map(([id, s]) => `&gt;${esc(id)}\n${esc(s)}…`).join('\n')}</pre>`;
    }
    return `<div class="sheet-head"><span>Terminal</span><span>bash</span></div><pre><span class="k">$</span> vhirprimerlab -i examples/family.fasta -o out \\
    --region cons:61-540 --region hv:701-1100 \\
    --amplicon-min 100 --amplicon-max 300 --allow-degenerate
<span class="d">INFO: Loaded 10 sequences from examples/family.fasta
INFO: Aligning 10 sequences with mafft
INFO: Region cons: 200 of 200 candidate pairs meet 95% coverage</span>
<span class="w">WARNING: Region hv: no pair reached 95% coverage (best 30.0%)
WARNING: Regions without a pair meeting --min-coverage: hv</span>
<span class="d">INFO: Wrote results to out</span>
<span class="k">$</span> echo $?
2</pre>`;
  }
  const sheet = $('#sheet');
  if (sheet && VP) {
    const note = $('#out-note');
    const show = (k, animate) => {
      const apply = () => { sheet.innerHTML = sheetHtml(k); sheet.classList.toggle('dark', k === 'exit'); };
      if (animate) { swap(sheet, apply); swap(note, () => (note.textContent = notes[k])); }
      else { apply(); note.textContent = notes[k]; }
      $$('.out-name').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.sheet === k)));
    };
    $$('.out-name').forEach((b) => b.addEventListener('click', () => { if (b.getAttribute('aria-pressed') !== 'true') show(b.dataset.sheet, true); }));
    show('tsv', false);
  }

  /* ---------- Options: filterable command list ---------- */
  const q = $('#opt-q');
  if (q) {
    const items = $$('#opt-list li');
    const select = (li) => items.forEach((x) => {
      x.classList.toggle('sel', x === li);
      $('button', x).setAttribute('aria-expanded', String(x === li));
    });
    const visible = () => items.filter((x) => !x.hidden);
    items.forEach((li) => {
      $('button', li).addEventListener('click', () => select(li));
      li.addEventListener('pointerenter', (e) => { if (e.pointerType === 'mouse') select(li); });
    });
    q.addEventListener('input', () => {
      const t = q.value.trim().toLowerCase();
      items.forEach((li) => (li.hidden = !!t && !li.textContent.toLowerCase().includes(t)));
      const v = visible();
      $('#opt-empty').hidden = v.length > 0;
      if (v.length && !v.some((x) => x.classList.contains('sel'))) select(v[0]);
    });
    q.addEventListener('keydown', (e) => {
      if (e.key !== 'ArrowDown' && e.key !== 'ArrowUp') return;
      e.preventDefault();
      const v = visible();
      if (!v.length) return;
      const i = v.findIndex((x) => x.classList.contains('sel'));
      select(v[(i + (e.key === 'ArrowDown' ? 1 : v.length - 1)) % v.length]);
    });
    select(items[0]);
  }

  /* ---------- FAQ: question index, one answer ---------- */
  const faq = $('#faq-q');
  if (faq) {
    const pane = $('#faq-pane');
    const btns = $$('button', faq);
    const wide = matchMedia('(min-width: 861px)');
    const pick = (b, animate) => {
      btns.forEach((x) => x.setAttribute('aria-expanded', String(x === b)));
      const text = $('#' + b.getAttribute('aria-controls')).textContent;
      if (animate) swap(pane, () => (pane.textContent = text));
      else pane.textContent = text;
    };
    btns.forEach((b) => {
      b.addEventListener('click', () => {
        if (!wide.matches && b.getAttribute('aria-expanded') === 'true') { b.setAttribute('aria-expanded', 'false'); return; }
        pick(b, wide.matches);
      });
      b.addEventListener('focus', () => { if (wide.matches && b.getAttribute('aria-expanded') !== 'true') pick(b, true); });
      b.addEventListener('pointerenter', (e) => { if (wide.matches && e.pointerType === 'mouse' && b.getAttribute('aria-expanded') !== 'true') pick(b, true); });
    });
    pick(btns[0], false);
  }

  /* ---------- Start: command builder ---------- */
  const form = $('#builder');
  if (form) {
    let aligner = 'mafft', mode = 'tile';
    const quote = (s) => (/^[\w./:@%+=,-]+$/.test(s) ? s : `'${s.replace(/'/g, `'\\''`)}'`);
    const regionRe = /^([A-Za-z0-9_.-]+:)?\d+-\d+$/;
    const err = $('#dlg-err');

    function build() {
      const f = form.elements;
      const problems = [];
      const mark = (el, bad) => (bad ? el.setAttribute('aria-invalid', 'true') : el.removeAttribute('aria-invalid'));
      const input = f.input.value.trim(), out = f.out.value.trim() || 'out';
      mark(f.input, !input); if (!input) problems.push('Name the input FASTA.');
      let parts;
      if (mode === 'tile') {
        const n = parseInt(f.n.value, 10);
        [f.primers, f.regions, f.amin, f.amax].forEach((el) => mark(el, false));
        mark(f.n, !(n > 0)); if (!(n > 0)) problems.push('The number of amplicons must be 1 or more.');
        parts = ['vhirprimerlab', '-i', quote(input || 'genomes.fasta'), '-n', String(n > 0 ? n : 10), '-o', quote(out)];
        if (aligner !== 'mafft') parts.push('--aligner', aligner);
        if (!f.deg.checked) parts.push('--no-degenerate');
      } else if (mode === 'check') {
        const primers = f.primers.value.trim();
        mark(f.primers, !primers); if (!primers) problems.push('Name the primers file.');
        [f.regions, f.amin, f.amax, f.n].forEach((el) => mark(el, false));
        parts = ['vhirprimerlab', 'check', '-i', quote(input || 'family.fasta'), '--primers', quote(primers || 'assays.tsv'), '-o', quote(out)];
        if (aligner !== 'mafft') parts.push('--aligner', aligner);
      } else {
        const regions = f.regions.value.trim().split(/\s+/).filter(Boolean);
        const amin = parseInt(f.amin.value, 10), amax = parseInt(f.amax.value, 10);
        mark(f.primers, false); mark(f.n, false);
        const badReg = regions.filter((r) => !regionRe.test(r));
        mark(f.regions, badReg.length > 0); if (badReg.length) problems.push(`Region “${badReg[0]}” should look like cons:61-540.`);
        const badAmp = !(amin > 0) || !(amax > 0) || amin > amax;
        mark(f.amin, badAmp); mark(f.amax, badAmp); if (badAmp) problems.push('Amplicon min and max must be positive, with min ≤ max.');
        parts = ['vhirprimerlab', '-i', quote(input || 'family.fasta'), '-o', quote(out)];
        if (aligner !== 'mafft') parts.push('--aligner', aligner);
        regions.filter((r) => regionRe.test(r)).forEach((r) => parts.push('--region', r));
        parts.push('--amplicon-min', String(amin || 100), '--amplicon-max', String(amax || 300));
        if (f.deg.checked) parts.push('--allow-degenerate');
      }
      err.textContent = problems.join(' ');
      $('#cmd').textContent = parts.join(' ');
      return problems.length === 0;
    }

    $$('[data-mode]', form).forEach((b) => b.addEventListener('click', () => {
      mode = b.dataset.mode;
      $$('[data-mode]', form).forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
      $$('[data-only]', form).forEach((el) => (el.hidden = !el.dataset.only.split(' ').includes(mode)));
      $('#dlg-title').textContent = { tile: 'New tiling', design: 'New region run', check: 'New check' }[mode];
      build();
    }));
    form.addEventListener('input', build);
    $$('[data-aligner]', form).forEach((b) => b.addEventListener('click', () => {
      aligner = b.dataset.aligner;
      $$('[data-aligner]', form).forEach((x) => x.setAttribute('aria-pressed', String(x === b)));
      build();
    }));
    let timer;
    form.addEventListener('submit', async (e) => {
      e.preventDefault();
      const msg = $('#copied');
      if (!build()) {
        msg.textContent = 'Fix the highlighted field, then copy again.';
        const bad = form.querySelector('[aria-invalid="true"]');
        if (bad) bad.focus();
        return;
      }
      const text = $('#cmd').textContent;
      try {
        await navigator.clipboard.writeText(text);
        msg.textContent = 'Copied. Paste it in your terminal.';
      } catch {
        const range = document.createRange();
        range.selectNodeContents($('#cmd'));
        const sel = getSelection(); sel.removeAllRanges(); sel.addRange(range);
        msg.textContent = 'Selected. Press Ctrl+C to copy.';
      }
      clearTimeout(timer);
      timer = setTimeout(() => (msg.textContent = ''), 4000);
    });
    build();
  }
})();
