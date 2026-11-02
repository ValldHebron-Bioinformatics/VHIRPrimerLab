/* VHIRPrimerLab web app: read sequences in the browser, start a job, follow it, draw the result. */
(() => {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const fmt = (n) => Number(n).toLocaleString('en-US');
  const pct = (n) => `${Number(n).toFixed(Number(n) % 1 ? 1 : 0)} %`;

  const LARGE_INPUT = 3000000; // total bases above which the browser engine asks first
  const STAGES = {
    aligning: 'Aligning your sequences…',
    'finding primers': 'Finding candidate primers along the genome…',
    scoring: 'Checking primers against every sequence…',
    'choosing amplicons': 'Choosing the best set of amplicons…',
    writing: 'Writing the files…',
    done: 'Finishing…',
  };
  const FILES = {
    tile: ['primers.tsv', 'amplicons.tsv', 'scheme.bed', 'report.txt', 'alignment.fasta'],
    regions: ['primers.tsv', 'pairs.tsv', 'report.txt', 'alignment.fasta'],
    check: ['check.tsv', 'report.txt', 'alignment.fasta'],
  };

  const state = { page: 'design', seqs: null, limits: null, job: null, map: null, result: null,
    engine: null, handle: null, hasResults: false, lastJob: null };

  // ---------- mode ----------
  const mode = () => (state.page === 'check' ? 'check' : $('use-regions').checked ? 'regions' : 'tile');

  function applyMode() {
    const m = mode();
    document.querySelectorAll('[data-modes]').forEach((el) => {
      el.hidden = !el.dataset.modes.split(' ').includes(m);
    });
    document.querySelectorAll('.modes a').forEach((a) => {
      if (a.dataset.mode === state.page) a.setAttribute('aria-current', 'page');
      else a.removeAttribute('aria-current');
    });
    const check = m === 'check';
    $('title').textContent = check ? 'Check my primers' : 'Design primers';
    $('lede').textContent = check
      ? 'Paste primer pairs you already have and see how many of your sequences each pair amplifies, with the same matching rules as design.'
      : 'Drop your genomes, choose how many amplicons you want, and get overlapping primer pairs in two pools, checked against every sequence.';
    $('run').textContent = check ? 'Check primers' : 'Design primers';
    document.title = `VHIRPrimerLab: ${check ? 'check primers' : 'tile a genome with primers'}`;
    if (!state.job) showIdle();
  }

  // ---------- sequences ----------
  function parseFasta(text) {
    const ids = [];
    const lens = [];
    let len = -1;
    for (const raw of text.split(/\r?\n/)) {
      const line = raw.trim();
      if (!line) continue;
      if (line.startsWith('>')) {
        if (len >= 0) lens.push(len);
        ids.push(line.slice(1).split(/\s+/)[0] || `sequence ${ids.length + 1}`);
        len = 0;
      } else if (len < 0) {
        throw new Error('This does not look like FASTA: each sequence needs a header line starting with ">".');
      } else {
        len += line.replace(/[^A-Za-z]/g, '').length;
      }
    }
    if (len >= 0) lens.push(len);
    if (!ids.length) throw new Error('No sequences found.');
    const sorted = [...lens].sort((a, b) => a - b);
    const mid = sorted.length >> 1;
    const median = sorted.length % 2 ? sorted[mid] : Math.round((sorted[mid - 1] + sorted[mid]) / 2);
    return { ids, lens, median };
  }

  function seqError(message) {
    $('seq-error').textContent = message || '';
    $('seq-error').hidden = !message;
  }

  function loadText(name, text) {
    let parsed;
    try {
      parsed = parseFasta(text);
    } catch (err) {
      seqError(err.message);
      return;
    }
    const L = state.limits || { max_seqs: 500, max_seq_len: 35000 };
    const longest = Math.max(...parsed.lens);
    if (parsed.ids.length > L.max_seqs) {
      seqError(`${fmt(parsed.ids.length)} sequences: the web app accepts up to ${fmt(L.max_seqs)}. Use the command line for larger sets.`);
      return;
    }
    if (longest > L.max_seq_len) {
      seqError(`The longest sequence is ${fmt(longest)} bp: the web app accepts up to ${fmt(L.max_seq_len)} bp. Use the command line for longer genomes.`);
      return;
    }
    seqError('');
    state.seqs = { name, text, ...parsed };
    $('file-name').textContent = name;
    $('file-stats').textContent = `${fmt(parsed.ids.length)} sequence${parsed.ids.length === 1 ? '' : 's'} · median ${fmt(parsed.median)} bp`;
    $('drop-empty').hidden = true;
    $('drop-full').hidden = false;
    const ref = $('reference_id');
    ref.length = 1;
    parsed.ids.forEach((id) => ref.add(new Option(id, id)));
    updateHint();
    if (!state.job) showIdle();
    if (state.engine && state.engine.warmup) state.engine.warmup();
  }

  function clearSeqs() {
    state.seqs = null;
    $('paste').value = '';
    $('file').value = '';
    $('drop-empty').hidden = false;
    $('drop-full').hidden = true;
    $('reference_id').length = 1;
    seqError('');
    updateHint();
    if (!state.job) showIdle();
  }

  function readFile(file) {
    if (!file) return;
    const maxBytes = (state.limits ? state.limits.max_upload_mb : 20) * 1024 * 1024;
    if (file.size > maxBytes) {
      seqError(`${file.name} is larger than ${state.limits ? state.limits.max_upload_mb : 20} MB. Use the command line for larger inputs.`);
      return;
    }
    file.text().then((text) => loadText(file.name, text));
  }

  // ---------- amplicon size hint and planned layout ----------
  const num = (id) => Number($(id).value);

  function maxAmplicons(G, overlap, lenMax) {
    let best = 0;
    for (let n = 1; n <= G; n += 1) {
      const L = Math.round((G + (n - 1) * overlap) / n);
      if (Math.floor((Math.floor(L * 0.8) - 1) / 2) < overlap + 2 * lenMax) break;
      best = n;
    }
    return best;
  }

  function plan() {
    if (!state.seqs) return null;
    const G = state.seqs.median;
    const n = Math.round(num('amplicons'));
    const overlap = num('overlap');
    if (!(n >= 1)) return { G, n, error: 'Enter how many amplicons you want.' };
    const most = maxAmplicons(G, overlap, num('primer_len_max'));
    if (n > most) return { G, n, error: `Too many for a ${fmt(G)} bp genome: at most ${fmt(most)} amplicons fit.` };
    const L = Math.round((G + (n - 1) * overlap) / n);
    return { G, n, L, overlap };
  }

  function updateHint() {
    const p = plan();
    const hint = $('amp-hint');
    hint.classList.toggle('error', Boolean(p && p.error));
    hint.classList.toggle('help', !(p && p.error));
    if (!p) hint.textContent = 'Add sequences to see the amplicon size.';
    else if (p.error) hint.textContent = p.error;
    else hint.textContent = `≈ ${fmt(p.L)} bp each · ${fmt(p.overlap)} bp overlap · ${p.n === 1 ? '1 pool' : '2 pools'}`;
  }

  // ---------- the genome map ----------
  function niceStep(G, maxTicks) {
    for (const s of [100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000, 50000]) {
      if (G / s <= maxTicks) return s;
    }
    return 100000;
  }
  const tickLabel = (v) => (v === 0 ? '0' : v >= 1000 ? `${+(v / 1000).toFixed(1)} kb` : `${v}`);

  function drawMap(items, G, lanes) {
    state.map = { items, G, lanes };
    $('map-figure').hidden = false; // measure the real width, not a hidden box
    const svg = $('map');
    const W = Math.max(280, Math.round(svg.parentElement.clientWidth || 800));
    const left = lanes.length > 1 || lanes[0] ? 64 : 0;
    const right = 8;
    const laneH = 30;
    const H = lanes.length * (laneH + 8) + 30;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.setAttribute('width', W);
    svg.setAttribute('height', H);
    const plotW = W - left - right;
    const x = (pos) => left + ((pos - 1) / G) * plotW;
    const NS = 'http://www.w3.org/2000/svg';
    const el = (tag, attrs, parent) => {
      const node = document.createElementNS(NS, tag);
      Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, v));
      (parent || svg).appendChild(node);
      return node;
    };
    svg.textContent = '';
    lanes.forEach((name, i) => {
      if (name) el('text', { x: 0, y: i * (laneH + 8) + laneH / 2 + 4, class: 'lane-label' }).textContent = name;
    });
    const axisY = lanes.length * (laneH + 8) + 2;
    el('line', { x1: left, x2: left + plotW, y1: axisY, y2: axisY, class: 'axis' });
    const step = niceStep(G, Math.max(3, Math.floor(plotW / 90)));
    for (let v = 0; v <= G; v += step) {
      const tx = left + (v / G) * plotW;
      el('line', { x1: tx, x2: tx, y1: axisY, y2: axisY + 5, class: 'tick' });
      el('text', { x: tx, y: axisY + 19, 'text-anchor': v === 0 ? 'start' : 'middle' }).textContent = tickLabel(v);
    }
    items.forEach((it) => {
      if (!(it.start > 0 && it.end >= it.start)) return;
      const lane = Math.max(0, lanes.length > 1 ? it.lane : 0);
      const g = el('g', {
        class: `amp${it.low ? ' low' : ''}${it.planned ? ' planned' : ''}`,
        tabindex: it.planned ? -1 : 0,
        'data-index': it.key,
        'aria-label': it.detail,
      });
      const bx = x(it.start);
      const bw = Math.max(3, x(it.end + 1) - bx);
      el('rect', { x: bx, y: lane * (laneH + 8), width: bw, height: laneH, rx: 3 }, g);
      el('title', {}, g).textContent = it.detail;
      if (bw > it.label.length * 7.5 + 10) {
        el('text', { x: bx + bw / 2, y: lane * (laneH + 8) + laneH / 2 + 4, 'text-anchor': 'middle' }, g).textContent = it.label;
      }
      if (!it.planned) {
        g.addEventListener('mouseenter', () => select(it.key));
        g.addEventListener('focus', () => select(it.key));
        g.addEventListener('click', () => select(it.key, true));
        g.addEventListener('keydown', (e) => {
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); select(it.key, true); }
        });
      }
    });
    $('map-figure').hidden = false;
  }

  function select(key, scroll) {
    const item = state.map && state.map.items.find((it) => it.key === key);
    if (!item) return;
    $('map-figure').classList.add('dim');
    document.querySelectorAll('#map .amp').forEach((g) => g.classList.toggle('sel', g.dataset.index === String(key)));
    document.querySelectorAll('#table tbody tr').forEach((tr) => tr.classList.toggle('sel', tr.dataset.index === String(key)));
    $('map-detail').textContent = item.detail;
    if (scroll) {
      const row = document.querySelector(`#table tbody tr[data-index="${key}"]`);
      if (row) row.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    }
  }

  window.addEventListener('resize', () => {
    if (state.map) drawMap(state.map.items, state.map.G, state.map.lanes);
  });

  // ---------- idle, progress and results views ----------
  function resetResults() {
    ['progress', 'figures', 'table-wrap', 'notes', 'downloads', 'steps'].forEach((id) => { $(id).hidden = true; });
    $('map-figure').hidden = true;
    $('map-figure').classList.remove('dim');
    $('map-detail').textContent = '';
    $('verdict').className = 'verdict';
    state.map = null;
  }

  function drawPlan() {
    const p = mode() === 'tile' ? plan() : null;
    if (!p || p.error) return null;
    const items = [];
    const room = p.G - p.L;
    for (let k = 0; k < p.n; k += 1) {
      const start = p.n === 1 ? 1 : Math.round(1 + (k * room) / (p.n - 1));
      items.push({ key: k + 1, label: String(k + 1), lane: k % 2, start, end: Math.min(p.G, start + p.L - 1), planned: true, detail: `Planned amplicon ${k + 1}` });
    }
    drawMap(items, p.G, p.n === 1 ? ['Pool 1'] : ['Pool 1', 'Pool 2']);
    $('map-detail').textContent = 'Grey bars show the planned spacing. The run moves each one onto primer sites that match your sequences.';
    return p;
  }

  const OUTCOME = {
    tile: 'When it finishes you get both primers of every amplicon, how many of your sequences each amplicon covers and why it misses any, and the files primers.tsv (names and sequences, ready to order), amplicons.tsv (every detail), scheme.bed (ARTIC primer scheme) and report.txt.',
    regions: 'When it finishes you get the best primer pairs in each region, ranked by how many of your sequences they amplify, with primers.tsv (the best pair of each region, ready to order), pairs.tsv (every pair with its details) and report.txt.',
    check: 'Add your sequences and paste your primer pairs. Each pair is placed on the alignment and scored on every sequence, with the reason for each miss.',
  };

  function showIdle() {
    resetResults();
    const m = mode();
    const verdict = $('verdict');
    $('empty').hidden = false;
    $('empty').querySelector('p').textContent = OUTCOME[m];
    $('steps').hidden = m === 'check';
    const p = drawPlan();
    if (p) {
      verdict.className = 'verdict plan';
      verdict.textContent = `Planned: ${p.n} amplicon${p.n === 1 ? '' : 's'} of about ${fmt(p.L)} bp across ${fmt(p.G)} bp.`;
    } else if (m === 'tile') {
      verdict.textContent = 'Your tiling will appear here.';
    } else if (m === 'regions') {
      verdict.textContent = 'The best pairs for each region will appear here.';
    } else {
      verdict.textContent = 'Coverage for each of your pairs will appear here.';
    }
  }

  function showProgress(ev) {
    resetResults();
    $('empty').hidden = true;
    drawPlan();
    $('progress').hidden = false;
    $('cancel').hidden = false;
    $('verdict').className = 'verdict plan';
    $('verdict').textContent = state.engine && state.engine.name === 'browser'
      ? 'Working on it on your computer. Keep this tab open.'
      : 'Working on it. You can leave this page open or come back to this address later.';
    let label;
    if (ev.phase === 'loading') label = `Preparing the analysis tools: ${ev.stage || 'starting'}…`;
    else if (ev.phase === 'queued') label = ev.queuePosition > 1 ? `Waiting for a free slot (${ev.queuePosition - 1} job${ev.queuePosition > 2 ? 's' : ''} ahead)` : 'Waiting for a free slot…';
    else label = STAGES[ev.stage] || 'Starting…';
    $('progress-label').textContent = label;
    const value = Math.round((ev.progress || 0) * 100);
    $('bar').setAttribute('aria-valuenow', value);
    $('bar-fill').style.width = `${Math.max(value, ev.phase === 'queued' ? 0 : 3)}%`;
  }

  function reasons(notCovered) {
    const counts = {};
    (notCovered.match(/\(([^)]*)\)/g) || []).forEach((m) => {
      m.slice(1, -1).split(';').forEach((r) => { const k = r.trim(); if (k) counts[k] = (counts[k] || 0) + 1; });
    });
    return Object.entries(counts).sort((a, b) => b[1] - a[1]);
  }

  // "21 nt · Tm 60.3 °C (59.4–61.2) · GC 47.6 % · 1 degenerate base (2 variants)"
  function primerFacts(r, side) {
    const seq = r[`${side}_seq`] || '';
    const parts = [`${r[`${side}_len`] || seq.length} nt`];
    const tm = r[`${side}_tm`];
    const range = r[`${side}_tm_range`];
    if (tm) parts.push(`Tm ${tm} °C${range ? ` (${range.replace('-', '–')})` : ''}`);
    if (r[`${side}_gc`]) parts.push(`GC ${r[`${side}_gc`]} %`);
    const nd = +(r[`${side}_degenerate_bases`] ?? [...seq].filter((b) => !'ACGT'.includes(b)).length);
    parts.push(nd
      ? `${nd} degenerate base${nd === 1 ? '' : 's'} (${fmt(+r[`${side}_variants`] || 0)} variants)`
      : 'no degenerate bases');
    return parts.join(' · ');
  }

  function primerCell(r, owner) {
    const td = document.createElement('td');
    td.className = 'wide primers';
    [['F', 'forward'], ['R', 'reverse']].forEach(([tag, what]) => {
      const seq = r[`${what}_seq`];
      const line = document.createElement('span');
      line.className = 'seq';
      const mark = document.createElement('span');
      mark.className = 'strand';
      mark.textContent = tag;
      const code = document.createElement('code');
      code.textContent = seq;
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'copy';
      btn.textContent = 'Copy';
      btn.setAttribute('aria-label', `Copy the ${what} primer of ${owner}`);
      btn.addEventListener('click', () => copy(seq, btn));
      line.append(mark, code, btn);
      const facts = document.createElement('span');
      facts.className = 'facts';
      facts.textContent = primerFacts(r, what);
      td.append(line, facts);
    });
    return td;
  }

  function leadCell(title, sub) {
    const td = document.createElement('td');
    td.className = 'lead';
    td.textContent = title;
    if (sub) {
      const small = document.createElement('span');
      small.className = 'sub';
      small.textContent = sub;
      td.append(small);
    }
    return td;
  }

  function copy(text, btn) {
    const done = () => {
      const old = btn.textContent;
      btn.textContent = 'Copied';
      btn.classList.add('done');
      setTimeout(() => { btn.textContent = old; btn.classList.remove('done'); }, 1400);
    };
    if (navigator.clipboard) navigator.clipboard.writeText(text).then(done, () => {});
  }

  function cell(text, cls, label) {
    const td = document.createElement('td');
    if (cls) td.className = cls;
    if (label) td.dataset.label = label;
    td.textContent = text;
    return td;
  }

  function figures(pairs) {
    const dl = $('figures');
    dl.textContent = '';
    pairs.forEach(([term, value, unit]) => {
      const div = document.createElement('div');
      const dt = document.createElement('dt');
      dt.textContent = term;
      const dd = document.createElement('dd');
      dd.textContent = value;
      if (unit) {
        const small = document.createElement('small');
        small.textContent = unit;
        dd.append(small);
      }
      div.append(dt, dd);
      dl.append(div);
    });
    dl.hidden = false;
  }

  function span(r) {
    return r.amplicon_start ? `${fmt(r.amplicon_start)}–${fmt(r.amplicon_end)}` : '—';
  }

  function renderTable(head, rows) {
    const table = $('table');
    table.textContent = '';
    const thead = table.createTHead();
    const tr = thead.insertRow();
    head.forEach(([text, cls]) => {
      const th = document.createElement('th');
      th.textContent = text;
      if (cls) th.className = cls;
      tr.append(th);
    });
    const tbody = table.createTBody();
    rows.forEach((row) => tbody.append(row));
    $('table-wrap').hidden = false;
  }

  function renderDownloads(outcome, orderList) {
    const row = $('dl-row');
    row.textContent = '';
    const link = ({ name, href }, label) => {
      const a = document.createElement('a');
      a.className = 'btn btn-quiet';
      a.href = href;
      a.download = name;
      a.textContent = label || name;
      row.append(a);
    };
    outcome.downloads.forEach((d) => link(d));
    if (outcome.zip) link(outcome.zip, 'All files (.zip)');
    if (orderList) {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'btn btn-quiet';
      btn.textContent = 'Copy order list';
      btn.addEventListener('click', () => copy(orderList, btn));
      row.append(btn);
    }
    $('downloads').hidden = false;
  }

  function renderNotes(warnings) {
    const ul = $('notes');
    ul.textContent = '';
    (warnings || []).filter((w) => !/^Pairs? |^Amplicons below|^Regions without/.test(w)).forEach((w) => {
      const li = document.createElement('li');
      li.textContent = w;
      ul.append(li);
    });
    ul.hidden = !ul.children.length;
  }

  function renderTile(st, res) {
    const rows = res.rows;
    const G = res.reference_length || Math.max(...rows.map((r) => +r.amplicon_end || 0));
    if (!state.seqs) { $('amplicons').value = rows.length; updateHint(); }
    const low = rows.filter((r) => r.below_threshold === 'True');
    const worst = Math.min(...rows.map((r) => +r.coverage_pct));
    const total = rows.length ? +rows[0].n_total : 0;
    const verdict = $('verdict');
    if (!low.length) {
      verdict.textContent = worst === 100
        ? `All ${rows.length} amplicons cover every sequence`
        : `All ${rows.length} amplicons meet the coverage target; the lowest covers ${pct(worst)}`;
    } else {
      const lowest = low.reduce((a, b) => (+a.coverage_pct <= +b.coverage_pct ? a : b));
      const why = reasons(lowest.not_covered)[0];
      verdict.innerHTML = '';
      const strong = document.createElement('span');
      strong.className = 'miss';
      strong.textContent = `${low.length} of ${rows.length} amplicons fall short`;
      verdict.append(strong, ` of the coverage target: ${low.map((r) => r.amplicon).join(', ')}.${why ? ` Most misses in amplicon ${lowest.amplicon}: ${why[0]}.` : ''}`);
    }
    figures([
      ['Amplicons', fmt(rows.length), rows.length > 1 ? ' in 2 pools' : ' in 1 pool'],
      ['Lowest coverage', pct(worst).replace(' %', ''), ' %'],
      ['Sequences', fmt(total)],
      ['Genome', fmt(G), ' bp'],
    ]);
    const items = rows.map((r) => ({
      key: +r.amplicon, label: r.amplicon, lane: +r.pool - 1, start: +r.amplicon_start, end: +r.amplicon_end,
      low: r.below_threshold === 'True',
      detail: `Amplicon ${r.amplicon} · pool ${r.pool} · ${span(r)} (${fmt(r.amplicon_len || 0)} bp) · covers ${r.n_covered} of ${r.n_total} sequences (${pct(r.coverage_pct)})`,
    }));
    drawMap(items, G, rows.length > 1 ? ['Pool 1', 'Pool 2'] : ['Pool 1']);
    $('map-detail').textContent = 'Point at an amplicon to see its position and coverage.';
    renderTable(
      [['Amplicon'], ['Primers (5′→3′)'], ['Position', 'num'], ['Length', 'num'], ['Coverage', 'num']],
      rows.map((r) => {
        const tr = document.createElement('tr');
        tr.dataset.index = r.amplicon;
        if (r.below_threshold === 'True') tr.className = 'low';
        tr.append(
          leadCell(r.amplicon, `pool ${r.pool}`),
          primerCell(r, `amplicon ${r.amplicon}`),
          cell(span(r), 'num pos'),
          cell(r.amplicon_len ? `${fmt(r.amplicon_len)} bp` : '—', 'num len'),
          cell(pct(r.coverage_pct), 'num cov'),
        );
        tr.addEventListener('mouseenter', () => select(+r.amplicon));
        return tr;
      }),
    );
    const order = rows.flatMap((r) => [`${r.forward_name}\t${r.forward_seq}`, `${r.reverse_name}\t${r.reverse_seq}`]).join('\n');
    renderDownloads(st, order);
  }

  function renderPairs(st, res, m) {
    const rows = res.rows;
    const G = res.reference_length || Math.max(1, ...rows.map((r) => +r.amplicon_end || 0));
    const nameOf = (r) => (m === 'check' ? r.name : `${r.region} #${r.rank}`);
    const best = m === 'check' ? rows : rows.filter((r) => r.rank === '1');
    const low = best.filter((r) => r.below_threshold === 'True');
    const noun = m === 'check' ? 'pairs' : 'regions';
    const verdict = $('verdict');
    verdict.textContent = '';
    if (!best.length) verdict.textContent = 'No pairs were found.';
    else if (!low.length) verdict.textContent = `All ${best.length} ${noun} meet the coverage target`;
    else {
      const strong = document.createElement('span');
      strong.className = 'miss';
      strong.textContent = `${low.length} of ${best.length} ${noun} fall short`;
      verdict.append(strong, ` of the coverage target: ${low.map((r) => (m === 'check' ? r.name : r.region)).join(', ')}.`);
    }
    const worst = best.length ? Math.min(...best.map((r) => +r.coverage_pct)) : 0;
    figures([
      [m === 'check' ? 'Pairs' : 'Regions', fmt(best.length)],
      ['Lowest coverage', pct(worst).replace(' %', ''), ' %'],
      ['Sequences', fmt(best.length ? best[0].n_total : 0)],
      ['Reference', fmt(G), ' bp'],
    ]);
    const items = best.map((r, i) => ({
      key: i + 1, label: m === 'check' ? r.name : r.region, lane: 0, start: +r.amplicon_start, end: +r.amplicon_end,
      low: r.below_threshold === 'True',
      detail: `${nameOf(r)} · ${span(r)} · covers ${r.n_covered} of ${r.n_total} sequences (${pct(r.coverage_pct)})${r.note ? ` · ${r.note}` : ''}`,
    }));
    drawMap(items, G, [m === 'check' ? 'Pairs' : 'Best pairs']);
    $('map-detail').textContent = 'Point at a pair to see its position and coverage.';
    const keyOf = new Map(best.map((r, i) => [r, i + 1]));
    const head = [['Name'], ['Primers (5′→3′)'], ['Position', 'num'], ['Coverage', 'num']];
    renderTable(head, rows.map((r) => {
      const tr = document.createElement('tr');
      if (keyOf.has(r)) tr.dataset.index = keyOf.get(r);
      if (r.below_threshold === 'True') tr.className = 'low';
      tr.append(
        leadCell(nameOf(r), m === 'check' ? r.note : `rank ${r.rank}`),
        primerCell(r, nameOf(r)),
        cell(span(r), 'num pos'),
        cell(pct(r.coverage_pct), 'num cov'),
      );
      if (keyOf.has(r)) tr.addEventListener('mouseenter', () => select(keyOf.get(r)));
      return tr;
    }));
    renderDownloads(st, null);
  }

  function renderResult(st) {
    resetResults();
    $('empty').hidden = true;
    const res = st.result || {};
    if (st.status === 'failed') {
      const verdict = $('verdict');
      verdict.className = 'verdict failed';
      verdict.textContent = `The run stopped: ${st.error || 'something went wrong.'}`;
      renderNotes(res.warnings);
      showRetry();
      return;
    }
    if (res.mode === 'tile') renderTile(st, res);
    else renderPairs(st, res, res.mode);
    renderNotes(res.warnings);
  }

  function showRetry() {
    const ul = $('notes');
    const li = document.createElement('li');
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-quiet';
    btn.textContent = 'Retry';
    btn.addEventListener('click', () => { if (state.lastJob) startRun(state.lastJob); });
    li.append(btn);
    ul.append(li);
    ul.hidden = false;
  }

  function showFailure(message) {
    resetResults();
    $('empty').hidden = true;
    const verdict = $('verdict');
    verdict.className = 'verdict failed';
    verdict.textContent = message;
    $('notes').textContent = '';
    showRetry();
  }

  // ---------- running a job ----------
  function settings(m) {
    const out = {};
    if (m === 'tile') out.amplicons = Math.round(num('amplicons'));
    if (m === 'regions') {
      out.regions = $('regions').value.split(/\r?\n/).map((s) => s.trim()).filter(Boolean);
      out.amplicon_min = num('amplicon_min');
      out.amplicon_max = num('amplicon_max');
      if ($('allow_degenerate').checked) out.allow_degenerate = true;
    }
    if (m === 'tile' && !$('degenerate').checked) out.degenerate = false;
    document.querySelectorAll('[data-setting]').forEach((input) => {
      const holder = input.closest('[data-modes]');
      if (!holder || !holder.dataset.modes.split(' ').includes(m)) return;
      const key = input.dataset.setting;
      if (key in out) return;
      const changed = input.tagName === 'SELECT' ? input.value !== '' && input.selectedIndex !== 0 : input.value !== input.defaultValue;
      if (!changed || input.value === '') return;
      out[key] = input.type === 'number' ? Number(input.value) : input.value;
    });
    return out;
  }

  function formError(message) {
    $('form-error').textContent = message || '';
    $('form-error').hidden = !message;
  }

  async function submit(event) {
    event.preventDefault();
    formError('');
    const m = mode();
    if (!state.seqs) {
      seqError('Add your sequences first: drop a FASTA file, paste sequences, or try the example.');
      $('paste').focus();
      return;
    }
    if (m === 'tile') {
      const p = plan();
      if (p && p.error) { formError(p.error); $('amplicons').focus(); return; }
    }
    if (m === 'check' && !$('primers').value.trim()) {
      formError('Paste your primer pairs first.');
      $('primers').focus();
      return;
    }
    let primers = '';
    if (m === 'check') {
      primers = `${$('primers').value.split(/\r?\n/).map((l) => l.trim()).filter(Boolean)
        .map((l) => l.split(/\s+/).join('\t')).join('\n')}\n`;
    }
    if (state.engine.name === 'browser') {
      const bases = state.seqs.lens.reduce((a, b) => a + b, 0);
      if (bases > LARGE_INPUT && !window.confirm('This may take several minutes and use a lot of memory on your computer. For large sets, the command line is faster.\n\nRun it here anyway?')) return;
    }
    startRun({ mode: m, settings: settings(m), fasta: state.seqs.text, primers });
  }

  function onEvent(ev) {
    if (ev.jobId && state.engine.name === 'server'
        && new URLSearchParams(location.search).get('job') !== ev.jobId) {
      history.replaceState(null, '', `?job=${ev.jobId}#${state.page}`);
    }
    showProgress(ev);
  }

  async function startRun(job, resumeId) {
    state.lastJob = job || state.lastJob;
    state.job = true;
    $('run').disabled = true;
    formError('');
    showProgress({ phase: state.engine.name === 'browser' ? 'loading' : 'queued', progress: 0 });
    const handle = resumeId ? state.engine.resume(resumeId, onEvent) : state.engine.run(job, onEvent);
    state.handle = handle;
    try {
      const outcome = await handle.promise;
      state.hasResults = state.engine.name === 'browser';
      renderResult(outcome);
      $('results').focus({ preventScroll: true });
    } catch (err) {
      if (err.cancelled) {
        history.replaceState(null, '', `${location.pathname}#${state.page}`);
        showIdle();
        $('verdict').className = 'verdict plan';
        $('verdict').textContent = 'Cancelled. Change the settings or run it again.';
      } else if (err.invalid) {
        showIdle();
        formError(err.message);
      } else {
        showFailure(err.message);
      }
    } finally {
      state.job = null;
      state.handle = null;
      $('run').disabled = false;
      $('cancel').hidden = true;
    }
  }

  // ---------- wiring ----------
  function setPage(page) {
    if (state.job) return;
    state.page = page;
    applyMode();
    updateHint();
  }

  document.querySelectorAll('.modes a').forEach((a) => {
    a.addEventListener('click', (e) => {
      e.preventDefault();
      history.replaceState(null, '', `#${a.dataset.mode}`);
      setPage(a.dataset.mode);
    });
  });
  $('use-regions').addEventListener('change', applyMode);
  $('form').addEventListener('submit', submit);
  $('file').addEventListener('change', (e) => readFile(e.target.files[0]));
  $('clear').addEventListener('click', clearSeqs);
  $('example').addEventListener('click', async () => {
    const r = await fetch(state.engine && state.engine.name === 'server' ? '/api/example' : 'example.fasta');
    if (r.ok) loadText('Example: 24 synthetic genomes', await r.text());
  });
  let pasteTimer;
  $('paste').addEventListener('input', () => {
    clearTimeout(pasteTimer);
    pasteTimer = setTimeout(() => {
      const text = $('paste').value;
      if (text.trim()) loadText('Pasted sequences', text);
    }, 500);
  });
  ['amplicons', 'overlap', 'primer_len_max'].forEach((id) => $(id).addEventListener('input', () => {
    updateHint();
    if (!state.job) showIdle();
  }));
  const drop = $('drop');
  ['dragenter', 'dragover'].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.add('over'); }));
  ['dragleave', 'drop'].forEach((t) => drop.addEventListener(t, () => drop.classList.remove('over')));
  drop.addEventListener('drop', (e) => { e.preventDefault(); readFile(e.dataTransfer.files[0]); });
  $('results').tabIndex = -1;

  $('cancel').addEventListener('click', () => { if (state.handle) state.handle.cancel(); });
  window.addEventListener('beforeunload', (e) => {
    if (state.job || state.hasResults) {
      e.preventDefault();
      e.returnValue = '';
    }
  });

  function applyLimits(limits) {
    state.limits = limits;
    document.querySelectorAll('[data-limit]').forEach((el) => { el.textContent = fmt(limits[el.dataset.limit]); });
    $('amplicons').max = limits.max_amplicons;
    Array.from($('aligner').options).forEach((o) => {
      o.disabled = !(limits.aligners || []).includes(o.value);
      if (o.disabled && !/not available/.test(o.textContent)) o.textContent += ' (not available here)';
    });
  }

  state.page = location.hash === '#check' ? 'check' : 'design';
  applyMode();
  $('run').disabled = true;
  window.VhirEngines.pickEngine().then(async (engine) => {
    state.engine = engine;
    applyLimits(await engine.limits());
    const note = $('privacy');
    if (engine.name === 'browser') {
      note.textContent = 'Runs on your computer — your sequences are not uploaded. The first run loads the analysis tools (about 20 MB); your browser keeps them for next time.';
    } else {
      note.textContent = `Files are deleted ${fmt(state.limits.retention_hours)} hours after your job finishes. For private or larger data, use the command line.`;
    }
    $('run').disabled = false;
    const jobId = new URLSearchParams(location.search).get('job');
    if (engine.name === 'server' && jobId && /^[0-9a-f]{32}$/.test(jobId)) startRun(null, jobId);
  });
})();
