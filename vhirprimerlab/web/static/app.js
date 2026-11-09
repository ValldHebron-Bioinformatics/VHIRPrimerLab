/* VHIRPrimerLab web app: read sequences in the browser, start a job, follow it, draw the result. */
import { buildAnchorMap } from './anchors.js';
import { fetchAnnotation, parseGenbank } from './annotation.js';
import { ampliconCount, blockFromFeature, mergeBlocks, smallestSpan, snapPosition } from './blocks.js';

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
    'checking primer dimers': 'Checking primer dimers…',
    writing: 'Writing the files…',
    done: 'Finishing…',
  };
  const FILES = {
    tile: ['primers.tsv', 'amplicons.tsv', 'scheme.bed', 'report.txt', 'alignment.fasta'],
    regions: ['primers.tsv', 'pairs.tsv', 'report.txt', 'alignment.fasta'],
    check: ['check.tsv', 'report.txt', 'alignment.fasta'],
  };

  const state = { page: 'design', seqs: null, limits: null, job: null, map: null, result: null,
    engine: null, handle: null, hasResults: false, lastJob: null,
    genome: { ann: null, map: null, refId: null, refLen: 0, blocks: [], edges: [], note: '',
      drag: null, pending: null } };

  // ---------- mode ----------
  const PAGES = ['design', 'genes', 'check'];
  const advanced = () => document.querySelector('input[name="level"]:checked').value === 'advanced';
  const mode = () => (state.page === 'check' ? 'check'
    : state.page === 'genes' ? 'genes'
    : advanced() && $('use-regions').checked ? 'regions' : 'tile');

  const COPY = {
    tile: {
      title: 'Design primers',
      lede: 'Drop your genomes, choose how many amplicons you want, and get overlapping primer pairs in two pools, checked against every sequence.',
      run: 'Design primers',
      tab: 'tile a genome with primers',
    },
    genes: {
      title: 'Pick by gene',
      lede: 'Drop your genomes, pick the organism, and click the genes you want: each one is tiled with overlapping primer pairs, checked against every sequence.',
      run: 'Design primers',
      tab: 'pick primers by gene',
    },
    check: {
      title: 'Check my primers',
      lede: 'Paste primer pairs you already have and see how many of your sequences each pair amplifies, with the same matching rules as design.',
      run: 'Check primers',
      tab: 'check primers',
    },
  };

  // ---------- basic / advanced ----------
  const LEVEL_KEY = 'vhirprimerlab.level';

  function applyLevel() {
    $('advanced').hidden = !advanced();
    $('basic-note').hidden = advanced();
    try { localStorage.setItem(LEVEL_KEY, advanced() ? 'advanced' : 'basic'); } catch (e) { /* private window */ }
    drawRule();
    applyMode();
    updateHint();
  }

  // One box per base of a 20-base primer, 5' to 3'; the last N must match exactly.
  function drawRule() {
    const strip = $('rule-strip');
    const mm = Math.max(0, Math.round(num('max_mismatches')) || 0);
    const raw = Math.round(num('three_prime_window'));
    const exact = Math.min(20, Math.max(1, raw || 1));
    strip.textContent = '';
    const end = (t) => { const e = document.createElement('span'); e.className = 'end'; e.textContent = t; return e; };
    strip.append(end('5′'));
    for (let i = 0; i < 20; i += 1) {
      const b = document.createElement('span');
      b.className = i >= 20 - exact ? 'b exact' : 'b';
      strip.append(b);
    }
    strip.append(end('3′'));
    const free = 20 - exact;
    const where = free
      ? `anywhere in the first ${free} bases of a 20-base primer`
      : 'nowhere in a 20-base primer';
    $('rule-text').textContent = `A primer binds a sequence when it has at most ${mm} mismatch${mm === 1 ? '' : 'es'}, ${where}. The last ${exact} base${exact === 1 ? '' : 's'} at the 3′ end (filled) must match exactly: the polymerase starts there.`;
  }

  function applyMode() {
    const m = mode();
    document.querySelectorAll('[data-modes]').forEach((el) => {
      el.hidden = !el.dataset.modes.split(' ').includes(m);
    });
    document.querySelectorAll('.modes a').forEach((a) => {
      if (a.dataset.mode === state.page) a.setAttribute('aria-current', 'page');
      else a.removeAttribute('aria-current');
    });
    const copy = COPY[m] || COPY.tile;   // regions is the design page with its own fields
    $('title').textContent = copy.title;
    $('lede').textContent = copy.lede;
    $('run').textContent = copy.run;
    document.title = `VHIRPrimerLab: ${copy.tab}`;
    if (!state.job) showIdle();
    if (m === 'genes') refreshGenome();
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
    refreshGenome(true);   // the genes diagram is drawn over the reference, which just changed
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
    refreshGenome(true);
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

  // ---------- pick by gene ----------
  // The annotation comes from NCBI (annotation.js), is placed on the user's own reference with
  // shared k-mers (anchors.js), and the spans the user clicks together are kept by blocks.js.
  // Only the drawing, the dragging and the wiring live here.

  // The records asked for most often; any other accession goes in the box below them.
  const ORGANISMS = [
    { label: 'Human metapneumovirus', accession: 'NC_039199.1' },
    { label: 'Respiratory syncytial virus A', accession: 'NC_038235.1' },
    { label: 'Respiratory syncytial virus B', accession: 'NC_001781.1' },
    { label: 'SARS-CoV-2', accession: 'NC_045512.2' },
    { label: 'Influenza A (H1N1) segment 4, HA', accession: 'NC_026433.1' },
    { label: 'Influenza A (H1N1) segment 6, NA', accession: 'NC_026434.1' },
  ];
  const LANES = [
    { label: 'Genes', cls: 'lane-gene' },
    { label: 'CDS', cls: 'lane-cds' },
    { label: 'Other', cls: 'lane-other' },
  ];
  const laneOf = (f) => (f.kind === 'gene' ? 0 : f.kind === 'CDS' ? 1 : 2);
  // A span's name becomes the prefix of every primer name in the order sheet, and the engine
  // accepts only letters, digits, _ . - in it. A feature named after its product has spaces.
  const spanName = (name) => name.replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '') || 'span';
  // The engine refuses two spans with the same name, which two features of one name would give.
  const freeName = (name, taken) => {
    if (!taken.has(name)) return name;
    for (let n = 2; ; n += 1) if (!taken.has(`${name}.${n}`)) return `${name}.${n}`;
  };
  const GRID = 10;                 // a drag that lands near no feature edge snaps to this
  const minSpan = () => smallestSpan(num('primer_len_max') || 25);
  let plot = null;                 // the last draw's geometry, so a drag can read bases off pixels

  function genomeError(message) {
    $('organism-error').textContent = message || '';
    $('organism-error').hidden = !message;
  }

  // The chosen reference without its gaps: the coordinates the engine and the server both use.
  function referenceBases() {
    if (!state.seqs) return '';
    const want = $('reference_id').value;
    const parts = [];
    let keep = false;
    for (const line of state.seqs.text.split(/\r?\n/)) {
      if (line.startsWith('>')) {
        if (keep) break;
        keep = !want || (line.slice(1).split(/\s+/)[0] || '') === want;   // no choice: the first
      } else if (keep) {
        parts.push(line.replace(/[^A-Za-z]/g, ''));
      }
    }
    return parts.join('').toUpperCase();
  }

  // A trusted map puts the diagram and the spans in the user's reference coordinates. An
  // untrusted one would invent positions, so both fall back to the record's own and the note
  // says so: whatever the case, the ruler and the spans are in the same coordinates.
  function axis() {
    const g = state.genome;
    if (g.map && g.map.trusted) return { convert: g.map.convert, G: g.refLen, record: false };
    return { convert: (pos) => pos, G: Math.max(1, (g.ann && g.ann.length) || g.refLen), record: true };
  }

  function featureDetail(f, start, end, off) {
    const parts = [f.name];
    if (f.product && f.product !== f.name) parts.push(f.product);
    parts.push(f.kind, `${fmt(start)}–${fmt(end)} (${fmt(end - start + 1)} bp)`, `${f.strand} strand`);
    if (off) parts.push('outside the part of the record your sequences cover');
    return parts.join(' · ');
  }

  function drawGenome(ann, map, blocks) {
    const { convert, G, record } = axis();
    $('genome-figure').hidden = false;   // measure the real width, not a hidden box
    const svg = $('genome');
    const W = Math.max(280, Math.round(svg.parentElement.clientWidth || 800));
    const left = 64;
    const right = 8;
    const band = 28;                     // room above the lanes for the block brackets
    const laneH = 24;
    const laneGap = 8;
    const rows = LANES.flatMap((lane, i) => (ann.features.some((f) => laneOf(f) === i) ? [i] : []));
    const rowOf = new Map(rows.map((i, r) => [i, r]));
    const H = band + Math.max(1, rows.length) * (laneH + laneGap) + 30;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.setAttribute('width', W);
    svg.setAttribute('height', H);
    const plotW = W - left - right;
    const x = (pos) => left + ((Math.min(G, Math.max(1, pos)) - 1) / G) * plotW;
    plot = { W, left, plotW, G };
    const NS = 'http://www.w3.org/2000/svg';
    const el = (tag, attrs, parent) => {
      const node = document.createElementNS(NS, tag);
      Object.entries(attrs).forEach(([k, v]) => node.setAttribute(k, v));
      (parent || svg).appendChild(node);
      return node;
    };
    svg.textContent = '';
    rows.forEach((i, r) => {
      el('text', { x: 0, y: band + r * (laneH + laneGap) + laneH / 2 + 4, class: 'lane-label' })
        .textContent = LANES[i].label;
    });
    const axisY = band + Math.max(1, rows.length) * (laneH + laneGap) + 2;
    el('line', { x1: left, x2: left + plotW, y1: axisY, y2: axisY, class: 'axis' });
    const step = niceStep(G, Math.max(3, Math.floor(plotW / 90)));
    for (let v = 0; v <= G; v += step) {
      const tx = left + (v / G) * plotW;
      el('line', { x1: tx, x2: tx, y1: axisY, y2: axisY + 5, class: 'tick' });
      el('text', { x: tx, y: axisY + 19, 'text-anchor': v === 0 ? 'start' : 'middle' }).textContent = tickLabel(v);
    }
    ann.features.forEach((f, i) => {
      const y = band + rowOf.get(laneOf(f)) * (laneH + laneGap);
      // Greying only means something under a trusted map; under the record's own ruler every
      // feature is equally unplaced, and the note already says so.
      const off = !record && (f.end < map.mappedRange[0] || f.start > map.mappedRange[1]);
      const start = convert(f.start);
      const end = convert(f.end);
      const detail = featureDetail(f, start, end, off);
      const g = el('g', {
        class: `feat ${LANES[laneOf(f)].cls}${off ? ' off' : ''}`,
        tabindex: off ? -1 : 0,
        'data-feature': i,
        'aria-label': detail,
      });
      const bx = x(start);
      const bw = Math.max(3, x(end + 1) - bx);
      if (f.intervals.length > 1) {   // a spliced feature: its pieces joined by a thin line
        el('line', { x1: bx, x2: bx + bw, y1: y + laneH / 2, y2: y + laneH / 2, class: 'join' }, g);
      }
      f.intervals.forEach(([a, b]) => {
        const ix = x(convert(a));
        el('rect', { x: ix, y, width: Math.max(3, x(convert(b) + 1) - ix), height: laneH, rx: 3 }, g);
      });
      el('title', {}, g).textContent = detail;
      if (bw > f.name.length * 7.5 + 10) {
        el('text', { x: bx + bw / 2, y: y + laneH / 2 + 4, 'text-anchor': 'middle' }, g).textContent = f.name;
      }
      if (off) return;
      const show = () => { $('genome-detail').textContent = detail; };
      g.addEventListener('mouseenter', show);
      g.addEventListener('focus', show);
      g.addEventListener('click', () => addBlock(f));
      g.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); addBlock(f); }
      });
    });
    blocks.forEach((b, i) => {
      const bx = x(b.start);
      const ex = x(Math.min(G, b.end + 1));
      const g = el('g', { class: 'block', 'data-block': i });
      el('path', { d: `M${bx} ${band - 6}V${band - 14}H${ex}V${band - 6}`, class: 'bracket' }, g);
      el('text', { x: (bx + ex) / 2, y: band - 18, 'text-anchor': 'middle' }, g).textContent = b.name;
      el('title', {}, g).textContent = `${b.name}: ${fmt(b.start)}–${fmt(b.end)}`;
      g.addEventListener('mouseenter', () => selectBlock(i));
      [['start', bx, b.start], ['end', ex, b.end]].forEach(([edge, hx, at]) => {
        const handle = el('rect', {
          class: 'handle', x: hx - 5, y: band - 18, width: 10, height: 16, rx: 2,
          tabindex: 0, 'data-block': i, 'data-edge': edge,
          'aria-label': `Move the ${edge} of ${b.name}, now ${fmt(at)}. Left and right arrows move it.`,
        }, g);
        handle.addEventListener('pointerdown', (e) => startDrag(e, i, edge));
        handle.addEventListener('keydown', (e) => nudge(e, i, edge));
        handle.addEventListener('focus', () => selectBlock(i));
      });
    });
  }

  function selectBlock(index) {
    const key = String(index);
    document.querySelectorAll('#genome .block').forEach((g) => g.classList.toggle('sel', g.dataset.block === key));
    document.querySelectorAll('#blocks li').forEach((li) => li.classList.toggle('sel', li.dataset.block === key));
    const b = state.genome.blocks[index];
    if (b) {
      $('genome-detail').textContent = `${b.name} · ${fmt(b.start)}–${fmt(b.end)} · `
        + `${fmt(b.end - b.start + 1)} bp. Drag either handle to adjust it.`;
    }
  }

  function redrawGenome() {
    const g = state.genome;
    if (!g.ann || !g.map) return;
    drawGenome(g.ann, g.map, g.blocks);
    renderBlocks();
  }

  // Spans that touch become one. The joined name becomes the prefix of every primer name in the
  // order sheet, so the note has to say which spans were joined.
  function commitBlocks(blocks) {
    const g = state.genome;
    const merged = mergeBlocks(blocks);
    const joined = merged
      .filter((m) => blocks.filter((b) => b.start >= m.start && b.end <= m.end).length > 1)
      .map((m) => m.name);
    g.note = joined.length ? `Spans that touch become one: ${joined.join(', ')}.` : '';
    g.blocks = merged;
    redrawGenome();
  }

  function addBlock(feature) {
    const g = state.genome;
    const name = freeName(spanName(feature.name), new Set(g.blocks.map((b) => b.name)));
    const block = blockFromFeature({ ...feature, name }, {
      convert: axis().convert, limit: g.refLen, minSpan: minSpan(),
    });
    commitBlocks([...g.blocks, block]);
    selectBlock(g.blocks.findIndex((b) => b.start <= block.start && b.end >= block.end));
  }

  function removeBlock(index) {
    const g = state.genome;
    g.blocks.splice(index, 1);
    g.note = '';
    redrawGenome();
  }

  function renderBlocks() {
    const g = state.genome;
    const ul = $('blocks');
    const ampLen = Math.round(num('amplicon_len'));
    const overlap = Math.round(num('gene_overlap'));
    const sane = ampLen > 0 && overlap >= 0 && ampLen > overlap;
    let total = 0;
    ul.textContent = '';
    g.blocks.forEach((b, i) => {
      const n = sane ? ampliconCount(b, ampLen, overlap) : 0;
      total += n;
      const li = document.createElement('li');
      li.dataset.block = String(i);
      const name = document.createElement('span');
      name.className = 'b-name';
      name.textContent = b.name;
      const where = document.createElement('span');
      where.className = 'b-span';
      where.textContent = `${fmt(b.start)}–${fmt(b.end)} · ${fmt(b.end - b.start + 1)} bp`;
      const count = document.createElement('span');
      count.className = 'b-count';
      count.textContent = sane ? `${fmt(n)} amplicon${n === 1 ? '' : 's'}` : '—';
      const remove = document.createElement('button');
      remove.type = 'button';
      remove.className = 'b-remove';
      remove.textContent = 'Remove';
      remove.setAttribute('aria-label', `Remove ${b.name}`);
      remove.addEventListener('click', () => removeBlock(i));
      li.append(name, where, count, remove);
      li.addEventListener('mouseenter', () => selectBlock(i));
      ul.append(li);
    });
    const cap = state.limits ? state.limits.max_amplicons : 0;
    const over = Boolean(sane && cap && total > cap);
    const lines = [];
    if (!g.blocks.length) lines.push('Click a gene above to cover it with amplicons.');
    else if (!sane) lines.push('The target amplicon length must be greater than the minimum overlap.');
    else {
      lines.push(`${fmt(g.blocks.length)} span${g.blocks.length === 1 ? '' : 's'} · `
        + `${fmt(total)} amplicon${total === 1 ? '' : 's'} of about ${fmt(ampLen)} bp.`);
    }
    if (over) lines.push(`The web app designs up to ${fmt(cap)} amplicons: lengthen the amplicons or pick fewer genes.`);
    if (g.note) lines.push(g.note);
    const bad = Boolean(g.blocks.length) && (!sane || over);
    const note = $('blocks-note');
    note.textContent = lines.join(' ');
    note.classList.toggle('error', bad);
    note.classList.toggle('help', !bad);
    $('block-controls').hidden = false;
  }

  // ---- dragging a block's handles ----
  // The listeners sit on the window, not on the handle: each move redraws the whole diagram, so
  // the element under the pointer is gone by the time the next event arrives.
  function posFromEvent(event) {
    const rect = $('genome').getBoundingClientRect();
    const scale = rect.width ? plot.W / rect.width : 1;
    return Math.round(1 + (((event.clientX - rect.left) * scale - plot.left) / plot.plotW) * plot.G);
  }

  function moveEdge(index, edge, pos) {
    const g = state.genome;
    const b = g.blocks[index];
    if (!b) return;
    const least = minSpan();
    if (edge === 'start') b.start = Math.max(1, Math.min(pos, b.end - least + 1));
    else b.end = Math.min(g.refLen, Math.max(pos, b.start + least - 1));
  }

  function onDrag(event) {
    const d = state.genome.drag;
    if (!d) return;
    moveEdge(d.index, d.edge, snapPosition(posFromEvent(event), state.genome.edges, GRID));
    redrawGenome();
    selectBlock(d.index);
  }

  function endDrag() {
    window.removeEventListener('pointermove', onDrag);
    window.removeEventListener('pointerup', endDrag);
    window.removeEventListener('pointercancel', endDrag);
    state.genome.drag = null;
    commitBlocks(state.genome.blocks);   // a drag can push two spans together
  }

  function startDrag(event, index, edge) {
    event.preventDefault();
    state.genome.drag = { index, edge };
    selectBlock(index);
    window.addEventListener('pointermove', onDrag);
    window.addEventListener('pointerup', endDrag);
    window.addEventListener('pointercancel', endDrag);
  }

  function nudge(event, index, edge) {
    const step = event.shiftKey ? 100 : GRID;
    const delta = event.key === 'ArrowLeft' ? -step : event.key === 'ArrowRight' ? step : 0;
    if (!delta) return;
    event.preventDefault();
    const b = state.genome.blocks[index];
    if (!b) return;
    moveEdge(index, edge, (edge === 'start' ? b.start : b.end) + delta);
    commitBlocks(state.genome.blocks);
    const same = document.querySelector(`#genome .handle[data-block="${index}"][data-edge="${edge}"]`);
    if (same) same.focus();
  }

  // ---- loading an annotation ----
  function useAnnotation(ann) {
    const g = state.genome;
    const what = ann.organism || ann.accession || 'This record';
    // The record is kept even when it cannot be drawn yet, so sequences dropped afterwards draw
    // it without picking the organism again. A map always belongs to one record and one
    // reference, so it goes now and is rebuilt below.
    g.ann = ann;
    g.map = null;
    g.blocks = [];
    g.note = '';
    const stop = (message) => {
      $('genome-figure').hidden = true;
      $('block-controls').hidden = true;
      $('genome-note').hidden = true;
      $('genome-verdict').textContent = `${what} is ready.`;
      genomeError(message);
    };
    const bases = referenceBases();
    if (!bases) {
      stop('Add your sequences first: the genes are drawn over your own reference.');
      return;
    }
    if (!ann.sequence) {
      stop(`${ann.accession || 'This record'} came without its sequence, so its genes cannot `
        + 'be placed on your reference. Load the GenBank file including the bases under ORIGIN.');
      return;
    }
    genomeError('');
    g.refId = $('reference_id').value;
    g.refLen = bases.length;
    g.map = buildAnchorMap(ann.sequence, bases);
    const { convert, record } = axis();
    g.edges = ann.features.flatMap((f) => [convert(f.start), convert(f.end)]);
    const notes = [];
    if (record) {
      notes.push('Your sequences and this record share too little exact sequence to place the genes '
        + 'on your reference, so the ruler below shows the record’s own positions and the spans '
        + 'follow them. Check them before you run, or pick another organism.');
    }
    if (ann.skipped.length) {
      const kinds = [...new Set(ann.skipped.map((s) => s.kind))].join(', ');
      notes.push(`${fmt(ann.skipped.length)} feature${ann.skipped.length === 1 ? '' : 's'} `
        + `could not be read and ${ann.skipped.length === 1 ? 'is' : 'are'} not drawn (${kinds}).`);
    }
    if (!ann.features.length) {
      notes.push('This record has no genes or proteins to click; cover the whole genome instead, '
        + 'or load a record that carries an annotation.');
    }
    const note = $('genome-note');
    note.textContent = notes.join(' ');
    note.hidden = !notes.length;
    note.classList.toggle('error', record);
    $('genome-verdict').textContent = `${what} · `
      + `${fmt(ann.features.length)} feature${ann.features.length === 1 ? '' : 's'} · `
      + `${fmt(g.map.count)} anchor${g.map.count === 1 ? '' : 's'} over ${pct(g.map.covered * 100)} `
      + `of the record, on your ${fmt(g.refLen)} bp reference`;
    $('genome-detail').textContent = 'Click a gene or protein to cover it with amplicons.';
    redrawGenome();
  }

  // The map holds for one reference only, so a new reference rebuilds it and drops the spans.
  // `force` is for new sequences, which can leave the reference's name unchanged.
  function refreshGenome(force) {
    const g = state.genome;
    if (!state.seqs) {
      g.map = null;
      g.blocks = [];
      $('genome-figure').hidden = true;
      $('block-controls').hidden = true;
      $('genome-note').hidden = true;
      $('genome-verdict').textContent = 'Add your sequences on the left, then pick the organism.';
      return;
    }
    if (!g.ann) {
      $('genome-verdict').textContent = 'Pick the organism your sequences come from to see its genes.';
      return;
    }
    if (!force && g.map && g.refId === $('reference_id').value) {
      redrawGenome();   // the panel may have changed width while it was hidden
      return;
    }
    const had = g.blocks.length;
    useAnnotation(g.ann);
    if (had && g.map) {
      g.note = 'Your reference changed, so the spans you had picked were cleared.';
      renderBlocks();
    }
  }

  async function loadAccession(accession) {
    const id = accession.trim();
    const g = state.genome;
    if (!id) return;
    if (!state.seqs) {
      genomeError('Add your sequences first: the genes are drawn over your own reference.');
      return;
    }
    genomeError('');
    $('genome-verdict').textContent = `Asking NCBI for ${id}…`;
    const token = {};
    g.pending = token;
    try {
      const ann = await fetchAnnotation(id);
      if (g.pending !== token) return;   // the user picked something else meanwhile
      useAnnotation(ann);
    } catch (err) {
      if (g.pending !== token) return;
      genomeError(err.message);
      $('genome-verdict').textContent = 'Pick the organism your sequences come from to see its genes.';
    }
  }

  async function loadGenbankFile(file) {
    if (!file) return;
    $('organism').value = '';
    state.genome.pending = null;
    try {
      useAnnotation(parseGenbank(await file.text()));
    } catch (err) {
      genomeError(err.message);
    }
    $('genbank').value = '';
  }

  ORGANISMS.forEach((o) => $('organism').add(new Option(`${o.label} — ${o.accession}`, o.accession)));
  $('organism').addEventListener('change', (e) => { $('accession').value = ''; loadAccession(e.target.value); });
  $('accession-load').addEventListener('click', () => { $('organism').value = ''; loadAccession($('accession').value); });
  $('accession').addEventListener('keydown', (e) => {
    if (e.key === 'Enter') { e.preventDefault(); $('accession-load').click(); }
  });
  $('genbank').addEventListener('change', (e) => loadGenbankFile(e.target.files[0]));
  ['amplicon_len', 'gene_overlap'].forEach((id) => $(id).addEventListener('input', () => {
    if (state.genome.map) renderBlocks();
  }));
  $('whole-genome').addEventListener('click', () => {
    if (state.genome.refLen) commitBlocks([{ name: 'genome', start: 1, end: state.genome.refLen }]);
  });
  $('reference_id').addEventListener('change', () => { if (mode() === 'genes') refreshGenome(); });
  window.addEventListener('resize', () => { if (mode() === 'genes') redrawGenome(); });

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
    genes: 'Click the genes you want on the diagram above. Each span you pick is tiled with overlapping amplicons named after it, so you get primers.tsv (ready to order), amplicons.tsv, scheme.bed and report.txt for those regions only.',
    check: 'Add your sequences and paste your primer pairs. Each pair is placed on the alignment and scored on every sequence, with the reason for each miss.',
  };

  function showIdle() {
    resetResults();
    const m = mode();
    const verdict = $('verdict');
    $('empty').hidden = false;
    $('empty').querySelector('p').textContent = OUTCOME[m];
    $('steps').hidden = m === 'check' || m === 'genes';   // the steps describe the whole-genome run
    const p = drawPlan();
    if (p) {
      verdict.className = 'verdict plan';
      verdict.textContent = `Planned: ${p.n} amplicon${p.n === 1 ? '' : 's'} of about ${fmt(p.L)} bp across ${fmt(p.G)} bp.`;
    } else if (m === 'tile') {
      verdict.textContent = 'Your tiling will appear here.';
    } else if (m === 'genes') {
      verdict.textContent = 'The amplicons covering the genes you pick will appear here.';
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
      const problems = r[`${what}_structure`];
      if (problems && problems !== 'ok') {
        const warn = document.createElement('span');
        warn.className = 'facts warn';
        warn.textContent = problems;
        td.append(warn);
      }
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

  const DIMER_KCAL = -9; // the same limit as vhirprimerlab.structure.DIMER_DG

  function checkNote(r) {
    const dg = r.pair_dimer_dg === '' || r.pair_dimer_dg == null ? null : +r.pair_dimer_dg;
    const dimer = dg !== null && dg < DIMER_KCAL ? `forward and reverse form a dimer (ΔG ${dg.toFixed(1).replace('-', '−')} kcal/mol)` : '';
    return [r.note, dimer].filter(Boolean).join(' · ');
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
        leadCell(nameOf(r), m === 'check' ? checkNote(r) : `rank ${r.rank}`),
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

  // ---------- primer files ----------
  // Any line of name/forward/reverse (or name/sequence) separated by tabs, commas, semicolons or
  // spaces becomes a tab-separated line; quotes from CSV files are dropped.
  function primerTable(text) {
    return `${text.split(/\r?\n/).map((l) => l.trim()).filter(Boolean)
      .map((l) => (l.startsWith('#') ? l : l.split(/[\t,;]|\s+/).map((f) => f.replace(/^"|"$/g, '').trim())
        .filter(Boolean).join('\t'))).join('\n')}\n`;
  }

  async function inflate(bytes, method) {
    if (method === 0) return bytes;
    if (method !== 8 || typeof DecompressionStream === 'undefined') throw new Error('unsupported');
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate-raw'));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }

  // The files inside a ZIP (an .xlsx is one), read from its central directory.
  async function unzip(buffer) {
    const v = new DataView(buffer);
    const bytes = new Uint8Array(buffer);
    let eocd = -1;
    for (let i = bytes.length - 22; i >= Math.max(0, bytes.length - 65557); i -= 1) {
      if (v.getUint32(i, true) === 0x06054b50) { eocd = i; break; }
    }
    if (eocd < 0) throw new Error('not a zip');
    const files = new Map();
    let p = v.getUint32(eocd + 16, true);
    const dec = new TextDecoder();
    for (let n = v.getUint16(eocd + 10, true); n > 0; n -= 1) {
      if (v.getUint32(p, true) !== 0x02014b50) throw new Error('bad zip');
      const method = v.getUint16(p + 10, true);
      const size = v.getUint32(p + 20, true);
      const nameLen = v.getUint16(p + 28, true);
      const skip = nameLen + v.getUint16(p + 30, true) + v.getUint16(p + 32, true);
      const local = v.getUint32(p + 42, true);
      const name = dec.decode(bytes.subarray(p + 46, p + 46 + nameLen));
      const start = local + 30 + v.getUint16(local + 26, true) + v.getUint16(local + 28, true);
      files.set(name, () => inflate(bytes.subarray(start, start + size), method).then((b) => dec.decode(b)));
      p += 46 + skip;
    }
    return files;
  }

  // Rows of the first worksheet of an .xlsx file, as tab-separated text.
  async function xlsxToText(buffer) {
    const files = await unzip(buffer);
    const xml = async (name) => new DOMParser().parseFromString(await files.get(name)(), 'application/xml');
    const tags = (doc, tag) => Array.from(doc.getElementsByTagNameNS('*', tag));
    const shared = files.has('xl/sharedStrings.xml')
      ? tags(await xml('xl/sharedStrings.xml'), 'si').map((si) => tags(si, 't').map((t) => t.textContent).join(''))
      : [];
    const book = await xml('xl/workbook.xml');
    const id = tags(book, 'sheet')[0].getAttributeNS('http://schemas.openxmlformats.org/officeDocument/2006/relationships', 'id');
    const rel = tags(await xml('xl/_rels/workbook.xml.rels'), 'Relationship').find((r) => r.getAttribute('Id') === id);
    const target = rel.getAttribute('Target');
    const sheet = await xml(target.startsWith('/') ? target.slice(1) : `xl/${target}`);
    return tags(sheet, 'row').map((row) => tags(row, 'c').map((c) => {
      const t = c.getAttribute('t');
      const val = tags(c, 'v')[0];
      if (t === 's' && val) return shared[+val.textContent] || '';
      if (t === 'inlineStr') return tags(c, 't').map((x) => x.textContent).join('');
      return val ? val.textContent : '';
    }).join('\t')).join('\n');
  }

  function primersError(message) {
    $('primers-error').textContent = message || '';
    $('primers-error').hidden = !message;
  }

  async function readPrimersFile(file) {
    primersError('');
    if (!file) return;
    if (/\.xls$/i.test(file.name)) {
      primersError(`${file.name} is an old Excel file. Save it as .xlsx or CSV and load it again.`);
      return;
    }
    try {
      const buffer = await file.arrayBuffer();
      const zip = new Uint8Array(buffer, 0, 2);
      const text = zip[0] === 0x50 && zip[1] === 0x4b ? await xlsxToText(buffer) : new TextDecoder().decode(buffer);
      $('primers').value = primerTable(text.replace(/^\uFEFF/, '')).trim();
    } catch (e) {
      primersError(`Could not read ${file.name}. Save it as .xlsx, CSV or TSV and load it again.`);
    }
    $('primers-file').value = '';
  }

  // ---------- running a job ----------
  function settings(m) {
    const out = {};
    if (m === 'tile') out.amplicons = Math.round(num('amplicons'));
    if (m === 'genes') {
      // The picked spans and their two controls are the mode itself, not advanced settings.
      out.tile_regions = state.genome.blocks.map((b) => `${b.name}:${b.start}-${b.end}`);
      out.amplicon_len = Math.round(num('amplicon_len'));
      out.overlap = Math.round(num('gene_overlap'));
    }
    if (!advanced()) return out; // Basic: the recommended settings
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
    if (m === 'genes') {
      if (!state.genome.blocks.length) {
        formError('Click a gene on the diagram to choose what to amplify.');
        return;
      }
      if (!(num('amplicon_len') > num('gene_overlap'))) {
        formError('The target amplicon length must be greater than the minimum overlap.');
        $('amplicon_len').focus();
        return;
      }
    }
    if (m === 'check' && !$('primers').value.trim()) {
      formError('Paste your primer pairs first.');
      $('primers').focus();
      return;
    }
    let primers = '';
    if (m === 'check') primers = primerTable($('primers').value);
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
  document.querySelectorAll('input[name="level"]').forEach((r) => r.addEventListener('change', applyLevel));
  ['max_mismatches', 'three_prime_window'].forEach((id) => $(id).addEventListener('input', drawRule));
  $('primers-file').addEventListener('change', (e) => readPrimersFile(e.target.files[0]));
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

  state.page = PAGES.includes(location.hash.slice(1)) ? location.hash.slice(1) : 'design';
  try {
    if (localStorage.getItem(LEVEL_KEY) === 'advanced') document.querySelector('input[name="level"][value="advanced"]').checked = true;
  } catch (e) { /* private window */ }
  applyLevel();
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
