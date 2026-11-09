// node --test tests/browser/settings.test.mjs
// settings() builds the payload a run is started with. It lives inside app.js's IIFE, where the
// page's DOM is its closure, so rather than fake a browser the test lifts the function's own
// source out of app.js and runs it against stubs. What is pinned is then the real code, including
// which branch returns before which: the Basic early return is the whole point of these tests.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { lossLimits, minAmpliconLen } from '../../vhirprimerlab/web/static/blocks.js';

const source = readFileSync(new URL('../../vhirprimerlab/web/static/app.js', import.meta.url), 'utf8');

// The source of one function of the IIFE, by counting braces from its opening one. Every brace in
// settings() is balanced (its template literals hold no stray one), so counting is enough.
function lift(name) {
  const head = source.indexOf(`  function ${name}(`);
  assert.notEqual(head, -1, `${name}() is not in app.js`);
  let depth = 0;
  for (let i = source.indexOf('{', head); i < source.length; i += 1) {
    if (source[i] === '{') depth += 1;
    else if (source[i] === '}') {
      depth -= 1;
      if (depth === 0) return source.slice(head, i + 1);
    }
  }
  throw new Error(`${name}() does not close`);
}

const FIELDS = {
  amplicons: '10', amplicon_len: '400', gene_overlap: '50', primer_len_max: '25',
  regions: '', reference_id: '',
};

// The page as settings() sees it: the fields it reads, the level, the picked spans, and the
// [data-setting] inputs of the Advanced panel.
function build({ level = 'basic', blocks = [], inputs = [], ...fields } = {}) {
  const values = { ...FIELDS, ...fields };
  const el = (id) => ({ value: values[id] ?? '', checked: Boolean(values[`${id}:checked`]) });
  return new Function('$', 'num', 'advanced', 'state', 'document', 'lossLimits',
    `${lift('settings')}; return settings;`)(
    el,
    (id) => Number(values[id] ?? ''),
    () => level === 'advanced',
    { genome: { blocks }, seqs: { ids: ['first.1', 'second.1'] } },
    { querySelectorAll: () => inputs },
    lossLimits,
  );
}

const BLOCKS = [{ name: 'F', start: 2912, end: 4718 }];

test('the genes payload names the reference even in Basic', () => {
  // The engine's own default reference is the first sequence left after length filtering, which
  // need not be the one the diagram was drawn over. Sending it is what keeps the two in step.
  const out = build({ level: 'basic', blocks: BLOCKS })('genes');
  assert.equal(out.reference_id, 'first.1');
  assert.deepEqual(out.tile_regions, ['F:2912-4718']);
});

test('the genes payload keeps the reference the user chose, in Basic too', () => {
  // Advanced lets them pick a non-first reference; switching back to Basic keeps the select's
  // value and the diagram, so the payload has to keep it as well.
  const out = build({ level: 'basic', blocks: BLOCKS, reference_id: 'second.1' })('genes');
  assert.equal(out.reference_id, 'second.1');
});

test('the genes payload names the reference in Advanced as well', () => {
  const out = build({ level: 'advanced', blocks: BLOCKS, reference_id: 'second.1' })('genes');
  assert.equal(out.reference_id, 'second.1');
});

test('an Advanced reference input cannot overwrite the reference already in the payload', () => {
  const input = {
    tagName: 'SELECT', value: 'second.1', selectedIndex: 1, defaultValue: '',
    dataset: { setting: 'reference_id' },
    closest: () => ({ dataset: { modes: 'tile regions check genes' } }),
  };
  const out = build({ level: 'advanced', blocks: BLOCKS, reference_id: 'second.1', inputs: [input] })('genes');
  assert.equal(out.reference_id, 'second.1');
});

test('the other modes do not gain a reference of their own', () => {
  assert.equal('reference_id' in build({ level: 'basic' })('tile'), false);
  assert.equal('reference_id' in build({ level: 'basic' })('check'), false);
});

test('the loss limits sent are the padding the spans achieved, not what they asked for', () => {
  // F keeps its full 100 nt on both sides; N is clamped at base 1 and has none at its start.
  const full = { name: 'F', start: 2812, end: 4818, feat: { start: 2912, end: 4718 } };
  const atStart = { name: 'N', start: 1, end: 1285, feat: { start: 1, end: 1185 } };
  const out = build({ level: 'basic', blocks: [atStart, full] })('genes');
  assert.equal(out.max_start_loss, 50);
  assert.equal(out.max_end_loss, 100);
  assert.deepEqual(out.tile_regions, ['N:1-1285', 'F:2812-4818']);
  const padded = build({ level: 'basic', blocks: [full] })('genes');
  assert.equal(padded.max_start_loss, 100);
  assert.equal(padded.max_end_loss, 100);
});

// ---- what the page says before a run, about the two block controls -------------------------

// blockProblem() reads the same two controls as the payload above, and the same primer length, so
// it is lifted the same way. fmt() is the page's thousands separator.
function buildProblem({ blocks = [], primer_len_max: lenMax = '25' } = {}) {
  return new Function('fmt', 'num', 'state', 'minAmpliconLen',
    `${lift('blockProblem')}; return blockProblem;`)(
    (n) => Number(n).toLocaleString('en-US'),
    (id) => Number(id === 'primer_len_max' ? lenMax : ''),
    { genome: { blocks } },
    minAmpliconLen,
  );
}

const F = [{ name: 'F', start: 2812, end: 4618 }];   // hMPV's F block: 1,807 bp

test('the default controls tile the F block, and say nothing', () => {
  assert.equal(buildProblem({ blocks: F })(400, 50), '');
});

test('an overlap too large for the amplicons names the amplicon length needed, not an N', () => {
  // The engine would say: Region 'F': N too large: amplicons would be ~361 bp; use at most 5.
  const msg = buildProblem({ blocks: F })(400, 120);
  assert.match(msg, /^400 bp amplicons cannot tile F \(1,807 bp\) at a 120 bp overlap/);
  assert.match(msg, /shortest that works here is 458 bp/);
  assert.equal(/\bN\b/.test(msg), false);
});

test('amplicons too short for any tiling name the shortest that works', () => {
  assert.match(buildProblem({ blocks: F })(200, 50), /shortest that works here is 270 bp/);
  assert.match(buildProblem({ blocks: F })(250, 50), /shortest that works here is 270 bp/);
  assert.equal(buildProblem({ blocks: F })(270, 50), '');
});

test('the binding span is the one that needs the longest amplicons', () => {
  const blocks = [{ name: 'N', start: 1, end: 1285 }, ...F];
  const msg = buildProblem({ blocks })(200, 50);
  assert.match(msg, /cannot tile F/);     // 1,807 bp needs 270; N's 1,285 bp needs only 192
  assert.equal(buildProblem({ blocks })(270, 50), '');
});

test('longer primers raise the amplicon length a span needs', () => {
  assert.match(buildProblem({ blocks: F, primer_len_max: '30' })(270, 50),
    /shortest that works here is 301 bp/);
});

test('an amplicon no longer than the overlap is still refused first', () => {
  assert.match(buildProblem({ blocks: F })(50, 50), /must be greater than the minimum overlap/);
  assert.match(buildProblem({ blocks: F })(0, 50), /must be greater than the minimum overlap/);
});

test('a span too short for one amplicon says so, instead of an amplicon length', () => {
  assert.match(buildProblem({ blocks: [{ name: 'tiny', start: 1, end: 40 }] })(400, 0),
    /^tiny is only 40 bp: too short for one amplicon with primers of up to 25 bases\.$/);
});
