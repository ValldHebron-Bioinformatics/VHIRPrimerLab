// node --test tests/browser/settings.test.mjs
// settings() builds the payload a run is started with. It lives inside app.js's IIFE, where the
// page's DOM is its closure, so rather than fake a browser the test lifts the function's own
// source out of app.js and runs it against stubs. What is pinned is then the real code, including
// which branch returns before which: the Basic early return is the whole point of these tests.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';

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
  return new Function('$', 'num', 'advanced', 'state', 'document', 'PADDING',
    `${lift('settings')}; return settings;`)(
    el,
    (id) => Number(values[id] ?? ''),
    () => level === 'advanced',
    { genome: { blocks }, seqs: { ids: ['first.1', 'second.1'] } },
    { querySelectorAll: () => inputs },
    100,
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
