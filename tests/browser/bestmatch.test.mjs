// node --test tests/browser/bestmatch.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { findBestMatch } from '../../vhirprimerlab/web/static/bestmatch.js';

// Stand-ins for the real fetch and anchor map: the ranking is what is under test.
const maps = {
  A: { trusted: true, covered: 0.6, count: 40 },
  B: { trusted: true, covered: 0.97, count: 120 },
  C: { trusted: false, covered: 0.99, count: 5 },
  D: { trusted: true, covered: 0.97, count: 300 },
};
const run = (accs, extra = {}) => findBestMatch(accs, 'ACGTACGT', {
  fetchAnnotation: async (acc) => {
    if (extra.fail?.includes(acc)) throw new Error('NCBI down');
    return { accession: acc, sequence: extra.noSeq?.includes(acc) ? '' : acc };
  },
  buildAnchorMap: (seq) => maps[seq],
  pause: 0,
});

test('the trusted map covering most of its record wins; anchors break a tie', async () => {
  const r = await run(['A', 'B', 'D']);
  assert.equal(r.best.accession, 'D');
  assert.deepEqual(r.ranked.map((c) => c.accession), ['D', 'B', 'A']);
});

test('an untrusted map never wins, however much it covers', async () => {
  const r = await run(['C']);
  assert.equal(r.best, null);
  assert.deepEqual(r.ranked, []);
});

test('fetch failures and records without bases are reported, not scored', async () => {
  const r = await run(['A', 'B', 'D'], { fail: ['B'], noSeq: ['D'] });
  assert.equal(r.best.accession, 'A');
  assert.deepEqual(r.failed.sort(), ['B']);
  assert.deepEqual(r.noSequence, ['D']);
});

test('the same accession given twice, or with another case, is fetched once', async () => {
  const seen = [];
  await findBestMatch(['A', 'a', 'A'], 'ACGT', {
    fetchAnnotation: async (acc) => { seen.push(acc); return { accession: acc, sequence: 'A' }; },
    buildAnchorMap: () => maps.A,
  });
  assert.equal(seen.length, 1);
});

test('the records are spaced out so NCBI is not asked for them all at once', async () => {
  const waits = [];
  await findBestMatch(['A', 'B', 'D'], 'ACGT', {
    fetchAnnotation: async (acc) => ({ accession: acc, sequence: acc }),
    buildAnchorMap: (seq) => maps[seq],
    pause: 400, wait: async (ms) => { waits.push(ms); },
  });
  assert.deepEqual(waits, [400, 400]);
});

test('nothing to score gives no best match and no error', async () => {
  const r = await run([]);
  assert.deepEqual(r, { best: null, ranked: [], failed: [], noSequence: [] });
});
