// node --test tests/browser/anchors.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { buildAnchorMap } from '../../vhirprimerlab/web/static/anchors.js';

// Deterministic pseudo-random ACGT from a seeded LCG. The high bits are used because the low
// bits of a power-of-two LCG cycle with a tiny period. At 2-3 kb a repeated 20-mer is
// astronomically unlikely, and the test below checks it.
const seq = (n, seed = 1) => {
  let state = (seed * 2654435761) >>> 0;
  let out = '';
  for (let i = 0; i < n; i++) {
    state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
    out += 'ACGT'[state >>> 30];
  }
  return out;
};

test('the generator gives non-degenerate sequences without repeated 20-mers', () => {
  const s = seq(3000);
  const kmers = new Set();
  for (let i = 0; i + 20 <= s.length; i++) kmers.add(s.slice(i, i + 20));
  assert.equal(kmers.size, s.length - 19);
  for (const base of 'ACGT') {
    const share = s.split(base).length - 1;
    assert.ok(share > 600 && share < 900, `${base} share ${share}`);
  }
  assert.notEqual(seq(100, 1), seq(100, 2));
});

test('an identical record maps one to one', () => {
  const s = seq(2000);
  const map = buildAnchorMap(s, s);
  assert.ok(map.trusted);
  assert.equal(map.count, map.anchors.length);
  assert.ok(map.covered > 0.9);
  assert.equal(map.convert(1), 1);
  assert.equal(map.convert(1500), 1500);
  assert.equal(map.convert(2000), 2000);
});

test('anchors are 1-based and strictly increasing in both coordinates', () => {
  const s = seq(2000);
  const map = buildAnchorMap(s, s.slice(0, 1000) + seq(30, 9) + s.slice(1000));
  assert.deepEqual(map.anchors[0], [1, 1]);
  for (let i = 1; i < map.anchors.length; i++) {
    assert.ok(map.anchors[i][0] > map.anchors[i - 1][0]);
    assert.ok(map.anchors[i][1] > map.anchors[i - 1][1]);
  }
});

test('an insertion in the reference shifts later positions', () => {
  const s = seq(2000);
  const reference = s.slice(0, 1000) + seq(30, 9) + s.slice(1000);
  const map = buildAnchorMap(s, reference);
  assert.ok(map.trusted);
  assert.equal(map.convert(500), 500);
  assert.ok(Math.abs(map.convert(1500) - 1530) <= 15);
});

test('a deletion in the reference shifts later positions back', () => {
  const s = seq(2000);
  const map = buildAnchorMap(s, s.slice(0, 1000) + s.slice(1040));
  assert.equal(map.convert(500), 500);
  assert.ok(Math.abs(map.convert(1500) - 1460) <= 15);
});

test('positions beyond the first and last anchor extrapolate by a constant offset', () => {
  const s = seq(2000);
  const reference = seq(50, 5) + s;                  // the reference has 50 extra bases in front
  const map = buildAnchorMap(s, reference);
  assert.equal(map.convert(1), 51);
  assert.equal(map.convert(2000), 2050);
});

test('a record longer than the reference reports the mapped range only', () => {
  const s = seq(3000);
  const map = buildAnchorMap(s, s.slice(800, 1800));   // the reference is an inner gene
  assert.ok(map.mappedRange[0] >= 790 && map.mappedRange[1] <= 1810);
  assert.ok(map.covered < 0.5);
  assert.equal(map.trusted, false);
});

test('an unrelated record is not trusted', () => {
  const map = buildAnchorMap(seq(2000, 1), seq(2000, 7));
  assert.equal(map.trusted, false);
  assert.ok(map.count < 10);
});

test('an empty or too short input gives an empty, untrusted map', () => {
  const map = buildAnchorMap('ACGT', seq(500));
  assert.equal(map.count, 0);
  assert.equal(map.trusted, false);
  assert.equal(map.convert(7), 7);
});
