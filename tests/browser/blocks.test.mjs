// node --test tests/browser/blocks.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  PAD, blockFromFeature, mergeBlocks, ampliconCount, snapPosition, smallestSpan,
} from '../../vhirprimerlab/web/static/blocks.js';

test('a block covers its feature with padding on both sides', () => {
  const b = blockFromFeature({ name: 'F', start: 3052, end: 4671 },
                             { convert: (p) => p - 40, limit: 13232 });
  assert.equal(b.name, 'F');
  assert.ok(b.start <= 3012 - 99 && b.end >= 4631 + 99);
});

test('padding is applied and clamped to the reference', () => {
  const b = blockFromFeature({ name: 'N', start: 1, end: 200 },
                             { convert: (p) => p, limit: 500 });
  assert.deepEqual([b.start, b.end], [1, 300]);   // clamped at the start, padded at the end
});

test('padding never leaves the reference', () => {
  const b = blockFromFeature({ name: 'N', start: 1, end: 600 },
                             { convert: (p) => p, limit: 500 });
  assert.deepEqual([b.start, b.end], [1, 500]);   // clamped at both ends
});

test('blocks that touch exactly are merged', () => {
  const merged = mergeBlocks([{ name: 'A', start: 1, end: 100 },
                              { name: 'B', start: 101, end: 200 }]);
  assert.equal(merged.length, 1);
  assert.deepEqual([merged[0].start, merged[0].end], [1, 200]);
  assert.equal(merged[0].name, 'A_B');
});

test('separate blocks stay separate and sorted', () => {
  const merged = mergeBlocks([{ name: 'B', start: 500, end: 600 },
                              { name: 'A', start: 1, end: 100 }]);
  assert.deepEqual(merged.map((b) => b.name), ['A', 'B']);
});

test('the amplicon count matches the engine', () => {
  assert.equal(ampliconCount({ start: 1, end: 400 }, 400, 50), 1);
  assert.equal(ampliconCount({ start: 1, end: 750 }, 400, 50), 2);
  assert.equal(ampliconCount({ start: 1, end: 751 }, 400, 50), 3);
});

test('dragging snaps to feature edges first, then to the grid', () => {
  assert.equal(snapPosition(3055, [3052, 4671]), 3052);
  assert.equal(snapPosition(3104, [3052, 4671]), 3100);
});

test('the smallest span is two primers plus one base at the default tolerance', () => {
  assert.equal(smallestSpan(400, 50, 25, 0.20), 64);
  assert.equal(PAD, 50);
});
