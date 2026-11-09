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

test('overlapping blocks are merged to their union', () => {
  const merged = mergeBlocks([{ name: 'A', start: 1, end: 150 },
                              { name: 'B', start: 100, end: 300 }]);
  assert.deepEqual(merged.map((b) => [b.name, b.start, b.end]), [['A_B', 1, 300]]);
});

test('a nested block does not shrink the block around it', () => {
  const merged = mergeBlocks([{ name: 'A', start: 1, end: 500 },
                              { name: 'B', start: 100, end: 200 }]);
  assert.deepEqual(merged.map((b) => [b.name, b.start, b.end]), [['A_B', 1, 500]]);
});

test('a chain of touching blocks merges into one, whatever the input order', () => {
  const merged = mergeBlocks([{ name: 'C', start: 201, end: 300 },
                              { name: 'A', start: 1, end: 100 },
                              { name: 'B', start: 101, end: 200 }]);
  assert.deepEqual(merged.map((b) => [b.name, b.start, b.end]), [['A_B_C', 1, 300]]);
});

test('separate blocks stay separate and sorted', () => {
  // Unsorted, the loop would reach A while B is current, and wrongly merge the two.
  const merged = mergeBlocks([{ name: 'B', start: 500, end: 600 },
                              { name: 'A', start: 1, end: 100 }]);
  assert.deepEqual(merged.map((b) => [b.name, b.start, b.end]), [['A', 1, 100], ['B', 500, 600]]);
});

test('an empty list and a single block come back unchanged', () => {
  assert.deepEqual(mergeBlocks([]), []);
  assert.deepEqual(mergeBlocks([{ name: 'A', start: 1, end: 5 }]), [{ name: 'A', start: 1, end: 5 }]);
});

test('the amplicon count matches the engine', () => {
  assert.equal(ampliconCount({ start: 1, end: 400 }, 400, 50), 1);
  assert.equal(ampliconCount({ start: 1, end: 750 }, 400, 50), 2);
  assert.equal(ampliconCount({ start: 1, end: 751 }, 400, 50), 3);
});

test('an amplicon no longer than its overlap is refused, as the engine refuses it', () => {
  assert.throws(() => ampliconCount({ start: 1, end: 400 }, 50, 50), /larger than overlap/);
  assert.throws(() => ampliconCount({ start: 1, end: 400 }, 40, 50), /larger than overlap/);
});

test('dragging snaps to feature edges first, then to the grid', () => {
  assert.equal(snapPosition(3055, [3052, 4671]), 3052);
  assert.equal(snapPosition(3104, [3052, 4671]), 3100);
});

test('a drag near two edges snaps to the nearer one', () => {
  assert.equal(snapPosition(3057, [3052, 3058]), 3058);
  assert.equal(snapPosition(3056, [3058, 3052]), 3058);
});

test('a block shorter than minSpan grows to it, taking from the side with room', () => {
  // The feature sits near the reference start, so the left side can give only 2 bases.
  const b = blockFromFeature({ name: 'S', start: 3, end: 5 },
                             { convert: (p) => p, pad: 0, endLoss: 0, minSpan: 64, limit: 200 });
  assert.deepEqual([b.start, b.end], [1, 64]);
});

test('a block stops at the reference bounds when minSpan cannot be reached', () => {
  const b = blockFromFeature({ name: 'S', start: 3, end: 5 },
                             { convert: (p) => p, pad: 0, endLoss: 0, minSpan: 64, limit: 40 });
  assert.deepEqual([b.start, b.end], [1, 40]);
});

test('the smallest span is two primers plus one base at the default tolerance', () => {
  assert.equal(smallestSpan(25, 0.20), 64);
  assert.equal(PAD, 50);
});
