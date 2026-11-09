// node --test tests/browser/blocks.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  PAD, achievedPad, ampliconCount, blockFromFeature, lossLimits, lossNotes, maxAmplicons,
  mergeBlocks, minAmpliconLen, smallestSpan, snapPosition, targetLength,
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

// ---- the padding a block actually achieved, and the loss limits that follow ----------------

test('a padded block records where its feature sits, so its padding can be read back', () => {
  const b = blockFromFeature({ name: 'F', start: 2912, end: 4718 },
                             { convert: (p) => p, limit: 13232 });
  assert.deepEqual([b.start, b.end], [2812, 4818]);
  assert.deepEqual(achievedPad(b), { start: 100, end: 100 });
});

test('a block clamped at the reference start achieved no padding there', () => {
  // hMPV's N gene maps to 1-1185 of a 13,232 bp reference: the span cannot start before base 1.
  const b = blockFromFeature({ name: 'N', start: 1, end: 1185 },
                             { convert: (p) => p, limit: 13232 });
  assert.deepEqual([b.start, b.end], [1, 1285]);
  assert.deepEqual(achievedPad(b), { start: 0, end: 100 });
});

test('a block with no feature behind it has no padding at all', () => {
  assert.deepEqual(achievedPad({ name: 'genome', start: 1, end: 13232 }), { start: 0, end: 0 });
});

test('moving a handle changes the padding, because it is measured live', () => {
  const b = blockFromFeature({ name: 'F', start: 2912, end: 4718 },
                             { convert: (p) => p, limit: 13232 });
  b.start = 2852;                       // dragged 40 bp inward, 60 of the padding left
  assert.deepEqual(achievedPad(b), { start: 60, end: 100 });
  b.start = 3012;                       // dragged 100 bp past the feature's own start
  assert.equal(achievedPad(b).start, -100);
});

test('the loss limits are the least padding any block achieved, capped and floored', () => {
  const full = { name: 'F', start: 2812, end: 4818, feat: { start: 2912, end: 4718 } };
  const atStart = { name: 'N', start: 1, end: 1285, feat: { start: 1, end: 1185 } };
  assert.deepEqual(lossLimits([full]), { start: 100, end: 100 });
  // The start is floored at the engine's own default: a limit of 0 is shorter than any primer
  // and no tiling could ever meet it.
  assert.deepEqual(lossLimits([atStart]), { start: 50, end: 100 });
  assert.deepEqual(lossLimits([full, atStart]), { start: 50, end: 100 });
  // Widening a short feature to the smallest span gives more padding than the engine is asked
  // for; asking for more than the padding would buy nothing, so it is capped.
  const wide = { name: 'S', start: 1, end: 400, feat: { start: 200, end: 202 } };
  assert.deepEqual(lossLimits([wide]), { start: 100, end: 100 });
});

test('a run with no blocks, or only the whole genome, asks for the engine default', () => {
  assert.deepEqual(lossLimits([]), { start: 50, end: 50 });
  assert.deepEqual(lossLimits([{ name: 'genome', start: 1, end: 13232 }]), { start: 50, end: 50 });
});

test('the loss notes name the blocks that can lose bases of their own feature', () => {
  const full = { name: 'F', start: 2812, end: 4818, feat: { start: 2912, end: 4718 } };
  const atStart = { name: 'N', start: 1, end: 1285, feat: { start: 1, end: 1185 } };
  const blocks = [atStart, full];
  const limits = lossLimits(blocks);
  // N has no padding at its start, so the 50 nt the engine may leave outside the inserts come
  // out of N itself. F keeps its full padding and loses nothing.
  assert.deepEqual(lossNotes(blocks, limits), [{ name: 'N', start: 50, end: 0 }]);
  assert.deepEqual(lossNotes([full], lossLimits([full])), []);
});

test('a feature filling its whole record loses bases at both ends', () => {
  // The influenza HA record: the CDS is 21-1721 of 1,733 nt, so neither end can be padded.
  const b = blockFromFeature({ name: 'HA', start: 21, end: 1721 },
                             { convert: (p) => p, limit: 1733 });
  assert.deepEqual(achievedPad(b), { start: 20, end: 12 });
  const limits = lossLimits([b]);
  assert.deepEqual(limits, { start: 50, end: 50 });
  assert.deepEqual(lossNotes([b], limits), [{ name: 'HA', start: 30, end: 38 }]);
});

test('the whole-genome block is in the limits but never in the notes', () => {
  const genome = [{ name: 'genome', start: 1, end: 13232 }];
  assert.deepEqual(lossNotes(genome, lossLimits(genome)), []);
});

test('merging keeps the outer blocks padding, which is what the merged span has', () => {
  const left = { name: 'A', start: 1, end: 1285, feat: { start: 1, end: 1185 } };
  const right = { name: 'B', start: 1286, end: 2000, feat: { start: 1386, end: 1900 } };
  const [m] = mergeBlocks([left, right]);
  assert.deepEqual([m.name, m.start, m.end], ['A_B', 1, 2000]);
  assert.deepEqual(achievedPad(m), { start: 0, end: 100 });
});

test('merging a nested block keeps the padding of the block that spans it', () => {
  const outer = { name: 'A', start: 1, end: 2000, feat: { start: 101, end: 1900 } };
  const inner = { name: 'B', start: 500, end: 700, feat: { start: 550, end: 650 } };
  const [m] = mergeBlocks([outer, inner]);
  assert.deepEqual([m.start, m.end], [1, 2000]);
  assert.deepEqual(achievedPad(m), { start: 100, end: 100 });
});

// ---- the amplicon length a span needs ------------------------------------------------------

test('the most amplicons that fit matches the engine, single amplicons included', () => {
  assert.equal(maxAmplicons(1807, 50, 25), 8);
  assert.equal(maxAmplicons(1807, 120, 25), 5);
  assert.equal(maxAmplicons(13232, 50, 25), 65);
  // A span with room for exactly one amplicon: the neighbour rule does not apply to it.
  assert.equal(maxAmplicons(64, 50, 25), 1);
  assert.equal(maxAmplicons(63, 50, 25), 0);
});

test('the target amplicon length rounds as the engine rounds', () => {
  assert.equal(targetLength(1807, 6, 50), 343);
  assert.equal(targetLength(1807, 8, 50), 270);
  assert.equal(targetLength(1807, 7, 120), 361);
  assert.equal(targetLength(1807, 2, 50), 928);   // 928.5: Python's round() is half to even
});

test('the shortest amplicon length that tiles a span is the one the engine would accept', () => {
  // hMPV's F block, 1,807 bp. At the default overlap 400 bp amplicons work and 250 bp do not;
  // the boundary is 270, and at an overlap of 120 nothing under 458 bp can tile it at all.
  const F = { name: 'F', start: 2812, end: 4618 };
  assert.equal(F.end - F.start + 1, 1807);
  assert.equal(minAmpliconLen(F, 50, 25), 270);
  assert.equal(minAmpliconLen(F, 120, 25), 458);
  assert.equal(minAmpliconLen(F, 0, 25), 130);
  // One below the boundary needs more amplicons than the geometry allows.
  assert.ok(ampliconCount(F, 270, 50) <= maxAmplicons(1807, 50, 25));
  assert.ok(ampliconCount(F, 269, 50) > maxAmplicons(1807, 50, 25));
});

test('a span too short for even one amplicon has no amplicon length that works', () => {
  assert.equal(minAmpliconLen({ start: 1, end: 40 }, 0, 25), 0);
  assert.equal(minAmpliconLen({ start: 1, end: 64 }, 50, 25), 64);
});
