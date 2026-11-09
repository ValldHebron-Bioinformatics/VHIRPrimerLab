/* Turn what the user clicks into the spans they want primers for. This is pure logic: the page
   draws and drags, and hands the results here. It touches no DOM, storage or network. */

// Padding each side of a clicked feature. The anchor map is accurate to about 13 nt over 13 kb,
// and the engine may leave up to --max-end-loss bases outside the insert at a span edge, so the
// feature itself must end up inside the amplicons rather than under a primer.
export const PAD = 50;
// The engine's own --max-start-loss / --max-end-loss default: the span bases it may leave outside
// the inserts at an edge, the primer there included. Each block carries this much padding on top
// of PAD so that the loss lands on padding rather than on the feature, and it is also the floor of
// what a run asks for, because a limit shorter than a primer can never be met. One definition: the
// page imports it from here.
export const END_LOSS = 50;
const SNAP_EDGE = 10;   // bp: a drag this close to a feature edge snaps to the edge

// The span a block must have so one amplicon holds two primers and a base of insert.
// Mirrors the engine's single-amplicon rule after its length tolerance.
export function smallestSpan(primerLenMax = 25, tolerance = 0.20) {
  return Math.ceil((2 * primerLenMax + 1) / (1 - tolerance));
}

// `convert` maps a position on the annotation's record to the user's reference. Padding is
// clamped to the reference. A block shorter than `minSpan` grows about its centre, and a side
// that hits the reference's edge hands its share of the growth to the other side.
export function blockFromFeature(feature, { convert, pad = PAD, endLoss = END_LOSS, minSpan = 0, limit }) {
  const a = convert(feature.start);
  const b = convert(feature.end);
  const margin = pad + endLoss;
  let start = Math.max(1, Math.min(a, b) - margin);
  let end = Math.min(limit, Math.max(a, b) + margin);
  const need = minSpan - (end - start + 1);
  if (need > 0) {
    const roomLeft = start - 1;
    const roomRight = limit - end;
    let growRight = Math.min(roomRight, Math.floor(need / 2));
    let growLeft = Math.min(roomLeft, need - growRight);
    growRight = Math.min(roomRight, need - growLeft);
    start -= growLeft;
    end += growRight;
  }
  // Where the feature itself landed, kept so the padding the block really achieved can be read
  // back later: clamping at a reference edge takes padding away, and so does dragging a handle.
  return { name: feature.name, start, end, feat: { start: Math.min(a, b), end: Math.max(a, b) } };
}

// The padding a block achieved at each end: the bases between the span's edge and the feature the
// user clicked. Measured from the block as it stands, so a dragged handle is accounted for, and
// negative where a drag went inside the feature. A block with no feature behind it (the
// whole-genome shortcut) has no padding: its span edge is the reference's own.
const featOf = (b) => b.feat || { start: b.start, end: b.end };

export function achievedPad(block) {
  const feat = featOf(block);
  return { start: feat.start - block.start, end: block.end - feat.end };
}

// The --max-start-loss / --max-end-loss a run of these blocks must ask for. The engine measures
// both from the span's edge, so a limit larger than the padding there lets it leave the feature's
// own bases outside the inserts: the limit is the least padding any block achieved. Asking for
// more than the full padding buys nothing, hence the cap; a limit under END_LOSS cannot be met by
// any primer, hence the floor, which is the engine's own default.
export function lossLimits(blocks, { pad = PAD, endLoss = END_LOSS } = {}) {
  const full = pad + endLoss;
  const least = (side) => Math.min(...blocks.map((b) => achievedPad(b)[side]));
  const limit = (side) => Math.min(full, Math.max(endLoss, least(side)));
  return blocks.length ? { start: limit('start'), end: limit('end') } : { start: endLoss, end: endLoss };
}

// Bases of the feature behind each block that the engine may leave outside the amplicons under
// these limits: what is lost is the limit beyond the padding. Since the limit is the least padding
// achieved (floored at END_LOSS), only a block with less than END_LOSS of padding at an end can
// lose anything there, and a fully padded block never appears here. The results panel never
// mentions this, so the block list must.
export function lossNotes(blocks, limits) {
  const notes = [];
  for (const b of blocks) {
    if (!b.feat) continue;   // the whole-genome shortcut covers no feature in particular
    const got = achievedPad(b);
    const start = Math.max(0, limits.start - got.start);
    const end = Math.max(0, limits.end - got.end);
    if (start || end) notes.push({ name: b.name, start, end });
  }
  return notes;
}

// Blocks that touch or overlap become one, named by joining the parts with an underscore: the
// name ends up in primer names (A_B_1_LEFT) and scheme.bed, where '+' is awkward.
export function mergeBlocks(blocks) {
  const sorted = [...blocks].sort((x, y) => x.start - y.start);
  const merged = [];
  for (const next of sorted) {
    const current = merged[merged.length - 1];
    if (current && next.start <= current.end + 1) {
      // The merged span has one pair of edges, so its padding is the outermost blocks': the first
      // block's at the start (they are sorted) and, at the end, that of whichever reaches furthest.
      if (next.end > current.end) {
        current.feat = { start: featOf(current).start, end: featOf(next).end };
        current.end = next.end;
      }
      current.name = `${current.name}_${next.name}`;
    } else {
      merged.push({ ...next });
    }
  }
  return merged;
}

// Same formula as amplicons_for() in vhirprimerlab/tiling.py, which is the other copy: keep
// them in step. The contract is the ampliconCount test in tests/browser/blocks.test.mjs.
export function ampliconCount(block, ampliconLen, overlap) {
  if (ampliconLen <= overlap) throw new Error('ampliconLen must be larger than overlap');
  const span = block.end - block.start + 1;
  return Math.max(1, Math.ceil((span - overlap) / (ampliconLen - overlap)));
}

// A drag snaps to the nearest feature edge within reach, and otherwise to the grid. It does not
// clamp, so the caller must keep the result inside the reference.
export function snapPosition(pos, edges, grid = 10) {
  let best;
  for (const e of edges) {
    const d = Math.abs(pos - e);
    if (d <= SNAP_EDGE && (best === undefined || d < Math.abs(pos - best))) best = e;
  }
  return best !== undefined ? best : Math.round(pos / grid) * grid;
}
