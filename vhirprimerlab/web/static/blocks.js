/* Turn what the user clicks into the spans they want primers for. This is pure logic: the page
   draws and drags, and hands the results here. It touches no DOM, storage or network. */

// Padding each side of a clicked feature. The anchor map is accurate to about 13 nt over 13 kb,
// and the engine may leave up to --max-end-loss bases outside the insert at a span edge, so the
// feature itself must end up inside the amplicons rather than under a primer.
export const PAD = 50;
const END_LOSS = 50;
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
  return { name: feature.name, start, end };
}

// Blocks that touch or overlap become one, named by joining the parts with an underscore: the
// name ends up in primer names (A_B_1_LEFT) and scheme.bed, where '+' is awkward.
export function mergeBlocks(blocks) {
  const sorted = [...blocks].sort((x, y) => x.start - y.start);
  const merged = [];
  for (const next of sorted) {
    const current = merged[merged.length - 1];
    if (current && next.start <= current.end + 1) {
      current.end = Math.max(current.end, next.end);
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
