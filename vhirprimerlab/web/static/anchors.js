/* Map a record's coordinates onto a reference with the exact k-mers the two sequences share.
   No aligner and no network: the browser needs this to put NCBI gene positions onto the user's
   own reference. */

const MIN_ANCHORS = 10;     // fewer anchors than this say little about the whole record
const MIN_COVERED = 0.5;    // and so does a map that spans less than half of it

export function buildAnchorMap(record, reference, { k = 20, step = 10 } = {}) {
  const rec = record.toUpperCase();
  const ref = reference.toUpperCase();

  // Only k-mers that occur once in the reference can anchor: a repeat would map anywhere.
  const where = new Map();
  for (let i = 0; i + k <= ref.length; i++) {
    const kmer = ref.slice(i, i + k);
    where.set(kmer, where.has(kmer) ? -1 : i);
  }
  const hits = [];
  for (let j = 0; j + k <= rec.length; j += step) {
    const i = where.get(rec.slice(j, j + k));
    if (i !== undefined && i >= 0) hits.push([j + 1, i + 1]);
  }

  // Longest strictly increasing run of reference positions (patience sorting). The record
  // positions already increase, so this drops hits that contradict the main diagonal.
  const tails = [];     // tails[p]: smallest reference position ending a chain of length p + 1
  const tailHit = [];   // the hit index behind each tail
  const back = new Array(hits.length).fill(-1);
  hits.forEach(([, r], h) => {
    let lo = 0;
    let hi = tails.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (tails[mid] < r) lo = mid + 1; else hi = mid;
    }
    tails[lo] = r;
    tailHit[lo] = h;
    back[h] = lo ? tailHit[lo - 1] : -1;
  });
  const anchors = [];
  for (let h = tails.length ? tailHit[tails.length - 1] : -1; h >= 0; h = back[h]) {
    anchors.push(hits[h]);
  }
  anchors.reverse();

  const count = anchors.length;
  const mappedRange = count
    ? [anchors[0][0], anchors[count - 1][0] + k - 1]
    : [0, 0];
  const covered = count && rec.length ? (mappedRange[1] - mappedRange[0] + 1) / rec.length : 0;

  function convert(pos) {
    if (!count) return pos;
    if (pos <= anchors[0][0]) return anchors[0][1] + (pos - anchors[0][0]);
    const last = anchors[count - 1];
    if (pos >= last[0]) return last[1] + (pos - last[0]);
    let lo = 0;           // the first anchor at or after pos
    let hi = count - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (anchors[mid][0] < pos) lo = mid + 1; else hi = mid;
    }
    const [r1, f1] = anchors[lo];
    if (r1 === pos) return f1;
    const [r0, f0] = anchors[lo - 1];
    return Math.round(f0 + ((pos - r0) / (r1 - r0)) * (f1 - f0));
  }

  return {
    anchors, count, covered,
    trusted: count >= MIN_ANCHORS && covered >= MIN_COVERED,
    mappedRange, convert,
  };
}
