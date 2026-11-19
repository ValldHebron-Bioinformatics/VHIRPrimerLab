// node --test tests/browser/example.test.mjs
// Pick by gene's "Find the best match" maps each catalogue record onto the first sequence of the
// user's file with exact shared 20-mers (anchors.js). The built-in example used to be synthetic
// genomes, which share nothing with any real virus, so the button always answered "no matches".
// These tests pin both halves offline: the real example finds NC_039199.1, the synthetic one does
// not (so nobody swaps them back), without calling NCBI.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { buildAnchorMap } from '../../vhirprimerlab/web/static/anchors.js';
import { findBestMatch } from '../../vhirprimerlab/web/static/bestmatch.js';

const read = (url) => readFileSync(new URL(url, import.meta.url), 'utf8');

function fasta(text) {
  return text.split('>').slice(1).map((block) => {
    const [header, ...lines] = block.split('\n');
    return { id: header.trim().split(/\s/)[0], seq: lines.join('').replace(/\s/g, '').toUpperCase() };
  });
}

const GENES_EXAMPLE = fasta(read('../../vhirprimerlab/web/static/example-genes.fasta'));
const SYNTHETIC_EXAMPLE = fasta(read('../../vhirprimerlab/web/static/example.fasta'));
// RefSeq NC_039199.1, the human metapneumovirus entry of organisms.json (13,350 nt).
const HMPV_REFSEQ = fasta(read('./fixtures/NC_039199.fa'))[0].seq;

// 13 kb that share no 20-mer with anything real: stands in for the other catalogue viruses.
function unrelated(length, seed) {
  let x = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0;
    out += 'ACGT'[x & 3];
  }
  return out;
}

const catalogue = async (accession) => ({
  accession,
  sequence: accession === 'NC_039199.1' ? HMPV_REFSEQ : unrelated(13000, accession.length * 7919 + 17),
});
const best = (reference) => findBestMatch(
  ['NC_039199.1', 'NC_038235.1', 'NC_001781.1'], reference,
  { fetchAnnotation: catalogue, buildAnchorMap, pause: 0 });

test('the Pick by gene example is a dozen full-length human metapneumovirus genomes', () => {
  assert.equal(GENES_EXAMPLE.length, 12);
  for (const { id, seq } of GENES_EXAMPLE) {
    assert.ok(seq.length > 13000 && seq.length < 13500, `${id}: ${seq.length} bp`);
    assert.match(seq, /^[ACGTRYKMSWBDHVN]+$/, id);
  }
});

test('every genome of the example maps onto the hMPV RefSeq with a trusted anchor map', () => {
  for (const { id, seq } of GENES_EXAMPLE) {
    const map = buildAnchorMap(seq, HMPV_REFSEQ);
    assert.ok(map.trusted, `${id}: ${map.count} anchors, ${map.covered.toFixed(2)} covered`);
  }
});

test('"find the best match" on the example picks the hMPV RefSeq out of the catalogue', async () => {
  const result = await best(GENES_EXAMPLE[0].seq);
  assert.equal(result.best?.accession, 'NC_039199.1');
  assert.deepEqual(result.ranked.map((c) => c.accession), ['NC_039199.1']);
});

test('the synthetic tiling example matches nothing, which is why it is not the genes example', async () => {
  const result = await best(SYNTHETIC_EXAMPLE[0].seq);
  assert.equal(result.best, null);
});
