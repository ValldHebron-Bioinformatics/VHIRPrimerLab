// node --test tests/browser/annotation.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import {
  parseLocation, parseInsdseq, parseGenbank, fetchAnnotation, fetchRecordLength, MAX_RECORD_LENGTH,
} from '../../vhirprimerlab/web/static/annotation.js';

const fixture = (name) => readFileSync(new URL(`./fixtures/${name}`, import.meta.url), 'utf8');
const hmpvXml = fixture('NC_039199.xml');
const giardiaXml = fixture('parasite-gene.xml');
const bigRecordXml = hmpvXml.replace('<INSDSeq_length>13350<', '<INSDSeq_length>600000<');
// A fetch that answers esummary with `slen` (or a 404 when undefined) and efetch with `efetch`.
const route = (efetch, slen) => async (url) => {
  if (url.includes('esummary')) {
    if (slen === undefined) return { ok: false, status: 404 };
    return { ok: true, status: 200, json: async () => ({ result: { uids: ['1'], 1: { slen } } }) };
  }
  return efetch(url);
};
const okFetch = (body, slen) => route(async () => ({ ok: true, status: 200, text: async () => body }), slen);
const noWait = async () => {};

test('a plain location', () => {
  assert.deepEqual(parseLocation('40..1224'), { strand: '+', intervals: [[40, 1224]] });
});

test('a single-base location', () => {
  assert.deepEqual(parseLocation('77'), { strand: '+', intervals: [[77, 77]] });
});

test('partial markers are ignored', () => {
  assert.deepEqual(parseLocation('<1..>300'), { strand: '+', intervals: [[1, 300]] });
});

test('a complement location keeps increasing intervals', () => {
  assert.deepEqual(parseLocation('complement(40..1224)'),
    { strand: '-', intervals: [[40, 1224]] });
});

test('a joined location', () => {
  assert.deepEqual(parseLocation('join(1..10,20..30)'),
    { strand: '+', intervals: [[1, 10], [20, 30]] });
});

test('a joined complement location', () => {
  const loc = parseLocation('join(complement(20..30),complement(1..10))');
  assert.deepEqual(loc, { strand: '-', intervals: [[1, 10], [20, 30]] });
});

test('a complement of a join', () => {
  assert.deepEqual(parseLocation('complement(join(1..10,20..30))'),
    { strand: '-', intervals: [[1, 10], [20, 30]] });
});

test('a malformed location throws', () => {
  assert.throws(() => parseLocation('nonsense'), /location/);
  assert.throws(() => parseLocation('join(1..10'), /location/);
  assert.throws(() => parseLocation('30..10'), /location/);
  assert.throws(() => parseLocation(''), /location/);
});

test('the hMPV record yields its genes and proteins', () => {
  const a = parseInsdseq(hmpvXml);
  assert.equal(a.accession, 'NC_039199.1');
  assert.equal(a.organism, 'human metapneumovirus');
  assert.equal(a.length, 13350);
  const f = a.features.find((x) => x.kind === 'CDS' && x.name === 'F');
  assert.deepEqual([f.start, f.end], [3052, 4671]);
  assert.equal(f.product, 'fusion protein');
  assert.equal(f.strand, '+');
  assert.deepEqual(f.intervals, [[3052, 4671]]);
  assert.ok(!a.features.some((x) => x.kind === 'source'));
  assert.equal(a.features.filter((x) => x.kind === 'gene').length, 8);
  assert.equal(a.features.filter((x) => x.kind === 'CDS').length, 9);
});

test('the influenza record keeps its signal and mature peptides', () => {
  const a = parseInsdseq(fixture('influenza-segment.xml'));
  assert.equal(a.length, 1701);
  const peptides = a.features.filter((x) => x.kind === 'mat_peptide');
  assert.deepEqual(peptides.map((p) => [p.product, p.start, p.end]),
    [['HA1', 52, 1032], ['HA2', 1033, 1698]]);
  const sig = a.features.find((x) => x.kind === 'sig_peptide');
  assert.deepEqual([sig.start, sig.end], [1, 51]);
});

test('an rRNA feature with no /gene is kept and named from its product', () => {
  const a = parseInsdseq(giardiaXml);
  assert.equal(a.organism, 'Giardia duodenalis');
  assert.equal(a.features.length, 1);
  const [r] = a.features;
  assert.equal(r.kind, 'rRNA');
  assert.deepEqual([r.start, r.end], [1, 1453]);
  assert.equal(r.name, '16S-like ribosomal RNA');
});

test('a feature without /gene is named from its locus_tag before its product', () => {
  const xml = giardiaXml.replace(
    '<INSDQualifier_name>product</INSDQualifier_name>',
    '<INSDQualifier_name>locus_tag</INSDQualifier_name><INSDQualifier_value>GL50803_r1</INSDQualifier_value></INSDQualifier><INSDQualifier><INSDQualifier_name>product</INSDQualifier_name>');
  const [r] = parseInsdseq(xml).features;
  assert.equal(r.name, 'GL50803_r1');
  assert.equal(r.product, '16S-like ribosomal RNA');
});

test('a GenBank file parses to the same shape', () => {
  const a = parseGenbank(fixture('hmpv.gb'));
  assert.equal(a.accession, 'NC_039199.1');
  assert.equal(a.organism, 'human metapneumovirus');
  assert.equal(a.length, 13350);
  const f = a.features.find((x) => x.kind === 'CDS' && x.name === 'F');
  assert.equal(f.start, 3052);
  assert.equal(f.end, 4671);
  assert.equal(f.product, 'fusion protein');
  assert.ok(!a.features.some((x) => x.kind === 'source'));
  assert.ok(a.features.every((x) => x.name));
});

test('a repeated, wrapped qualifier does not leak into the location', () => {
  const gb = `LOCUS       X1    100 bp    DNA   linear   VRL 01-JAN-2000
VERSION     X1.1
  ORGANISM  test virus
FEATURES             Location/Qualifiers
     CDS             1..20
                     /inference="first"
                     /inference="second that wraps
                     onto another line"
                     /gene="W"
ORIGIN
//
`;
  const a = parseGenbank(gb);
  assert.deepEqual(a.features[0].intervals, [[1, 20]]);
  assert.equal(a.features[0].name, 'W');
});

test('a GenBank join that wraps onto a second line is read whole', () => {
  const gb = `LOCUS       X1    100 bp    DNA   linear   VRL 01-JAN-2000
VERSION     X1.1
  ORGANISM  test virus
FEATURES             Location/Qualifiers
     CDS             join(complement(60..80),
                     complement(1..20))
                     /gene="W"
ORIGIN
//
`;
  const [f] = parseGenbank(gb).features;
  assert.equal(f.strand, '-');
  assert.deepEqual(f.intervals, [[1, 20], [60, 80]]);
  assert.equal(f.name, 'W');
});

test('a record over the cap is refused', async () => {
  assert.equal(MAX_RECORD_LENGTH, 500_000);
  await assert.rejects(fetchAnnotation('NC_000000.1', { fetch: okFetch(bigRecordXml), storage: null }),
    /500,?000.*gene record/s);
});

test('esummary refuses an oversize record before any efetch', async () => {
  const urls = [];
  const fetch = route(async (u) => { urls.push(u); return { ok: true, text: async () => hmpvXml }; },
    5_000_000);
  await assert.rejects(fetchAnnotation('NC_000001.11', { fetch, storage: null }),
    /5,000,000.*500,000.*gene record/s);
  assert.deepEqual(urls, []);
});

test('an esummary that fails still lets efetch proceed to the backstop', async () => {
  const fetch = route(async () => ({ ok: true, text: async () => bigRecordXml }), undefined);
  await assert.rejects(fetchAnnotation('NC_000000.1', { fetch, storage: null }), /500,000/);
  const throwing = async (u) => {
    if (u.includes('esummary')) throw new Error('offline');
    return { ok: true, text: async () => hmpvXml };
  };
  assert.ok(await fetchAnnotation('NC_039199.1', { fetch: throwing, storage: null }));
});

test('an esummary that omits the length is ignored', async () => {
  const fetch = async (u) => (u.includes('esummary')
    ? { ok: true, json: async () => ({ result: { uids: ['1'], 1: {} } }) }
    : { ok: true, text: async () => hmpvXml });
  assert.ok(await fetchAnnotation('NC_039199.1', { fetch, storage: null }));
});

test('the efetch request goes to the accession', async () => {
  const urls = [];
  const fetch = route(async (u) => { urls.push(u); return { ok: true, text: async () => hmpvXml }; });
  await fetchAnnotation('NC_039199.1', { fetch, storage: null });
  assert.deepEqual(urls, ['https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    + '?db=nuccore&rettype=gbc&retmode=xml&id=NC_039199.1']);
});

test('the esummary request asks for json', async () => {
  let url;
  const fetch = async (u) => { url = u; return { ok: true, json: async () => ({ result: { uids: ['7'], 7: { slen: 1453 } } }) }; };
  assert.equal(await fetchRecordLength('M54878.1', { fetch }), 1453);
  assert.equal(url, 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    + '?db=nuccore&retmode=json&id=M54878.1');
});

test('a second request is served from the cache, whatever the accession case', async () => {
  const data = new Map();
  const storage = { getItem: (k) => data.get(k) ?? null, setItem: (k, v) => data.set(k, v) };
  let calls = 0;
  const fetch = route(async () => { calls++; return { ok: true, text: async () => hmpvXml }; });
  const first = await fetchAnnotation('nc_039199.1', { fetch, storage });
  const second = await fetchAnnotation('NC_039199.1', { fetch, storage });
  assert.equal(calls, 1);
  assert.deepEqual(second, first);
  assert.equal(data.size, 1);
  assert.ok([...data.keys()][0].endsWith('NC_039199.1'));
});

test('a corrupt cache entry is ignored', async () => {
  const storage = { getItem: () => '{not json', setItem() {} };
  assert.equal((await fetchAnnotation('NC_039199.1', { fetch: okFetch(hmpvXml), storage })).organism,
    'human metapneumovirus');
});

test('storage that throws does not break the fetch', async () => {
  const storage = { getItem() { throw new Error('private'); },
                    setItem() { throw new Error('private'); } };
  assert.ok(await fetchAnnotation('NC_039199.1', { fetch: okFetch(hmpvXml), storage }));
});

test('a 429 and a 503 are retried after 1 s and 3 s', async () => {
  const replies = [429, 503, 200];
  const fetch = route(async () => {
    const status = replies.shift();
    return { ok: status === 200, status, text: async () => hmpvXml };
  }, 13350);
  const waits = [];
  const a = await fetchAnnotation('NC_039199.1',
    { fetch, storage: null, wait: async (ms) => { waits.push(ms); } });
  assert.equal(a.organism, 'human metapneumovirus');
  assert.deepEqual(waits, [1000, 3000]);
});

test('a persistent 500 gives up after two retries', async () => {
  let calls = 0;
  const fetch = route(async () => { calls++; return { ok: false, status: 500, text: async () => '' }; });
  await assert.rejects(fetchAnnotation('NC_039199.1', { fetch, storage: null, wait: noWait }), /500/);
  assert.equal(calls, 3);
});

test('a 404 is not retried', async () => {
  let calls = 0;
  const fetch = route(async () => { calls++; return { ok: false, status: 404, text: async () => '' }; });
  await assert.rejects(fetchAnnotation('NC_000000.1', { fetch, storage: null, wait: noWait }), /404/);
  assert.equal(calls, 1);
});

test('an unreadable location skips that feature and reports it', () => {
  const xml = hmpvXml.replaceAll('<INSDFeature_location>3052..4671</INSDFeature_location>',
    '<INSDFeature_location>garbage</INSDFeature_location>');
  const a = parseInsdseq(xml);
  assert.equal(a.skipped.length, 2);
  assert.deepEqual(a.skipped[0], { kind: 'gene', location: 'garbage' });
  assert.ok(!a.features.some((x) => x.name === 'F'));
  assert.ok(a.features.some((x) => x.name === 'N'));
  assert.deepEqual(parseInsdseq(hmpvXml).skipped, []);
});

test('a feature with no gene, locus_tag or product is named after its kind', () => {
  const xml = giardiaXml.replace(/<INSDQualifier>\s*<INSDQualifier_name>product<\/INSDQualifier_name>[\s\S]*?<\/INSDQualifier>/, '');
  assert.equal(parseInsdseq(xml).features[0].name, 'rRNA');
});

test('genes carry no product while their CDS siblings do', () => {
  const { features } = parseInsdseq(hmpvXml);
  assert.ok(features.filter((x) => x.kind === 'gene').every((x) => x.product === ''));
  assert.ok(features.filter((x) => x.kind === 'CDS').every((x) => x.product !== ''));
});

test('numeric character references are decoded', () => {
  const xml = giardiaXml.replace('16S-like ribosomal RNA</INSDQualifier_value>',
    'beta &#946; &#x3b2; &amp;amp</INSDQualifier_value>');
  assert.equal(parseInsdseq(xml).features[0].product, 'beta \u03b2 \u03b2 &amp');
});

test('a record hands over its own bases, which the anchor map needs', () => {
  const xml = hmpvXml.replace('</INSDSeq_feature-table>',
    '</INSDSeq_feature-table><INSDSeq_sequence>gtataaatta\ngattccaaaa</INSDSeq_sequence>');
  assert.equal(parseInsdseq(xml).sequence, 'GTATAAATTAGATTCCAAAA');
  assert.equal(parseInsdseq(hmpvXml).sequence, '');   // a record without them is not an error
});

test('a GenBank file hands over the bases under ORIGIN', () => {
  const gb = `LOCUS       X1    40 bp    DNA   linear   VRL 01-JAN-2000
VERSION     X1.1
  ORGANISM  test virus
FEATURES             Location/Qualifiers
     CDS             1..20
                     /gene="W"
ORIGIN      
        1 gtataaatta gattccaaaa
       21 aaatatggga caagtgaaaa
//
`;
  const a = parseGenbank(gb);
  assert.equal(a.sequence, 'GTATAAATTAGATTCCAAAAAAATATGGGACAAGTGAAAA');
  assert.equal(a.features[0].name, 'W');
  assert.equal(parseGenbank(fixture('hmpv.gb')).sequence, '');   // its ORIGIN block is empty
});
