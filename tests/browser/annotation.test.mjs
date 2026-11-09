// node --test tests/browser/annotation.test.mjs
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import {
  parseLocation, parseInsdseq, parseGenbank, fetchAnnotation, MAX_RECORD_LENGTH,
} from '../../vhirprimerlab/web/static/annotation.js';

const fixture = (name) => readFileSync(new URL(`./fixtures/${name}`, import.meta.url), 'utf8');
const hmpvXml = fixture('NC_039199.xml');
const giardiaXml = fixture('parasite-gene.xml');
const bigRecordXml = hmpvXml.replace('<INSDSeq_length>13350<', '<INSDSeq_length>600000<');
const okFetch = (body) => async () => ({ ok: true, status: 200, text: async () => body });
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
  await assert.rejects(fetchAnnotation('NC_000000.1', { fetch: okFetch(bigRecordXml) }),
    /500,?000.*gene record/s);
});

test('the request goes to efetch for the accession', async () => {
  let url;
  const fetch = async (u) => { url = u; return { ok: true, text: async () => hmpvXml }; };
  await fetchAnnotation('NC_039199.1', { fetch });
  assert.equal(url, 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    + '?db=nuccore&rettype=gbc&retmode=xml&id=NC_039199.1');
});

test('a second request is served from the cache', async () => {
  const data = new Map();
  const storage = { getItem: (k) => data.get(k) ?? null, setItem: (k, v) => data.set(k, v) };
  let calls = 0;
  const fetch = async () => { calls++; return { ok: true, text: async () => hmpvXml }; };
  const first = await fetchAnnotation('NC_039199.1', { fetch, storage });
  const second = await fetchAnnotation('NC_039199.1', { fetch, storage });
  assert.equal(calls, 1);
  assert.deepEqual(second, first);
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
  const fetch = async () => {
    const status = replies.shift();
    return { ok: status === 200, status, text: async () => hmpvXml };
  };
  const waits = [];
  const a = await fetchAnnotation('NC_039199.1', { fetch, wait: async (ms) => { waits.push(ms); } });
  assert.equal(a.organism, 'human metapneumovirus');
  assert.deepEqual(waits, [1000, 3000]);
});

test('a persistent 500 gives up after two retries', async () => {
  let calls = 0;
  const fetch = async () => { calls++; return { ok: false, status: 500, text: async () => '' }; };
  await assert.rejects(fetchAnnotation('NC_039199.1', { fetch, wait: noWait }), /500/);
  assert.equal(calls, 3);
});

test('a 404 is not retried', async () => {
  let calls = 0;
  const fetch = async () => { calls++; return { ok: false, status: 404, text: async () => '' }; };
  await assert.rejects(fetchAnnotation('NC_000000.1', { fetch, wait: noWait }), /404/);
  assert.equal(calls, 1);
});
