/* Read a reference annotation (genes, proteins, UTRs) from NCBI or from a GenBank file.
   This is the only module that talks to the network, and only from the browser: the Python
   engine never does. */

export const MAX_RECORD_LENGTH = 500_000;

const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
// The version is part of the key so a change to what is kept of a record (the kinds, the
// sequence) cannot be masked by old entries.
const CACHE_ROOT = 'vhirprimerlab.annotation.';
const CACHE_PREFIX = `${CACHE_ROOT}v3.`;
// A record's bases dominate the size of an entry, and the quota is shared with the rest of the
// app (the Basic/Advanced preference is in the same storage), so a big record is re-fetched
// instead of cached.
const CACHE_MAX_BASES = 100_000;
const RETRY_WAITS = [1000, 3000];   // ms before the first and second retry

// What a user can click to ask for primers. rRNA and the other RNAs are here because a
// parasite record is a marker locus (18S, 16S) that carries nothing else.
const KEPT_KINDS = new Set([
  'gene', 'CDS', 'mat_peptide', 'sig_peptide', "5'UTR", "3'UTR", 'stem_loop', 'misc_feature',
  'rRNA', 'tRNA', 'ncRNA',
]);

// ---- locations --------------------------------------------------------------------------

// Split on the commas that are not inside parentheses.
function splitTop(text) {
  const parts = [];
  let depth = 0;
  let from = 0;
  for (let i = 0; i < text.length; i++) {
    if (text[i] === '(') depth++;
    else if (text[i] === ')') depth--;
    else if (text[i] === ',' && depth === 0) { parts.push(text.slice(from, i)); from = i + 1; }
  }
  parts.push(text.slice(from));
  return parts;
}

function parseExpr(text, whole) {
  const bad = (why) => new Error(`Cannot read the location "${whole}": ${why}`);
  const s = text.trim();
  const wrapped = /^(complement|join|order)\((.*)\)$/s.exec(s);
  if (wrapped) {
    const [, op, inner] = wrapped;
    if (op === 'complement') {
      const sub = parseExpr(inner, whole);
      return { strand: sub.strand === '+' ? '-' : '+', intervals: sub.intervals };
    }
    const subs = splitTop(inner).map((p) => parseExpr(p, whole));
    if (subs.some((p) => p.strand !== subs[0].strand)) throw bad('parts on both strands');
    return { strand: subs[0].strand, intervals: subs.flatMap((p) => p.intervals) };
  }
  // `<` and `>` mark an end that lies beyond the record; the position is still the bound.
  const range = /^<?(\d+)(?:\.\.>?(\d+))?$/.exec(s);
  if (!range) throw bad('unrecognised syntax');
  const start = Number(range[1]);
  const end = range[2] === undefined ? start : Number(range[2]);
  if (start < 1 || end < start) throw bad('a range must run upward from 1');
  return { strand: '+', intervals: [[start, end]] };
}

// A join that wraps the origin of a circular record, such as join(1000..1200,1..50), is read as
// plain intervals: its span is 1..1200, not the 1000..50 it covers on the molecule.
export function parseLocation(text) {
  const loc = parseExpr(String(text), text);
  loc.intervals.sort((a, b) => a[0] - b[0]);
  return loc;
}

// ---- features ---------------------------------------------------------------------------

// raw: [{ key, location, quals: { name: value } }] in record order. A feature whose location
// cannot be read is skipped and reported, so one odd feature does not lose the whole record.
function buildFeatures(raw) {
  const features = [];
  const skipped = [];
  for (const { key, location, quals } of raw) {
    if (!KEPT_KINDS.has(key)) continue;
    let loc;
    try { loc = parseLocation(location); } catch { skipped.push({ kind: key, location }); continue; }
    const { strand, intervals } = loc;
    const product = quals.product || '';
    features.push({
      kind: key,
      // Not every record names its genes: a marker locus has only a product.
      name: quals.gene || quals.locus_tag || product || key,
      product,
      strand,
      intervals,
      start: intervals[0][0],
      end: Math.max(...intervals.map((iv) => iv[1])),
    });
  }
  return { features, skipped };
}

const ENTITIES = { lt: '<', gt: '>', quot: '"', apos: "'", amp: '&' };
const unescapeXml = (s) => s.replace(/&(?:#(\d+)|#x([0-9a-fA-F]+)|(lt|gt|quot|apos|amp));/g,
  (_, dec, hex, name) => (name ? ENTITIES[name] : String.fromCodePoint(dec ? +dec : parseInt(hex, 16))));

// One reader for the page, a worker and Node, so none depends on DOMParser. NCBI's
// XML is machine-generated and regular (no attributes, no nesting of the tags read here).
function insdRaw(xmlText) {
  if (!/<INSDSeq>/.test(xmlText)) return null;
  const one = (src, tag) => {
    const m = new RegExp(`<${tag}>([\\s\\S]*?)</${tag}>`).exec(src);
    return m ? unescapeXml(m[1]).trim() : '';
  };
  const raw = [];
  for (const [, body] of xmlText.matchAll(/<INSDFeature>([\s\S]*?)<\/INSDFeature>/g)) {
    const quals = {};
    const qualRe = /<INSDQualifier_name>([\s\S]*?)<\/INSDQualifier_name>(?:\s*<INSDQualifier_value>([\s\S]*?)<\/INSDQualifier_value>)?/g;
    for (const [, name, value] of body.matchAll(qualRe)) {
      if (!Object.hasOwn(quals, name)) quals[name] = unescapeXml(value ?? '').trim();
    }
    raw.push({ key: one(body, 'INSDFeature_key'), location: one(body, 'INSDFeature_location'), quals });
  }
  // The header fields sit before the feature table; the source feature repeats "organism".
  const table = xmlText.indexOf('<INSDSeq_feature-table>');
  const head = table < 0 ? xmlText : xmlText.slice(0, table);
  return {
    accession: one(head, 'INSDSeq_accession-version'), organism: one(head, 'INSDSeq_organism'),
    length: Number(one(head, 'INSDSeq_length')),
    // The record's own bases, which the anchor map needs to place its features on the user's
    // reference. It sits after the feature table, and a record may come without it.
    sequence: one(xmlText, 'INSDSeq_sequence').replace(/[^A-Za-z]/g, '').toUpperCase(),
    raw,
  };
}

export function parseInsdseq(xmlText) {
  const rec = insdRaw(xmlText);
  if (!rec) throw new Error('NCBI returned no sequence record (is the accession right?)');
  return {
    accession: rec.accession, organism: rec.organism, length: rec.length,
    sequence: rec.sequence, ...buildFeatures(rec.raw),
  };
}

// The bases under ORIGIN, if the file carries them.
function genbankSequence(lines) {
  const bases = [];
  let inOrigin = false;
  for (const line of lines) {
    if (!inOrigin) inOrigin = /^ORIGIN\b/.test(line);
    else if (/^\/\//.test(line)) break;
    else bases.push(line.replace(/[^A-Za-z]/g, ''));
  }
  return bases.join('').toUpperCase();
}

export function parseGenbank(text) {
  const lines = text.split(/\r?\n/);
  let accession = '';
  let organism = '';
  let length = 0;
  const raw = [];
  let inFeatures = false;
  let feature = null;
  let qual = null;      // the qualifier being continued, if any
  let skipping = false; // inside a repeated qualifier, whose wrapped lines are dropped

  for (const line of lines) {
    if (!inFeatures) {
      let m;
      if ((m = /^LOCUS\s+\S+\s+(\d+)\s+bp/.exec(line))) length = Number(m[1]);
      else if ((m = /^VERSION\s+(\S+)/.exec(line))) accession = m[1];
      else if ((m = /^\s+ORGANISM\s+(.+?)\s*$/.exec(line))) organism = organism || m[1];
      else if (/^FEATURES\b/.test(line)) inFeatures = true;
      continue;
    }
    if (/^\S/.test(line)) break;                     // ORIGIN, //, CONTIG: the table is over
    let m;
    if ((m = /^ {5}(\S+)\s+(\S.*?)\s*$/.exec(line))) {
      feature = { key: m[1], location: m[2], quals: {} };
      raw.push(feature);
      qual = null;
      skipping = false;
    } else if (feature && (m = /^ {21}\/([\w-]+)(?:=(.*?))?\s*$/.exec(line))) {
      qual = m[1];
      skipping = Object.hasOwn(feature.quals, qual);  // keep the first of a repeated qualifier
      if (skipping) qual = null;
      else feature.quals[qual] = m[2] ?? '';
    } else if (feature && (m = /^ {21}(\S.*?)\s*$/.exec(line))) {
      // A long location (join) or a long qualifier value wraps onto the next line.
      if (qual) feature.quals[qual] += ` ${m[1]}`;
      else if (!skipping) feature.location += m[1];
    }
  }
  for (const f of raw) {
    for (const k of Object.keys(f.quals)) f.quals[k] = f.quals[k].replace(/^"|"$/g, '').trim();
  }
  if (!accession && !raw.length) throw new Error('This does not look like a GenBank file');
  return { accession, organism, length, sequence: genbankSequence(lines), ...buildFeatures(raw) };
}

// ---- fetching ---------------------------------------------------------------------------

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function defaultStorage() {
  try { return globalThis.localStorage ?? null; } catch { return null; }   // blocked storage throws
}

const swept = new WeakSet();   // one sweep per storage, however many records are read

// Entries written by an earlier version of this cache can never be read again, and they are the
// bulky ones, so they go the first time the cache is used.
function sweepOldEntries(storage) {
  try {
    if (!storage || swept.has(storage)) return;
    swept.add(storage);
    const stale = [];
    for (let i = 0; i < storage.length; i++) {
      const key = storage.key(i);
      if (key && key.startsWith(CACHE_ROOT) && !key.startsWith(CACHE_PREFIX)) stale.push(key);
    }
    stale.forEach((key) => storage.removeItem(key));
  } catch { /* the cache is optional */ }
}

function readCache(storage, key) {
  try {
    sweepOldEntries(storage);
    const hit = JSON.parse(storage?.getItem(key));
    return hit && Array.isArray(hit.features) ? hit : null;
  } catch { return null; }
}

function writeCache(storage, key, annotation) {
  // A cached record without its bases could not be mapped, so such a record is not cached at all.
  if (annotation.sequence.length > CACHE_MAX_BASES) return;
  try { storage?.setItem(key, JSON.stringify(annotation)); } catch { /* the cache is optional */ }
}

const tooLong = (accession, length) => new Error(
  `${accession} is ${length.toLocaleString('en-US')} nt, over the `
  + `${MAX_RECORD_LENGTH.toLocaleString('en-US')} nt limit. Use a gene record instead of a chromosome.`);

// GET with the two retries on a 429 or a 5xx; returns the last response whatever its status.
async function getWithRetry(fetch, url, wait) {
  for (let attempt = 0; ; attempt++) {
    const response = await fetch(url);
    const retryable = response.status === 429 || response.status >= 500;
    if (response.ok || !retryable || attempt >= RETRY_WAITS.length) return response;
    await wait(RETRY_WAITS[attempt]);
  }
}

// The length of a record, from the small esummary reply, so a chromosome is refused before its
// sequence is downloaded. null when NCBI cannot say; the caller then relies on its own checks.
export async function fetchRecordLength(accession, { fetch = globalThis.fetch, wait = sleep } = {}) {
  try {
    const url = `${EUTILS}/esummary.fcgi?db=nuccore&retmode=json&id=${encodeURIComponent(accession)}`;
    const response = await getWithRetry(fetch, url, wait);
    if (!response.ok) return null;
    const { result } = await response.json();
    const length = Number(result?.[result.uids?.[0]]?.slen);
    return Number.isFinite(length) && length > 0 ? length : null;
  } catch { return null; }
}

// `wait` exists so the tests need not sit through the backoff.
export async function fetchAnnotation(accession, {
  fetch = globalThis.fetch, storage = defaultStorage(), wait = sleep,
} = {}) {
  const key = CACHE_PREFIX + accession.toUpperCase();
  const cached = readCache(storage, key);
  if (cached) return cached;

  const known = await fetchRecordLength(accession, { fetch, wait });
  if (known > MAX_RECORD_LENGTH) throw tooLong(accession, known);

  const url = `${EUTILS}/efetch.fcgi?db=nuccore&rettype=gbc&retmode=xml&id=${encodeURIComponent(accession)}`;
  const response = await getWithRetry(fetch, url, wait);
  if (!response.ok) throw new Error(`NCBI answered ${response.status} for ${accession}`);

  const annotation = parseInsdseq(await response.text());
  // Backstop for when esummary was unavailable or left the length out.
  if (annotation.length > MAX_RECORD_LENGTH) throw tooLong(accession, annotation.length);
  writeCache(storage, key, annotation);
  return annotation;
}
