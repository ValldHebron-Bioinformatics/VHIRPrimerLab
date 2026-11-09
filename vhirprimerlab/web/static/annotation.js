/* Read a reference annotation (genes, proteins, UTRs) from NCBI or from a GenBank file.
   This is the only module that talks to the network, and only from the browser: the Python
   engine never does. */

export const MAX_RECORD_LENGTH = 500_000;

const EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
const CACHE_PREFIX = 'vhirprimerlab.annotation.';
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

export function parseLocation(text) {
  const loc = parseExpr(String(text), text);
  loc.intervals.sort((a, b) => a[0] - b[0]);
  return loc;
}

// ---- features ---------------------------------------------------------------------------

// raw: [{ key, location, quals: { name: value } }] in record order.
function buildFeatures(raw) {
  const features = [];
  for (const { key, location, quals } of raw) {
    if (!KEPT_KINDS.has(key)) continue;
    const { strand, intervals } = parseLocation(location);
    const product = quals.product || '';
    features.push({
      kind: key,
      // Not every record names its genes: a marker locus has only a product.
      name: quals.gene || quals.locus_tag || product,
      product,
      strand,
      intervals,
      start: intervals[0][0],
      end: Math.max(...intervals.map((iv) => iv[1])),
    });
  }
  return features;
}

const ENTITIES = { lt: '<', gt: '>', quot: '"', apos: "'", amp: '&' };
const unescapeXml = (s) => s.replace(/&(lt|gt|quot|apos|amp);/g, (_, e) => ENTITIES[e]);

function insdRawDom(xmlText) {
  const doc = new DOMParser().parseFromString(xmlText, 'application/xml');
  const seq = doc.getElementsByTagName('INSDSeq')[0];
  if (!seq) return null;
  const text = (tag) => (seq.getElementsByTagName(tag)[0]?.textContent ?? '').trim();
  const raw = [...seq.getElementsByTagName('INSDFeature')].map((f) => {
    const quals = {};
    for (const q of f.getElementsByTagName('INSDQualifier')) {
      const name = q.getElementsByTagName('INSDQualifier_name')[0]?.textContent;
      const value = q.getElementsByTagName('INSDQualifier_value')[0]?.textContent;
      if (name && !(name in quals)) quals[name] = (value ?? '').trim();
    }
    return {
      key: f.getElementsByTagName('INSDFeature_key')[0]?.textContent.trim(),
      location: f.getElementsByTagName('INSDFeature_location')[0]?.textContent.trim() ?? '',
      quals,
    };
  });
  return {
    accession: text('INSDSeq_accession-version'), organism: text('INSDSeq_organism'),
    length: Number(text('INSDSeq_length')), raw,
  };
}

// Node has no DOMParser, and the tests should not need a DOM library. NCBI's XML is regular
// enough (no attributes, no nesting of the tags read here) for a small reader.
function insdRawRegex(xmlText) {
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
      if (!(name in quals)) quals[name] = unescapeXml(value ?? '').trim();
    }
    raw.push({ key: one(body, 'INSDFeature_key'), location: one(body, 'INSDFeature_location'), quals });
  }
  // The header fields sit before the feature table; the source feature repeats "organism".
  const table = xmlText.indexOf('<INSDSeq_feature-table>');
  const head = table < 0 ? xmlText : xmlText.slice(0, table);
  return {
    accession: one(head, 'INSDSeq_accession-version'), organism: one(head, 'INSDSeq_organism'),
    length: Number(one(head, 'INSDSeq_length')), raw,
  };
}

export function parseInsdseq(xmlText) {
  const rec = typeof DOMParser === 'function' ? insdRawDom(xmlText) : insdRawRegex(xmlText);
  if (!rec) throw new Error('NCBI returned no sequence record (is the accession right?)');
  return {
    accession: rec.accession, organism: rec.organism, length: rec.length,
    features: buildFeatures(rec.raw),
  };
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
    } else if (feature && (m = /^ {21}\/([\w-]+)(?:=(.*?))?\s*$/.exec(line))) {
      qual = m[1];
      if (!(qual in feature.quals)) feature.quals[qual] = m[2] ?? '';
      else qual = null;                              // keep the first of a repeated qualifier
    } else if (feature && (m = /^ {21}(\S.*?)\s*$/.exec(line))) {
      // A long location (join) or a long qualifier value wraps onto the next line.
      if (qual) feature.quals[qual] += ` ${m[1]}`;
      else feature.location += m[1];
    }
  }
  for (const f of raw) {
    for (const k of Object.keys(f.quals)) f.quals[k] = f.quals[k].replace(/^"|"$/g, '').trim();
  }
  if (!accession && !raw.length) throw new Error('This does not look like a GenBank file');
  return { accession, organism, length, features: buildFeatures(raw) };
}

// ---- fetching ---------------------------------------------------------------------------

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function defaultStorage() {
  try { return globalThis.localStorage ?? null; } catch { return null; }   // blocked storage throws
}

function readCache(storage, key) {
  try {
    const hit = JSON.parse(storage?.getItem(key));
    return hit && Array.isArray(hit.features) ? hit : null;
  } catch { return null; }
}

function writeCache(storage, key, annotation) {
  try { storage?.setItem(key, JSON.stringify(annotation)); } catch { /* the cache is optional */ }
}

// `wait` exists so the tests need not sit through the backoff.
export async function fetchAnnotation(accession, {
  fetch = globalThis.fetch, storage = defaultStorage(), wait = sleep,
} = {}) {
  const key = CACHE_PREFIX + accession;
  const cached = readCache(storage, key);
  if (cached) return cached;

  const url = `${EFETCH}?db=nuccore&rettype=gbc&retmode=xml&id=${encodeURIComponent(accession)}`;
  let response;
  for (let attempt = 0; ; attempt++) {
    response = await fetch(url);
    const retryable = response.status === 429 || response.status >= 500;
    if (response.ok || !retryable || attempt >= RETRY_WAITS.length) break;
    await wait(RETRY_WAITS[attempt]);
  }
  if (!response.ok) throw new Error(`NCBI answered ${response.status} for ${accession}`);

  const annotation = parseInsdseq(await response.text());
  if (annotation.length > MAX_RECORD_LENGTH) {
    throw new Error(`${accession} is ${annotation.length.toLocaleString('en-US')} nt, over the `
      + `${MAX_RECORD_LENGTH.toLocaleString('en-US')} nt limit. Use a gene record `
      + 'instead of a chromosome.');
  }
  writeCache(storage, key, annotation);
  return annotation;
}
