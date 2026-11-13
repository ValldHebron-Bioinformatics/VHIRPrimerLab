/* Which candidate record fits the user's reference best. Each candidate's sequence is mapped onto
   the reference with the shared k-mers (anchors.js); the map that spans the largest share of its
   record wins. Pure apart from the two functions it is handed, so it is tested without a page. */

const baseOf = (accession) => accession.trim().toUpperCase();

// NCBI refuses more than about three requests a second, and each record costs two, so the
// records are fetched one after another with `pause` ms between them (tests pass 0).
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

// A record whose map is not trusted never wins: it would place genes at invented positions.
export async function findBestMatch(accessions, reference, { fetchAnnotation, buildAnchorMap, pause = 400, wait = sleep }) {
  const seen = new Set();
  const ranked = [];
  const failed = [];
  const noSequence = [];
  for (const raw of accessions) {
    const accession = String(raw).trim();
    const key = baseOf(accession);
    if (!accession || seen.has(key)) continue;
    if (seen.size && pause) await wait(pause);
    seen.add(key);
    let annotation;
    try {
      annotation = await fetchAnnotation(accession);
    } catch {
      failed.push(accession);
      continue;
    }
    if (!annotation.sequence) {
      noSequence.push(accession);
      continue;
    }
    const map = buildAnchorMap(annotation.sequence, reference);
    if (map.trusted) ranked.push({ accession, annotation, map });
  }
  ranked.sort((a, b) => b.map.covered - a.map.covered || b.map.count - a.map.count);
  return { best: ranked[0] ?? null, ranked, failed, noSequence };
}
