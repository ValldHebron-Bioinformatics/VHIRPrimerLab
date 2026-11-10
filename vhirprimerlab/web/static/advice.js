/* Edge advice: what the engine says giving up a few more bases at a span edge would buy, put into
   words, and the settings a re-run with that allowance is started with. Pure, so it is tested
   without a page. The engine only advises; the person decides, so nothing here runs by itself. */

// Percentages as the results show them (one decimal only when there is one), except that a value
// that merely rounds up to 100 stays 99.9: "100 %" would promise every sequence.
function pct(n) {
  const x = Number(n);
  const shown = x < 100 && x >= 99.95 ? 99.9 : Number(x.toFixed(1));
  return `${shown} %`;
}

const word = (end) => (end === 'start' ? 'first' : 'last');

export function adviceSentence(a) {
  const of = a.span ? ` of ${a.span}` : '';
  return `The ${word(a.end)} amplicon${of} covers ${pct(a.edge_before)} of your sequences. `
    + `Giving up ${a.extra_loss} more ${a.extra_loss === 1 ? 'base' : 'bases'} at its ${a.end} `
    + `would reach ${pct(a.edge_after)}.`;
}

// The engine's own warning for the same advice. The page shows the sentence above instead, so it
// must not also list this one, which ends in a command-line flag.
export function isAdviceWarning(w) {
  return /^The (first|last) amplicon .* add --max-(start|end)-loss \d+$/.test(w);
}

// The engine takes one limit per end for the whole run, so the largest advised for an end is the
// only one that satisfies every span that asked about that end.
export function adviceLimits(advice) {
  const out = {};
  (advice || []).forEach((a) => {
    if (!a || (a.end !== 'start' && a.end !== 'end') || !Number.isFinite(Number(a.limit))) return;
    const key = `max_${a.end}_loss`;
    out[key] = Math.max(out[key] ?? -Infinity, Number(a.limit));
  });
  return out;
}

// The job again, with the advised ends raised. A limit the run already had is never lowered, and
// the job it was made from is left as it was (Retry still needs it).
export function rerunJob(job, advice) {
  const settings = { ...job.settings };
  Object.entries(adviceLimits(advice)).forEach(([key, limit]) => {
    settings[key] = Math.max(Number(settings[key]) || 0, limit);
  });
  return { ...job, settings };
}

export function rerunLabel(advice) {
  if (advice.length !== 1) return 'Re-run with these allowances';
  const { extra_loss: n, end } = advice[0];
  return `Re-run, giving up ${n} more ${n === 1 ? 'base' : 'bases'} at the ${end}`;
}

// Said whenever the allowance reaches spans that did not ask for it: several advices, or one
// advice in a run of several spans. Empty when the run has one span and one advice.
export function sharedNote(advice, spans) {
  if (advice.length < 2 && spans < 2) return '';
  const limits = adviceLimits(advice);
  const parts = [];
  if ('max_start_loss' in limits) parts.push(`${limits.max_start_loss} bases at the start`);
  if ('max_end_loss' in limits) parts.push(`${limits.max_end_loss} at the end`);
  return 'The design takes one allowance for each end, so it applies to every span, not only the '
    + `ones listed: up to ${parts.join(' and ')} of any span may fall outside its amplicons.`;
}
