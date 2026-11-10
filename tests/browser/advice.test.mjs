// node --test tests/browser/advice.test.mjs
// What the page says about edge advice and which limits a re-run sends. The limits are the part
// most likely to be wrong, so they are tested here directly, not through the page.
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  adviceLimits, adviceSentence, isAdviceWarning, rerunJob, rerunLabel, sharedNote,
} from '../../vhirprimerlab/web/static/advice.js';

const F_START = { span: 'F', end: 'start', limit: 159, extra_loss: 59, edge_before: 43.01075268817204, edge_after: 97.84946236559139 };
const G_START = { span: 'G', end: 'start', limit: 120, extra_loss: 20, edge_before: 60, edge_after: 90 };
const G_END = { span: 'G', end: 'end', limit: 140, extra_loss: 40, edge_before: 55.5, edge_after: 99 };

test('the sentence is composed from the fields, with percentages as the results show them', () => {
  assert.equal(adviceSentence(F_START),
    'The first amplicon of F covers 43 % of your sequences. Giving up 59 more bases at its start would reach 97.8 %.');
});

test('the last amplicon is told by its end, and a whole-genome run has no span name', () => {
  assert.equal(adviceSentence({ ...G_END, span: null, extra_loss: 1, edge_before: 12.5, edge_after: 80 }),
    'The last amplicon covers 12.5 % of your sequences. Giving up 1 more base at its end would reach 80 %.');
});

test('a percentage that rounds up to 100 is not shown as full coverage', () => {
  assert.match(adviceSentence({ ...F_START, edge_after: 99.96 }), /would reach 99\.9 %\.$/);
  assert.match(adviceSentence({ ...F_START, edge_after: 100 }), /would reach 100 %\.$/);
});

test('no advice sends no limits', () => {
  assert.deepEqual(adviceLimits([]), {});
  assert.deepEqual(adviceLimits(undefined), {});
});

test('one advice raises only its own end', () => {
  assert.deepEqual(adviceLimits([F_START]), { max_start_loss: 159 });
  assert.deepEqual(adviceLimits([G_END]), { max_end_loss: 140 });
});

test('two advices at the same end take the larger limit, whatever the order', () => {
  assert.deepEqual(adviceLimits([F_START, G_START]), { max_start_loss: 159 });
  assert.deepEqual(adviceLimits([G_START, F_START]), { max_start_loss: 159 });
});

test('advices at different ends raise each end on its own', () => {
  assert.deepEqual(adviceLimits([F_START, G_END]), { max_start_loss: 159, max_end_loss: 140 });
});

test('malformed advice is ignored instead of sending NaN', () => {
  assert.deepEqual(adviceLimits([{ end: 'start', limit: 'x' }, { end: 'middle', limit: 5 }, null]), {});
});

const JOB = {
  mode: 'genes', fasta: '>a\nAC\n', primers: '',
  settings: { tile_regions: ['F:1-9'], max_start_loss: 100, max_end_loss: 100, overlap: 50 },
};

test('a re-run keeps everything and raises only the advised end', () => {
  const next = rerunJob(JOB, [F_START]);
  assert.deepEqual(next.settings, { ...JOB.settings, max_start_loss: 159 });
  assert.equal(next.fasta, JOB.fasta);
  assert.equal(next.mode, 'genes');
});

test('a re-run does not touch the job it was made from', () => {
  rerunJob(JOB, [F_START, G_END]);
  assert.equal(JOB.settings.max_start_loss, 100);
  assert.equal(JOB.settings.max_end_loss, 100);
});

test('a re-run never lowers a limit the run already had', () => {
  const next = rerunJob({ ...JOB, settings: { ...JOB.settings, max_start_loss: 300 } }, [F_START]);
  assert.equal(next.settings.max_start_loss, 300);
});

test('in tile mode, where the limits may be absent, the advised one is added', () => {
  const next = rerunJob({ mode: 'tile', settings: { amplicons: 10 } }, [F_START, G_START]);
  assert.deepEqual(next.settings, { amplicons: 10, max_start_loss: 159 });
});

test('the button names what it gives up', () => {
  assert.equal(rerunLabel([F_START]), 'Re-run, giving up 59 more bases at the start');
  assert.equal(rerunLabel([F_START, G_END]), 'Re-run with these allowances');
});

test('several advices, or several spans, say the allowance is for every span', () => {
  assert.equal(sharedNote([F_START], 1), '');
  assert.match(sharedNote([F_START], 3), /every span/);
  const several = sharedNote([F_START, G_START, G_END], 2);
  assert.match(several, /every span/);
  assert.match(several, /159 bases at the start/);
  assert.match(several, /140 at the end/);
});

test('only the engine\'s own advice line is recognised as a warning to replace', () => {
  assert.equal(isAdviceWarning('The first amplicon of F covers 43.0 % of your sequences; giving up 59 more bases at its start would reach 97.8 %: add --max-start-loss 159'), true);
  assert.equal(isAdviceWarning('The last amplicon of the genome covers 4.0 % of your sequences; giving up 3 more bases at its end would reach 9.0 %: add --max-end-loss 53'), true);
  assert.equal(isAdviceWarning('Amplicons below --min-coverage: 3'), false);
});
