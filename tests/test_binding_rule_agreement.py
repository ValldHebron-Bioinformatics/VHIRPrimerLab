"""Coverage scoring and the specificity scan must call the same thing a binding site.

"Does this primer bind here" is implemented twice: scoring.match_primer decides coverage of the
aligned sequences, and specificity._Scanner finds secondary sites in raw ones. They are written
differently (one walks alignment columns, the other uses numpy seed codes), so a change to either
can make the report say a primer covers a sequence at a site the scan never sees, or the reverse.
The documented difference is N and other ambiguity codes in the sequence: coverage forgives a few,
the scan never forms a site on one. Everything else must agree, so these trials use plain A/C/G/T
targets.
"""

from __future__ import annotations

import random

import pytest

from vhirprimerlab.models import ScoringParams
from vhirprimerlab.scoring import match_primer
from vhirprimerlab.specificity import _Coded, _Scanner

FLANK = 20
PARAMS = ScoringParams(max_mismatches=2, window=5)


def _trial(rng: random.Random) -> tuple[str, str]:
    """A 20-mer primer, one base of it perhaps degenerate away from the 3' end, and a target that
    carries it with 0 to 4 substitutions among plain bases."""
    bases = [rng.choice("ACGT") for _ in range(20)]
    primer = list(bases)
    if rng.random() < 0.5:
        i = rng.randrange(0, 15)
        primer[i] = {"A": "R", "G": "R", "C": "Y", "T": "Y"}[bases[i]]
    target = list(bases)
    for i in rng.sample(range(20), rng.choice([0, 0, 1, 2, 2, 3, 4])):
        target[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    flank = "".join(rng.choice("ACGT") for _ in range(FLANK))
    return "".join(primer), flank + "".join(target) + flank[::-1]


def test_the_scan_and_the_scorer_agree_on_every_site() -> None:
    rng = random.Random(20260519)
    for n in range(300):
        primer, row = _trial(rng)
        scored = match_primer(primer, list(range(FLANK, FLANK + 20)), row, "right", PARAMS)
        scanner = _Scanner(primer, PARAMS.window, PARAMS.max_mismatches)
        found = dict(scanner.scan(_Coded(row)))
        also_by_text = dict(scanner.scan_text(row))
        assert found == also_by_text, f"trial {n}: numpy scan and text scan differ for {primer} / {row}"
        if scored.ok:
            assert found.get(FLANK) == scored.mismatches, (
                f"trial {n}: scored as covered ({scored.mismatches} mismatches) but the scan found "
                f"{found.get(FLANK)}: {primer} / {row[FLANK:FLANK + 20]}")
        else:
            assert FLANK not in found, (
                f"trial {n}: scored as not covered ({scored.reason}) but the scan reports a site: "
                f"{primer} / {row[FLANK:FLANK + 20]}")


@pytest.mark.parametrize("code", ["N", "R"])
def test_an_ambiguity_code_in_the_sequence_costs_the_scan_what_coverage_forgives(code: str) -> None:
    """The documented difference: coverage forgives up to max_ambiguous unknown bases (N, or R
    where the primer has A or G), the scan counts each as a full mismatch. Both see the same site;
    the scan's count is the scorer's mismatches plus the bases it forgave."""
    primer = "ACGTACGTACGTACGTACGT"
    row = "GATC" * 5 + primer[:8] + code + primer[9:] + "CATG" * 5   # position 9 of the primer is 'C'
    scored = match_primer(primer, list(range(FLANK, FLANK + 20)), row, "right", PARAMS)
    found = dict(_Scanner(primer, PARAMS.window, PARAMS.max_mismatches).scan(_Coded(row)))
    assert scored.ok and found[FLANK] == scored.mismatches + scored.ambiguous == 1


def test_an_ambiguity_code_in_the_three_prime_window_forms_no_site_and_is_not_covered() -> None:
    primer = "ACGTACGTACGTACGTACGT"
    row = "GATC" * 5 + primer[:18] + "N" + primer[19:] + "CATG" * 5
    scored = match_primer(primer, list(range(FLANK, FLANK + 20)), row, "right", PARAMS)
    found = dict(_Scanner(primer, PARAMS.window, PARAMS.max_mismatches).scan(_Coded(row)))
    assert not scored.ok and FLANK not in found
