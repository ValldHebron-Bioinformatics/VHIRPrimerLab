"""Candidate primer pair design with primer3 on the consensus sequence."""

from __future__ import annotations

import logging
from typing import Any

import primer3

from vhirprimerlab.structure import MAX_END_GC
from vhirprimerlab.thermo import SALT_CORRECTIONS, TM_FORMULA

from vhirprimerlab.models import (
    CandidatePair,
    Consensus,
    DesignParams,
    DesignResult,
    PipelineError,
    Primer,
)

log = logging.getLogger(__name__)

EXPLAIN_KEYS = ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN")


def _global_args(params: DesignParams) -> dict[str, Any]:
    return {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_MIN_TM": params.tm_min,
        "PRIMER_OPT_TM": params.tm_opt,
        "PRIMER_MAX_TM": params.tm_max,
        "PRIMER_PAIR_MAX_DIFF_TM": params.max_tm_diff,
        "PRIMER_MIN_SIZE": params.len_min,
        "PRIMER_OPT_SIZE": (params.len_min + params.len_max) // 2,
        "PRIMER_MAX_SIZE": params.len_max,
        "PRIMER_MIN_GC": params.gc_min,
        "PRIMER_MAX_GC": params.gc_max,
        "PRIMER_PRODUCT_SIZE_RANGE": [[params.amplicon_min, params.amplicon_max]],
        "PRIMER_NUM_RETURN": params.candidates,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_END_GC": MAX_END_GC,  # primer3's default (5) sets no limit
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_SALT_MONOVALENT": params.thermo.mv_conc,
        "PRIMER_SALT_DIVALENT": params.thermo.dv_conc,
        "PRIMER_DNTP_CONC": params.thermo.dntp_conc,
        "PRIMER_DNA_CONC": params.thermo.dna_conc,
        "PRIMER_TM_FORMULA": TM_FORMULA[params.thermo.tm_method],
        "PRIMER_SALT_CORRECTIONS": SALT_CORRECTIONS[params.thermo.salt_method],
    }


def run_primer3(
    consensus: Consensus, included: tuple[int, int], params: DesignParams
) -> DesignResult:
    """Ask primer3 for up to params.candidates pairs inside the included region."""
    seq_args = {
        "SEQUENCE_ID": "consensus",
        "SEQUENCE_TEMPLATE": consensus.seq,
        "SEQUENCE_INCLUDED_REGION": list(included),
    }
    try:
        res = primer3.bindings.design_primers(seq_args, _global_args(params))
    except (OSError, ValueError) as exc:
        raise PipelineError(f"primer3 failed: {exc}") from exc

    pairs: list[CandidatePair] = []
    for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
        lstart, llen = res[f"PRIMER_LEFT_{i}"]
        r3, rlen = res[f"PRIMER_RIGHT_{i}"]
        fwd = Primer(
            seq=res[f"PRIMER_LEFT_{i}_SEQUENCE"], start=lstart, length=llen,
            tm=res[f"PRIMER_LEFT_{i}_TM"], gc=res[f"PRIMER_LEFT_{i}_GC_PERCENT"],
        )
        rev = Primer(
            seq=res[f"PRIMER_RIGHT_{i}_SEQUENCE"], start=r3 - rlen + 1, length=rlen,
            tm=res[f"PRIMER_RIGHT_{i}_TM"], gc=res[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
        )
        pairs.append(CandidatePair(fwd, rev, res[f"PRIMER_PAIR_{i}_PENALTY"]))

    explain = {k: str(res[k]) for k in EXPLAIN_KEYS if k in res}
    log.debug("primer3 returned %d pairs; %s", len(pairs), explain)
    return DesignResult(pairs, explain)
