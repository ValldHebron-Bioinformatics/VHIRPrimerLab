from __future__ import annotations

from tests.synthetic import make_genome_family
from vhirprimerlab.tiling import geometry, max_amplicons, target_length


def test_make_genome_family_shape() -> None:
    rows = make_genome_family()
    assert len(rows) == 12 and {len(s) for _, s in rows} == {6000}
    assert rows[0][0] == "g01" and set("".join(s for _, s in rows)) <= set("ACGT")


def test_target_length() -> None:
    assert target_length(6000, 5, 50) == 1240
    assert target_length(1000, 1, 50) == 1000


def test_geometry_values() -> None:
    g = geometry(6000, 5, 50, 0.20, 25)
    assert g is not None
    assert (g.target_len, g.len_min, g.len_max) == (1240, 992, 1488)
    assert g.max_outer_overlap == 495
    assert g.tolerance == 0.20


def test_geometry_infeasible_when_amplicons_too_short() -> None:
    assert geometry(6000, 40, 50, 0.20, 25) is None  # L = 199


def test_max_amplicons() -> None:
    n = max_amplicons(6000, 50, 25)
    assert n > 1
    assert geometry(6000, n, 50, 0.20, 25) is not None
    assert geometry(6000, n + 1, 50, 0.20, 25) is None
