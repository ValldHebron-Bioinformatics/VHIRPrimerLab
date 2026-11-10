from __future__ import annotations

from tests.synthetic import make_genome_family
from vhirprimerlab.tiling import amplicons_for, geometry, max_amplicons, target_length


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


def test_amplicons_for() -> None:
    assert amplicons_for(400, 400, 50) == 1
    assert amplicons_for(150, 400, 50) == 1          # shorter than one amplicon
    assert amplicons_for(750, 400, 50) == 2          # (750 - 50) / 350 = 2
    assert amplicons_for(751, 400, 50) == 3
    assert amplicons_for(6018, 400, 50) == 18        # the hMPV L gene


def test_geometry_allows_one_short_amplicon() -> None:
    # 150 bp holds two 25 bp primers and an insert; with n == 1 there is no neighbour,
    # so the overlap condition must not apply
    g = geometry(150, 1, 50, 0.20, 25)
    assert g is not None and g.len_min == 120
    assert geometry(40, 1, 50, 0.20, 25) is None     # too short for two primers


def test_geometry_still_refuses_short_amplicons_with_neighbours() -> None:
    assert geometry(300, 2, 50, 0.20, 25) is None
