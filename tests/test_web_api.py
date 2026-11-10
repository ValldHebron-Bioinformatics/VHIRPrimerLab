from __future__ import annotations

import io
import json
import time
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from tests.conftest import needs_mafft
from tests.synthetic import make_genome_family
from vhirprimerlab.web.app import create_app
from vhirprimerlab.web.params import Limits

GENOMES = "".join(f">{i}\n{s}\n" for i, s in make_genome_family())


def _post(client: TestClient, mode: str, params: dict, fasta: str = GENOMES, **extra):
    data = {"mode": mode, "params": json.dumps(params), "fasta_text": fasta, **extra}
    return client.post("/api/jobs", data=data)


def _wait(client: TestClient, job_id: str, timeout: float = 180) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        st = client.get(f"/api/jobs/{job_id}").json()
        if st["status"] in ("done", "failed"):
            return st
        time.sleep(0.3)
    raise AssertionError("job did not finish")


@pytest.fixture
def client(tmp_path: Path):
    with TestClient(create_app(tmp_path / "jobs", Limits(workers=2))) as c:
        yield c


def test_page_example_and_limits(client: TestClient) -> None:
    assert client.get("/").status_code == 200
    assert "<html" in client.get("/").text.lower()
    ex = client.get("/api/example")
    assert ex.status_code == 200 and ex.text.startswith(">g01")
    limits = client.get("/api/limits").json()
    assert limits["max_seqs"] == 500 and limits["max_amplicons"] == 200
    assert limits["retention_hours"] == 24


@needs_mafft
def test_tile_job_end_to_end(client: TestClient) -> None:
    r = _post(client, "tile", {"amplicons": 4})
    assert r.status_code == 200, r.text
    job = r.json()["id"]
    st = _wait(client, job)
    assert st["status"] == "done", st
    res = st["result"]
    assert res["mode"] == "tile" and len(res["rows"]) == 4
    assert res["reference_id"] == "g01" and res["reference_length"] == 6000
    tsv = client.get(f"/api/jobs/{job}/files/primers.tsv")
    assert tsv.status_code == 200 and tsv.text.startswith("name\tsequence\nVHIRPrimerLab_1_LEFT\t")
    z = zipfile.ZipFile(io.BytesIO(client.get(f"/api/jobs/{job}/zip").content))
    assert {"primers.tsv", "amplicons.tsv", "scheme.bed", "report.txt", "alignment.fasta"} <= set(z.namelist())


@needs_mafft
def test_check_job(client: TestClient) -> None:
    r = _post(client, "tile", {"amplicons": 2})
    rows = _wait(client, r.json()["id"])["result"]["rows"]
    primers = "".join(f"p{x['amplicon']}\t{x['forward_seq']}\t{x['reverse_seq']}\n" for x in rows)
    r = _post(client, "check", {}, primers_text=primers)
    assert r.status_code == 200, r.text
    st = _wait(client, r.json()["id"])
    assert st["status"] == "done" and [x["name"] for x in st["result"]["rows"]] == ["p1", "p2"]


def test_unknown_job_and_file_names(client: TestClient) -> None:
    assert client.get("/api/jobs/not-a-job").status_code == 404
    assert client.get(f"/api/jobs/{'a' * 32}").status_code == 404
    assert client.get(f"/api/jobs/{'a' * 32}/files/report.txt").status_code == 404
    assert client.get(f"/api/jobs/{'a' * 32}/files/..%2F..%2Fetc%2Fpasswd").status_code == 404
    assert client.get(f"/api/jobs/{'a' * 32}/zip").status_code == 404


@pytest.mark.parametrize("mode, params, extra, msg", [
    ("tile", {"amplicons": 500}, {}, "up to 200 amplicons"),
    ("tile", {"amplicons": 40}, {}, "use at most"),
    ("tile", {"amplicons": 5, "outdir": "/etc"}, {}, "unknown setting: outdir"),
    ("tile", {"amplicons": 5, "tm_min": 65}, {}, "--tm-min"),
    ("tile", {}, {}, "use -n N"),
    ("check", {}, {}, "primer"),
    ("check", {}, {"primers_text": "p1\tACGTXACGTACGTACGTA\tACGTACGTACGTACGTAC\n"}, "invalid character"),
])
def test_bad_requests(client: TestClient, mode, params, extra, msg) -> None:
    r = _post(client, mode, params, **extra)
    assert r.status_code == 400, r.text
    assert msg in r.json()["detail"]


def test_bad_fasta_and_too_many(tmp_path: Path) -> None:
    with TestClient(create_app(tmp_path / "jobs", Limits(max_seqs=5))) as c:
        r = _post(c, "tile", {"amplicons": 4})
        assert r.status_code == 400 and "Too many sequences: 12" in r.json()["detail"]
        r = _post(c, "tile", {"amplicons": 4}, fasta="just some text")
        assert r.status_code == 400
        r = _post(c, "tile", {"amplicons": 4}, fasta="")
        assert r.status_code == 400 and "sequences" in r.json()["detail"].lower()
        r = c.post("/api/jobs", data={"mode": "tile", "params": "{not json", "fasta_text": ">a\nACGT\n"})
        assert r.status_code == 400


def test_vhirprimerlab_web_help(capsys) -> None:
    from vhirprimerlab.cli import main

    with pytest.raises(SystemExit) as exc:
        main(["web", "--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "--port" in out and "--no-browser" in out


def test_primer_length_capped(client: TestClient) -> None:
    r = _post(client, "tile", {"amplicons": 1, "primer_len_max": 15000})
    assert r.status_code == 400 and "--primer-len-max" in r.json()["detail"]


def test_busy_queue_returns_503(tmp_path: Path) -> None:
    from tests.test_web_jobs import fake_job

    app = create_app(tmp_path / "jobs", Limits(workers=1, max_queue=1), target=fake_job)
    app.state.manager.start = lambda: None  # keep jobs queued
    with TestClient(app) as c:
        assert _post(c, "tile", {"amplicons": 4}).status_code == 200
        r = _post(c, "tile", {"amplicons": 4})
        assert r.status_code == 503 and "busy" in r.json()["detail"]


def test_oversized_body_rejected_before_parsing(tmp_path: Path) -> None:
    with TestClient(create_app(tmp_path / "jobs", Limits(max_upload_mb=1))) as c:
        big = ">a\n" + "A" * (3 * 1024 * 1024) + "\n"
        r = c.post("/api/jobs", files={"fasta": ("big.fa", big.encode(), "text/plain")},
                   data={"mode": "tile", "params": '{"amplicons": 2}'})
        assert r.status_code == 413


def test_only_installed_aligners(client: TestClient, monkeypatch) -> None:
    import shutil

    real = shutil.which
    monkeypatch.setattr(shutil, "which", lambda name: None if name == "muscle" else real(name))
    assert client.get("/api/limits").json()["aligners"] == ["mafft"]
    r = _post(client, "tile", {"amplicons": 4, "aligner": "muscle"})
    assert r.status_code == 400 and "MUSCLE is not installed" in r.json()["detail"]


def test_page_files_at_the_root(client: TestClient) -> None:
    for name in ("app.js", "engine.js", "app.css", "engine-worker.js", "favicon.svg"):
        assert client.get(f"/{name}").status_code == 200, name
    assert client.get("/nope.js").status_code == 404
    assert 'src="engine.js"' in client.get("/").text


def test_served_page_marks_the_server(client: TestClient) -> None:
    assert '<meta name="vhirprimerlab-server" content="1">' in client.get("/").text
    static = (Path(__file__).resolve().parent.parent / "vhirprimerlab/web/static/index.html").read_text()
    assert "vhirprimerlab-server" not in static  # on GitHub Pages the page runs in the browser


def test_delete_cancels_a_job(tmp_path: Path) -> None:
    from tests.test_web_jobs import fake_job

    app = create_app(tmp_path / "jobs", Limits(workers=1), target=fake_job)
    app.state.manager.start = lambda: None  # keep the job queued
    with TestClient(app) as c:
        job = _post(c, "tile", {"amplicons": 4}).json()["id"]
        assert c.delete(f"/api/jobs/{job}").status_code == 200
        assert c.get(f"/api/jobs/{job}").json()["status"] == "failed"
        assert c.delete(f"/api/jobs/{'a' * 32}").status_code == 404
        assert c.delete("/api/jobs/not-a-job").status_code == 404


def test_a_genes_job_reads_the_amplicon_table(tmp_path: Path) -> None:
    from vhirprimerlab.web.result import build_result

    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "amplicons.tsv").write_text("amplicon\tstart\n1\t1\n")
    res = build_result(tmp_path, ["-i", "input.fasta", "--tile-region", "F:1-900"], 0)
    assert res["mode"] == "tile"
    assert len(res["rows"]) == 1


def test_build_result_carries_the_edge_advice(tmp_path: Path) -> None:
    import json

    from vhirprimerlab.web.result import OUTPUT_FILES, build_result

    advice = [{"span": "F", "end": "start", "limit": 159, "extra_loss": 59,
               "edge_before": 43.01075268817204, "edge_after": 97.84946236559139}]
    (tmp_path / "out").mkdir()
    (tmp_path / "out" / "amplicons.tsv").write_text("amplicon\tstart\n1\t1\n")
    (tmp_path / "out" / "advice.json").write_text(json.dumps(advice))
    res = build_result(tmp_path, ["-i", "input.fasta", "--tile-region", "F:1-900"], 0)
    assert res["advice"] == advice
    assert "advice.json" not in OUTPUT_FILES  # machine-readable, not a download


def test_build_result_without_advice_has_an_empty_list(tmp_path: Path) -> None:
    from vhirprimerlab.web.result import build_result

    (tmp_path / "out").mkdir()
    assert build_result(tmp_path, ["check"], 0)["advice"] == []
    assert build_result(tmp_path, ["-n", "4"], 1)["advice"] == []
