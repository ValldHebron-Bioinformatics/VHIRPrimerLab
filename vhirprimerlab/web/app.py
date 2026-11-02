"""The VHIRPrimerLab web app: a page, and a small job API around the command-line tool."""

from __future__ import annotations

import contextlib
import io
import json
import os
import re
import shutil
import tempfile
import zipfile
from pathlib import Path

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi import Request
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, PlainTextResponse, Response
from fastapi.staticfiles import StaticFiles

from vhirprimerlab import __version__
from vhirprimerlab.web.jobs import OUTPUT_FILES, JobManager, QueueFull, Target, run_job
from vhirprimerlab.web.params import Limits
from vhirprimerlab.web.prepare import prepare_job

STATIC = Path(__file__).parent / "static"
JOB_ID = re.compile(r"^[0-9a-f]{32}$")


def _aligners() -> list[str]:
    return [name for name in ("mafft", "muscle") if shutil.which(name)]


def _bad(message: str) -> HTTPException:
    return HTTPException(status_code=400, detail=message)


def create_app(
    root: Path | None = None, limits: Limits | None = None, target: Target = run_job
) -> FastAPI:
    limits = limits or Limits.from_env()
    root = Path(root or os.environ.get("VHIRPRIMERLAB_JOBS_DIR")
                or Path(tempfile.gettempdir()) / "vhirprimerlab-jobs")
    manager = JobManager(root, limits, target)

    @contextlib.asynccontextmanager
    async def lifespan(_app: FastAPI):
        manager.cleanup()
        manager.start()
        yield
        manager.stop()

    app = FastAPI(title="VHIRPrimerLab", version=__version__, lifespan=lifespan)
    app.state.manager = manager
    app.mount("/static", StaticFiles(directory=STATIC), name="static")
    static_files = {p.name for p in STATIC.iterdir() if p.is_file() and p.name != "index.html"}
    max_body = int(limits.max_upload_mb * 1024 * 1024) + 256 * 1024  # file plus form fields

    @app.middleware("http")
    async def refuse_large_bodies(request: Request, call_next):
        size = request.headers.get("content-length", "")
        if request.method == "POST" and size.isdigit() and int(size) > max_body:
            return JSONResponse(status_code=413, content={"detail": (
                f"The upload is larger than {limits.max_upload_mb:g} MB; "
                "use the command line for larger inputs.")})
        return await call_next(request)

    page = (STATIC / "index.html").read_text().replace(
        '<meta charset="utf-8">',
        '<meta charset="utf-8">\n<meta name="vhirprimerlab-server" content="1">', 1)

    @app.get("/", include_in_schema=False)
    def index() -> HTMLResponse:
        """The page, marked so it uses this server instead of running in the browser."""
        return HTMLResponse(page)

    @app.get("/{name}", include_in_schema=False)
    def page_file(name: str) -> FileResponse:
        """The page's own files at the root too, so its relative links work here and on Pages."""
        if name not in static_files:
            raise HTTPException(status_code=404, detail="not found")
        return FileResponse(STATIC / name)

    @app.get("/api/example", response_class=PlainTextResponse)
    def example() -> str:
        return (STATIC / "example.fasta").read_text()

    @app.get("/api/limits")
    def get_limits() -> dict[str, object]:
        return {
            "max_upload_mb": limits.max_upload_mb,
            "max_seqs": limits.max_seqs,
            "max_seq_len": limits.max_seq_len,
            "max_amplicons": limits.max_amplicons,
            "retention_hours": limits.retention_s / 3600,
            "timeout_minutes": limits.timeout_s / 60,
            "aligners": _aligners(),
        }

    @app.post("/api/jobs")
    async def create_job(
        mode: str = Form("tile"),
        params: str = Form("{}"),
        fasta_text: str = Form(""),
        primers_text: str = Form(""),
        fasta: UploadFile | None = File(None),
    ) -> dict[str, str]:
        try:
            settings = json.loads(params)
        except ValueError:
            raise _bad("settings are not valid JSON") from None
        max_bytes = int(limits.max_upload_mb * 1024 * 1024)
        content = await fasta.read(max_bytes + 1) if fasta is not None else fasta_text.encode()
        try:
            argv, files = prepare_job(mode, settings, content, primers_text, limits, _aligners())
        except ValueError as exc:
            raise _bad(str(exc)) from None
        try:
            return {"id": manager.submit(files, argv)}
        except QueueFull:
            raise HTTPException(status_code=503, detail=(
                "The server is busy: too many jobs are waiting. Try again in a few minutes, "
                "or use the command line.")) from None

    def _job(job_id: str) -> None:
        if not JOB_ID.match(job_id):
            raise HTTPException(status_code=404, detail="no such job")

    @app.get("/api/jobs/{job_id}")
    def job_status(job_id: str) -> dict:
        _job(job_id)
        try:
            return manager.status(job_id)
        except KeyError:
            raise HTTPException(status_code=404, detail="no such job (results are deleted "
                                f"after {limits.retention_s / 3600:g} hours)") from None

    @app.delete("/api/jobs/{job_id}")
    def cancel_job(job_id: str) -> dict[str, bool]:
        _job(job_id)
        try:
            return {"cancelled": manager.cancel(job_id)}
        except KeyError:
            raise HTTPException(status_code=404, detail="no such job") from None

    @app.get("/api/jobs/{job_id}/files/{name}")
    def job_file(job_id: str, name: str) -> FileResponse:
        _job(job_id)
        path = manager.output(job_id, name)
        if path is None:
            raise HTTPException(status_code=404, detail="no such file")
        return FileResponse(path, media_type="text/plain", filename=name)

    @app.get("/api/jobs/{job_id}/zip")
    def job_zip(job_id: str) -> Response:
        _job(job_id)
        paths = [p for name in OUTPUT_FILES if (p := manager.output(job_id, name))]
        if not paths:
            raise HTTPException(status_code=404, detail="no such job")
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            for p in paths:
                zf.write(p, p.name)
        return Response(buf.getvalue(), media_type="application/zip", headers={
            "Content-Disposition": f'attachment; filename="vhirprimerlab-{job_id[:8]}.zip"'})

    return app


app = create_app()
