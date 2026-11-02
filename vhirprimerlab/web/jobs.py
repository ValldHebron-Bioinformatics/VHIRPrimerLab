"""Background jobs for the web app: a queue, a few worker processes, and cleanup.

Each job runs `vhirprimerlab` (cli.main) in its own process inside its own directory, so
logs never mix and a job that runs too long can be stopped.
"""

from __future__ import annotations

import json
import multiprocessing
import os
import shutil
import signal
import sys
import threading
import time
import traceback
import uuid
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from vhirprimerlab.web.params import Limits
from vhirprimerlab.web.result import OUTPUT_FILES, build_result

Target = Callable[[Path, list[str]], None]


def run_job(job_dir: Path, argv: list[str]) -> None:
    """Child-process entry point: run vhirprimerlab in job_dir and write result.json."""
    from vhirprimerlab.cli import main  # imported in the child

    job_dir = Path(job_dir)
    os.chdir(job_dir)
    with (job_dir / "log.txt").open("w") as log:
        sys.stdout = sys.stderr = log
        try:
            code = main([*argv, "-o", "out", "--progress-file", "progress.json"])
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
        except Exception:  # noqa: BLE001 - reported to the user as a failed job
            traceback.print_exc()
            code = 1
        log.flush()
    result = build_result(job_dir, argv, code)
    (job_dir / "result.json").write_text(json.dumps(result))


class QueueFull(Exception):
    """Too many jobs are waiting."""


def _launch(target: Target, job_dir: Path, argv: list[str]) -> None:
    """Run the job in its own process group, so a timeout can stop the aligner too."""
    os.setpgrp()
    target(job_dir, argv)


def _kill_group(process: Any) -> None:
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        process.kill()


@dataclass
class _Job:
    id: str
    dir: Path
    argv: list[str]
    state: str = "queued"
    process: Any = None
    started: float | None = None
    finished: float | None = None
    error: str | None = None
    result: dict[str, Any] | None = field(default=None, repr=False)


class JobManager:
    def __init__(self, root: Path, limits: Limits, target: Target = run_job) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.limits = limits
        self.target = target
        self._jobs: dict[str, _Job] = {}
        self._queue: deque[str] = deque()
        self._lock = threading.RLock()
        self._ctx = multiprocessing.get_context("spawn")
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def submit(self, files: dict[str, bytes], argv: list[str]) -> str:
        with self._lock:
            if len(self._queue) >= self.limits.max_queue:
                raise QueueFull(len(self._queue))
        job_id = uuid.uuid4().hex
        job_dir = self.root / job_id
        job_dir.mkdir()
        for name, data in files.items():
            (job_dir / name).write_bytes(data)
        with self._lock:
            self._jobs[job_id] = _Job(job_id, job_dir, list(argv))
            self._queue.append(job_id)
        return job_id

    def _finish(self, job: _Job, error: str | None = None) -> None:
        job.finished = time.time()
        result_file = job.dir / "result.json"
        if error is None and result_file.is_file():
            job.result = json.loads(result_file.read_text())
            if job.result.get("exit_code") == 1:
                error = job.result.get("error") or "the job failed"
        elif error is None:
            error = "the job stopped unexpectedly"
        job.state, job.error = ("failed", error) if error else ("done", None)

    def tick(self) -> None:
        """Reap finished or overdue jobs and start queued ones."""
        with self._lock:
            now = time.time()
            running = [j for j in self._jobs.values() if j.state == "running"]
            for job in running:
                if not job.process.is_alive():
                    job.process.join()
                    self._finish(job)
                elif now - job.started > self.limits.timeout_s:
                    _kill_group(job.process)
                    job.process.join(5)
                    minutes = self.limits.timeout_s / 60
                    self._finish(job, f"took longer than {minutes:g} minutes and was stopped; "
                                      "use the command line for large inputs")
            free = self.limits.workers - sum(j.state == "running" for j in self._jobs.values())
            while free > 0 and self._queue:
                job = self._jobs[self._queue.popleft()]
                job.process = self._ctx.Process(target=_launch,
                                                args=(self.target, job.dir, job.argv), daemon=True)
                job.process.start()
                job.state, job.started = "running", time.time()
                free -= 1

    def status(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self._jobs[job_id]  # KeyError if unknown or expired
            stage, progress = None, 0.0
            progress_file = job.dir / "progress.json"
            if job.state == "running" and progress_file.is_file():
                try:
                    data = json.loads(progress_file.read_text())
                    stage, progress = data.get("stage"), float(data.get("progress", 0))
                except (ValueError, OSError):
                    pass
            if job.state in ("done", "failed"):
                progress = 1.0
            return {
                "id": job.id,
                "status": job.state,
                "stage": stage,
                "progress": progress,
                "queue_position": (list(self._queue).index(job_id) + 1
                                   if job.state == "queued" else None),
                "result": job.result,
                "error": job.error,
            }

    def cancel(self, job_id: str) -> bool:
        """Stop a queued or running job (and its aligner); False if it already finished."""
        with self._lock:
            job = self._jobs[job_id]  # KeyError if unknown
            if job.state == "queued":
                self._queue.remove(job_id)
            elif job.state == "running":
                _kill_group(job.process)
                job.process.join(5)
            else:
                return False
            self._finish(job, "Cancelled")
            return True

    def output(self, job_id: str, name: str) -> Path | None:
        """Path of a finished job's output file, or None."""
        with self._lock:
            job = self._jobs.get(job_id)
        if job is None or name not in OUTPUT_FILES:
            return None
        path = job.dir / "out" / name
        return path if path.is_file() else None

    def cleanup(self, now: float | None = None) -> int:
        """Delete job directories finished more than retention_s ago (and stray old ones)."""
        now = time.time() if now is None else now
        removed = 0
        with self._lock:
            for path in list(self.root.iterdir()):
                if not path.is_dir():
                    continue
                job = self._jobs.get(path.name)
                if job is not None:
                    if job.finished is None or now - job.finished <= self.limits.retention_s:
                        continue
                    del self._jobs[path.name]
                elif now - path.stat().st_mtime <= self.limits.retention_s:
                    continue
                shutil.rmtree(path, ignore_errors=True)
                removed += 1
        return removed

    def _loop(self) -> None:
        last_cleanup = 0.0
        while not self._stop.wait(0.5):
            self.tick()
            if time.time() - last_cleanup > 600:
                self.cleanup()
                last_cleanup = time.time()

    def start(self) -> None:
        if self._thread is None:
            self._stop.clear()
            self._thread = threading.Thread(target=self._loop, daemon=True)
            self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(2)
            self._thread = None
        with self._lock:
            for job in self._jobs.values():
                if job.state == "running" and job.process.is_alive():
                    _kill_group(job.process)
