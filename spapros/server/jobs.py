"""Job store and background queue for the spapros server.

Every job lives in its own directory under ``<data_dir>/jobs/<job_id>/``::

    job.json        metadata and status, written only by the web process
    input.h5ad      the uploaded dataset
    marker_list.csv optional marker list
    progress.json   current stage, written by the worker process
    selection/      ProbesetSelector save_dir (lets a rerun resume)
    evaluation/     ProbesetEvaluator results_dir (lets a rerun resume)
    results/        probeset.csv, evaluation_summary.csv, confusion matrices, verdict.json, report.html, ...
    result.json     summary returned by the pipeline once the job succeeded
    error.txt       traceback if the job failed
    log.txt         stdout/stderr of the worker process

Each job runs in its own spawned process, so a crash or an out-of-memory kill only takes down that job, and a job can
be cancelled by terminating its process. Jobs that were queued or running when the server stopped are requeued on the
next start and resume from their ``selection/`` and ``evaluation/`` checkpoints.
"""

import importlib
import json
import multiprocessing
import os
import queue
import shutil
import sys
import threading
import traceback
import uuid
from datetime import datetime
from datetime import timezone
from enum import Enum
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional

DEFAULT_PIPELINE = "spapros.server.pipeline:run_pipeline"


class JobStatus(str, Enum):
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"


FINISHED = {JobStatus.SUCCEEDED, JobStatus.FAILED, JobStatus.CANCELLED}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _read_json(path: Path) -> Optional[Any]:
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def _write_json(path: Path, data: Any) -> None:
    """Write atomically so readers never see a half-written file."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2))
    os.replace(tmp, path)


class JobStore:
    """Job metadata on disk, safe to use from several threads of one process."""

    def __init__(self, data_dir: os.PathLike):
        self.root = Path(data_dir) / "jobs"
        self.root.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def path(self, job_id: str) -> Path:
        return self.root / job_id

    def create(self, options: Dict[str, Any], input_filename: str) -> Dict[str, Any]:
        """Create the job directory and metadata. The caller writes the input files before submitting the job."""
        job_id = uuid.uuid4().hex[:12]
        self.path(job_id).mkdir()
        job = {
            "id": job_id,
            "status": JobStatus.QUEUED.value,
            "created_at": _now(),
            "started_at": None,
            "finished_at": None,
            "input_filename": input_filename,
            "options": options,
            "attempts": 0,
            "error": None,
        }
        _write_json(self.path(job_id) / "job.json", job)
        return job

    def get(self, job_id: str) -> Optional[Dict[str, Any]]:
        if not job_id.isalnum():
            return None
        job = _read_json(self.path(job_id) / "job.json")
        if job is not None:
            progress = _read_json(self.path(job_id) / "progress.json") or {}
            job["stage"] = progress.get("stage")
            result = _read_json(self.path(job_id) / "result.json") or {}
            job["verdict"] = result.get("verdict")
        return job

    def list(self) -> List[Dict[str, Any]]:
        jobs = [self.get(p.name) for p in self.root.iterdir() if p.is_dir()]
        return sorted((j for j in jobs if j is not None), key=lambda j: j["created_at"])

    def update(self, job_id: str, **fields: Any) -> Dict[str, Any]:
        with self._lock:
            path = self.path(job_id) / "job.json"
            job = _read_json(path)
            if job is None:
                raise KeyError(job_id)
            job.update({k: (v.value if isinstance(v, JobStatus) else v) for k, v in fields.items()})
            _write_json(path, job)
            return job

    def result(self, job_id: str) -> Optional[Dict[str, Any]]:
        return _read_json(self.path(job_id) / "result.json")

    def delete(self, job_id: str) -> None:
        shutil.rmtree(self.path(job_id))


def _load_callable(dotted: str) -> Callable:
    module, _, name = dotted.partition(":")
    return getattr(importlib.import_module(module), name)


def _run_in_worker(job_dir: str, pipeline: str, options: Dict[str, Any]) -> None:
    """Entry point of the worker process for one job."""
    job_path = Path(job_dir)
    log = open(job_path / "log.txt", "a", buffering=1)
    # Redirect at the file descriptor level so output from C extensions and joblib children lands in the log too.
    os.dup2(log.fileno(), 1)
    os.dup2(log.fileno(), 2)
    sys.stdout = sys.stderr = log

    def report_stage(stage: str) -> None:
        print(f"[{_now()}] stage: {stage}", flush=True)
        _write_json(job_path / "progress.json", {"stage": stage, "updated_at": _now()})

    try:
        result = _load_callable(pipeline)(job_path, options, report_stage)
        _write_json(job_path / "result.json", result)
    except Exception:
        (job_path / "error.txt").write_text(traceback.format_exc())
        traceback.print_exc()
        log.flush()
        os._exit(1)
    log.flush()
    # Exit right away: a normal interpreter exit waits for idle joblib/loky workers, which linger for minutes.
    os._exit(0)


class JobQueue:
    """Runs queued jobs in background threads, one spawned process per job.

    Args:
        store: Where jobs are kept.
        pipeline: ``"module:function"`` run for each job, see :func:`spapros.server.pipeline.run_pipeline`.
        workers: How many jobs run at the same time. spapros already uses all CPUs per job (``n_jobs=-1``), so 1 is a
            sensible default.
    """

    def __init__(self, store: JobStore, pipeline: str = DEFAULT_PIPELINE, workers: int = 1, poll_interval: float = 0.5):
        self.store = store
        self.pipeline = pipeline
        self.workers = workers
        self.poll_interval = poll_interval
        self._queue: "queue.Queue[Optional[str]]" = queue.Queue()
        self._threads: List[threading.Thread] = []
        self._processes: Dict[str, multiprocessing.process.BaseProcess] = {}
        self._cancel_requested: set = set()
        self._stopping = threading.Event()
        self._lock = threading.Lock()
        self._ctx = multiprocessing.get_context("spawn")

    def start(self) -> None:
        """Requeue jobs left unfinished by a previous server run, then start the worker threads."""
        for job in self.store.list():
            if job["status"] in (JobStatus.QUEUED.value, JobStatus.RUNNING.value):
                self.store.update(job["id"], status=JobStatus.QUEUED)
                self._queue.put(job["id"])
        self._stopping.clear()
        for i in range(self.workers):
            t = threading.Thread(target=self._work, name=f"spapros-worker-{i}", daemon=True)
            t.start()
            self._threads.append(t)

    def stop(self, timeout: float = 10) -> None:
        """Stop the worker threads. Running jobs are terminated but stay ``running`` so the next start resumes them."""
        self._stopping.set()
        with self._lock:
            for proc in self._processes.values():
                proc.terminate()
        for _ in self._threads:
            self._queue.put(None)
        for t in self._threads:
            t.join(timeout)
        self._threads = []

    def submit(self, job_id: str) -> None:
        self._queue.put(job_id)

    def cancel(self, job_id: str) -> Dict[str, Any]:
        """Cancel a queued or running job. Returns the job as it is after the request."""
        with self._lock:
            job = self.store.get(job_id)
            if job is None:
                raise KeyError(job_id)
            if job["status"] == JobStatus.QUEUED.value:
                return self.store.update(job_id, status=JobStatus.CANCELLED, finished_at=_now())
            if job["status"] == JobStatus.RUNNING.value:
                self._cancel_requested.add(job_id)
                proc = self._processes.get(job_id)
                if proc is not None:
                    proc.terminate()
        return job

    def retry(self, job_id: str) -> Dict[str, Any]:
        """Requeue a failed or cancelled job. It resumes from the checkpoints the previous attempt left behind."""
        job = self.store.get(job_id)
        if job is None:
            raise KeyError(job_id)
        if job["status"] not in (JobStatus.FAILED.value, JobStatus.CANCELLED.value):
            raise ValueError(f"only failed or cancelled jobs can be retried, job is {job['status']}")
        for name in ("error.txt", "result.json", "progress.json"):
            (self.store.path(job_id) / name).unlink(missing_ok=True)
        job = self.store.update(job_id, status=JobStatus.QUEUED, error=None, finished_at=None)
        self.submit(job_id)
        return job

    def _work(self) -> None:
        while not self._stopping.is_set():
            job_id = self._queue.get()
            if job_id is None:
                return
            job = self.store.get(job_id)
            # Skip jobs cancelled while queued, deleted, or submitted twice.
            if job is None or job["status"] != JobStatus.QUEUED.value:
                continue
            self._run(job)

    def _run(self, job: Dict[str, Any]) -> None:
        job_id = job["id"]
        job_dir = self.store.path(job_id)
        with self._lock:
            # Re-check under the lock: the job may have been cancelled since it was taken off the queue.
            current = self.store.get(job_id)
            if self._stopping.is_set() or current is None or current["status"] != JobStatus.QUEUED.value:
                return
            proc = self._ctx.Process(
                target=_run_in_worker, args=(str(job_dir), self.pipeline, job["options"]), name=f"spapros-job-{job_id}"
            )
            self.store.update(job_id, status=JobStatus.RUNNING, started_at=_now(), attempts=job["attempts"] + 1)
            proc.start()
            self._processes[job_id] = proc

        while proc.is_alive():
            proc.join(self.poll_interval)

        with self._lock:
            self._processes.pop(job_id, None)
            cancelled = job_id in self._cancel_requested
            self._cancel_requested.discard(job_id)
            if self._stopping.is_set() and not cancelled:
                return  # Server shutdown: leave the job "running" so it resumes on the next start.

            if cancelled:
                self.store.update(job_id, status=JobStatus.CANCELLED, finished_at=_now())
            elif proc.exitcode == 0 and (job_dir / "result.json").exists():
                self.store.update(job_id, status=JobStatus.SUCCEEDED, finished_at=_now())
            else:
                self.store.update(job_id, status=JobStatus.FAILED, finished_at=_now(), error=self._error(job_dir, proc))

    @staticmethod
    def _error(job_dir: Path, proc: multiprocessing.process.BaseProcess) -> str:
        error_file = job_dir / "error.txt"
        if error_file.exists():
            lines = [line for line in error_file.read_text().splitlines() if line.strip()]
            return lines[-1] if lines else "unknown error"
        if proc.exitcode is not None and proc.exitcode < 0:
            return f"worker process was killed by signal {-proc.exitcode} (possibly out of memory)"
        return f"worker process exited with code {proc.exitcode}"
