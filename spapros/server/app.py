"""HTTP API of the spapros server (FastAPI).

Endpoints:

- ``POST /jobs``: upload an ``.h5ad`` (``file``), optional marker list CSV (``marker_list``) and :class:`RunOptions` as
  a JSON string (``options``). Returns the queued job.
- ``GET /jobs`` and ``GET /jobs/{id}``: job status, including the current ``stage`` while it runs.
- ``GET /jobs/{id}/results``: selected genes and summary metrics per probe set, once the job succeeded.
- ``GET /jobs/{id}/files/{name}``: download one result file (probeset.csv, evaluation_summary.csv, ...).
- ``GET /jobs/{id}/log``: the end of the job's log.
- ``POST /jobs/{id}/cancel``, ``POST /jobs/{id}/retry``, ``DELETE /jobs/{id}``.
"""

import os
import shutil
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Annotated
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from fastapi import FastAPI
from fastapi import File
from fastapi import Form
from fastapi import HTTPException
from fastapi import UploadFile
from fastapi.responses import FileResponse
from fastapi.responses import PlainTextResponse
from pydantic import ValidationError

from spapros.server.jobs import DEFAULT_PIPELINE
from spapros.server.jobs import FINISHED
from spapros.server.jobs import JobQueue
from spapros.server.jobs import JobStatus
from spapros.server.jobs import JobStore
from spapros.server.options import RunOptions
from spapros.server.pipeline import INPUT_FILE
from spapros.server.pipeline import MARKER_FILE

DATA_DIR_ENV = "SPAPROS_SERVER_DATA_DIR"
WORKERS_ENV = "SPAPROS_SERVER_WORKERS"


def validate_input(h5ad_path: Path, options: RunOptions) -> None:
    """Check that the uploaded dataset fits the options, so mistakes fail at upload and not hours into a run.

    Raises:
        ValueError: With a message meant for the user.
    """
    import anndata

    try:
        adata = anndata.read_h5ad(h5ad_path, backed="r")
    except Exception as e:
        raise ValueError(f"could not read the file as .h5ad: {e}") from e
    try:
        if options.celltype_key not in adata.obs.columns:
            raise ValueError(
                f"cell type column {options.celltype_key!r} not found in adata.obs (columns: {list(adata.obs.columns)})"
            )
        n_celltypes = adata.obs[options.celltype_key].nunique()
        if n_celltypes < 2:
            raise ValueError(f"cell type column {options.celltype_key!r} needs at least 2 cell types")
        n_candidates = adata.n_vars
        if options.genes_key is not None:
            if options.genes_key not in adata.var.columns:
                raise ValueError(
                    f"gene subset column {options.genes_key!r} not found in adata.var; compute it (e.g. "
                    "sc.pp.highly_variable_genes) or set genes_key to null to select from all genes"
                )
            if adata.var[options.genes_key].dtype != bool:
                raise ValueError(f"adata.var[{options.genes_key!r}] must be boolean")
            n_candidates = int(adata.var[options.genes_key].sum())
        if options.n > n_candidates:
            raise ValueError(f"n={options.n} is larger than the {n_candidates} candidate genes")
        missing = [g for g in options.preselected_genes + options.prior_genes if g not in adata.var_names]
        if missing:
            raise ValueError(f"genes not found in adata.var_names: {missing[:20]}")
    finally:
        adata.file.close()


def _save_upload(upload: UploadFile, dest: Path) -> None:
    with open(dest, "wb") as f:
        shutil.copyfileobj(upload.file, f, length=16 * 1024 * 1024)


def create_app(
    data_dir: Optional[os.PathLike] = None,
    pipeline: str = DEFAULT_PIPELINE,
    workers: Optional[int] = None,
    start_queue: bool = True,
) -> FastAPI:
    """Build the app.

    Args:
        data_dir: Where uploads and job results are stored. Defaults to ``$SPAPROS_SERVER_DATA_DIR`` or
            ``./spapros_server_data``.
        pipeline: ``"module:function"`` each job runs; tests swap in a fast fake.
        workers: Jobs run at the same time. Defaults to ``$SPAPROS_SERVER_WORKERS`` or 1.
        start_queue: Start the background queue with the app. Disable to only accept and inspect jobs.
    """
    data_dir = Path(data_dir or os.environ.get(DATA_DIR_ENV, "spapros_server_data"))
    store = JobStore(data_dir)
    job_queue = JobQueue(store, pipeline=pipeline, workers=workers or int(os.environ.get(WORKERS_ENV, "1")))

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        if start_queue:
            job_queue.start()
        yield
        if start_queue:
            job_queue.stop()

    app = FastAPI(title="spapros server", lifespan=lifespan)
    app.state.store = store
    app.state.queue = job_queue

    def get_job(job_id: str) -> Dict[str, Any]:
        job = store.get(job_id)
        if job is None:
            raise HTTPException(404, f"job {job_id} not found")
        return job

    @app.get("/health")
    def health() -> Dict[str, str]:
        return {"status": "ok"}

    @app.post("/jobs", status_code=201)
    def create_job(
        file: Annotated[UploadFile, File(description="AnnData .h5ad with raw counts and a cell type column.")],
        options: Annotated[str, Form(description="RunOptions as JSON.")],
        marker_list: Annotated[Optional[UploadFile], File(description="Optional marker list CSV.")] = None,
    ) -> Dict[str, Any]:
        try:
            run_options = RunOptions.model_validate_json(options)
        except ValidationError as e:
            raise HTTPException(422, e.errors(include_url=False, include_context=False)) from e
        filename = file.filename or ""
        if not filename.endswith(".h5ad"):
            raise HTTPException(422, "file must be an .h5ad file")

        job = store.create(run_options.model_dump(), input_filename=filename)
        job_dir = store.path(job["id"])
        try:
            _save_upload(file, job_dir / INPUT_FILE)
            if marker_list is not None and marker_list.filename:
                _save_upload(marker_list, job_dir / MARKER_FILE)
            validate_input(job_dir / INPUT_FILE, run_options)
        except ValueError as e:
            store.delete(job["id"])
            raise HTTPException(422, str(e)) from e
        except BaseException:
            store.delete(job["id"])
            raise
        job_queue.submit(job["id"])
        return get_job(job["id"])

    @app.get("/jobs")
    def list_jobs() -> List[Dict[str, Any]]:
        return store.list()

    @app.get("/jobs/{job_id}")
    def read_job(job_id: str) -> Dict[str, Any]:
        return get_job(job_id)

    @app.get("/jobs/{job_id}/results")
    def read_results(job_id: str) -> Dict[str, Any]:
        job = get_job(job_id)
        if job["status"] != JobStatus.SUCCEEDED.value:
            raise HTTPException(409, f"job is {job['status']}, results exist only for succeeded jobs")
        return {"job": job, **(store.result(job_id) or {})}

    @app.get("/jobs/{job_id}/files/{name}")
    def read_file(job_id: str, name: str) -> FileResponse:
        get_job(job_id)
        results_dir = store.path(job_id) / "results"
        files = {p.name: p for p in results_dir.iterdir() if p.is_file()} if results_dir.exists() else {}
        if name not in files:
            raise HTTPException(404, f"no result file {name}")
        return FileResponse(files[name], filename=name)

    @app.get("/jobs/{job_id}/log", response_class=PlainTextResponse)
    def read_log(job_id: str, lines: int = 200) -> str:
        get_job(job_id)
        log = store.path(job_id) / "log.txt"
        if not log.exists():
            return ""
        return "\n".join(log.read_text(errors="replace").splitlines()[-lines:])

    @app.post("/jobs/{job_id}/cancel")
    def cancel_job(job_id: str) -> Dict[str, Any]:
        job = get_job(job_id)
        if job["status"] in {s.value for s in FINISHED}:
            raise HTTPException(409, f"job is already {job['status']}")
        job_queue.cancel(job_id)
        return get_job(job_id)

    @app.post("/jobs/{job_id}/retry")
    def retry_job(job_id: str) -> Dict[str, Any]:
        get_job(job_id)
        try:
            job_queue.retry(job_id)
        except ValueError as e:
            raise HTTPException(409, str(e)) from e
        return get_job(job_id)

    @app.delete("/jobs/{job_id}", status_code=204)
    def delete_job(job_id: str) -> None:
        job = get_job(job_id)
        if job["status"] not in {s.value for s in FINISHED}:
            raise HTTPException(409, "cancel the job before deleting it")
        store.delete(job_id)

    return app
