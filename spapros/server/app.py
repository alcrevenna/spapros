"""HTTP API and browser GUI of the spapros server (FastAPI).

The GUI is served at ``/``. Endpoints:

- ``POST /uploads``: upload an ``.h5ad`` (``file``) and get back what the run form needs (cell type columns with their
  counts, gene subsets). ``GET /uploads``, ``GET /uploads/{id}``, ``DELETE /uploads/{id}``.
- ``POST /jobs``: the dataset as ``upload_id`` (from ``POST /uploads``) or as a new ``file``, optional marker list CSV
  (``marker_list``) and :class:`RunOptions` as a JSON string (``options``). Returns the queued job.
- ``GET /jobs`` and ``GET /jobs/{id}``: job status, including the current ``stage`` while it runs.
- ``GET /jobs/{id}/results``: selected genes and summary metrics per probe set, once the job succeeded.
- ``GET /jobs/{id}/files/{name}``: download one result file (probeset.csv, evaluation_summary.csv, ...).
- ``GET /jobs/{id}/report``: the feasibility report (HTML); ``?download=1`` serves it as an attachment.
- ``GET /jobs/{id}/log``: the end of the job's log.
- ``POST /jobs/{id}/cancel``, ``POST /jobs/{id}/retry``, ``DELETE /jobs/{id}``.
"""

import base64
import os
import secrets
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
from fastapi import Request
from fastapi import UploadFile
from fastapi.responses import FileResponse
from fastapi.responses import PlainTextResponse
from fastapi.responses import Response
from fastapi.staticfiles import StaticFiles
from pydantic import ValidationError

from spapros.server.jobs import DEFAULT_PIPELINE
from spapros.server.jobs import FINISHED
from spapros.server.jobs import JobQueue
from spapros.server.jobs import JobStatus
from spapros.server.jobs import JobStore
from spapros.server.options import RunOptions
from spapros.server.pipeline import INPUT_FILE
from spapros.server.pipeline import MARKER_FILE
from spapros.server.uploads import UploadStore

STATIC_DIR = Path(__file__).parent / "static"

DATA_DIR_ENV = "SPAPROS_SERVER_DATA_DIR"
WORKERS_ENV = "SPAPROS_SERVER_WORKERS"
PASSWORD_ENV = "SPAPROS_SERVER_PASSWORD"


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
        celltypes = set(adata.obs[options.celltype_key].astype(str))
        unknown = [ct for ct in options.critical_celltypes if ct not in celltypes]
        if unknown:
            raise ValueError(f"critical cell types not found in {options.celltype_key!r}: {unknown}")
    finally:
        adata.file.close()


def _save_upload(upload: UploadFile, dest: Path) -> None:
    with open(dest, "wb") as f:
        shutil.copyfileobj(upload.file, f, length=16 * 1024 * 1024)


def _check_password(request: Request, password: str) -> bool:
    """HTTP basic auth: any user name, the configured password."""
    header = request.headers.get("authorization", "")
    if not header.lower().startswith("basic "):
        return False
    try:
        _, _, given = base64.b64decode(header[6:]).decode().partition(":")
    except Exception:
        return False
    return secrets.compare_digest(given.encode(), password.encode())


def create_app(
    data_dir: Optional[os.PathLike] = None,
    pipeline: str = DEFAULT_PIPELINE,
    workers: Optional[int] = None,
    start_queue: bool = True,
    password: Optional[str] = None,
) -> FastAPI:
    """Build the app.

    Args:
        data_dir: Where uploads and job results are stored. Defaults to ``$SPAPROS_SERVER_DATA_DIR`` or
            ``./spapros_server_data``.
        pipeline: ``"module:function"`` each job runs; tests swap in a fast fake.
        workers: Jobs run at the same time. Defaults to ``$SPAPROS_SERVER_WORKERS`` or 1.
        start_queue: Start the background queue with the app. Disable to only accept and inspect jobs.
        password: Require this password (HTTP basic auth, any user name) for every request. Defaults to
            ``$SPAPROS_SERVER_PASSWORD``; unset means no authentication.
    """
    data_dir = Path(data_dir or os.environ.get(DATA_DIR_ENV, "spapros_server_data"))
    store = JobStore(data_dir)
    uploads = UploadStore(data_dir)
    password = password if password is not None else os.environ.get(PASSWORD_ENV) or None
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
    app.state.uploads = uploads
    app.state.queue = job_queue

    if password:

        @app.middleware("http")
        async def require_password(request: Request, call_next):
            if request.url.path != "/health" and not _check_password(request, password):
                return Response(status_code=401, headers={"WWW-Authenticate": 'Basic realm="spapros"'})
            return await call_next(request)

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")

    @app.get("/", include_in_schema=False)
    def gui() -> FileResponse:
        return FileResponse(STATIC_DIR / "index.html")

    def get_job(job_id: str) -> Dict[str, Any]:
        job = store.get(job_id)
        if job is None:
            raise HTTPException(404, f"job {job_id} not found")
        return job

    @app.get("/health")
    def health() -> Dict[str, str]:
        return {"status": "ok"}

    @app.post("/uploads", status_code=201)
    def create_upload(
        file: Annotated[UploadFile, File(description="AnnData .h5ad with raw counts and a cell type column.")],
    ) -> Dict[str, Any]:
        filename = file.filename or ""
        if not filename.endswith(".h5ad"):
            raise HTTPException(422, "file must be an .h5ad file")
        try:
            return uploads.create(lambda dest: _save_upload(file, dest), filename)
        except ValueError as e:
            raise HTTPException(422, str(e)) from e

    @app.get("/uploads")
    def list_uploads() -> List[Dict[str, Any]]:
        return uploads.list()

    def get_upload(upload_id: str) -> Dict[str, Any]:
        meta = uploads.get(upload_id)
        if meta is None:
            raise HTTPException(404, f"upload {upload_id} not found")
        return meta

    @app.get("/uploads/{upload_id}")
    def read_upload(upload_id: str) -> Dict[str, Any]:
        return get_upload(upload_id)

    @app.delete("/uploads/{upload_id}", status_code=204)
    def delete_upload(upload_id: str) -> None:
        get_upload(upload_id)
        uploads.delete(upload_id)

    @app.post("/jobs", status_code=201)
    def create_job(
        options: Annotated[str, Form(description="RunOptions as JSON.")],
        file: Annotated[
            Optional[UploadFile], File(description="AnnData .h5ad with raw counts and a cell type column.")
        ] = None,
        upload_id: Annotated[Optional[str], Form(description="A dataset uploaded with POST /uploads.")] = None,
        marker_list: Annotated[Optional[UploadFile], File(description="Optional marker list CSV.")] = None,
    ) -> Dict[str, Any]:
        try:
            run_options = RunOptions.model_validate_json(options)
        except ValidationError as e:
            raise HTTPException(422, e.errors(include_url=False, include_context=False)) from e
        if (file is None) == (upload_id is None):
            raise HTTPException(422, "send the dataset either as file or as upload_id")
        if file is not None:
            filename = file.filename or ""
            if not filename.endswith(".h5ad"):
                raise HTTPException(422, "file must be an .h5ad file")
        else:
            filename = get_upload(upload_id)["filename"]

        job = store.create(run_options.model_dump(), input_filename=filename)
        job_dir = store.path(job["id"])
        try:
            if file is not None:
                _save_upload(file, job_dir / INPUT_FILE)
            else:
                uploads.link_into(upload_id, job_dir / INPUT_FILE)
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

    @app.get("/jobs/{job_id}/report", include_in_schema=False)
    def read_report(job_id: str, download: bool = False) -> FileResponse:
        get_job(job_id)
        report = store.path(job_id) / "results" / "report.html"
        if not report.exists():
            raise HTTPException(404, "this job has no report (yet)")
        if download:
            return FileResponse(report, filename=f"feasibility_report_{job_id}.html")
        return FileResponse(report, media_type="text/html")

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
