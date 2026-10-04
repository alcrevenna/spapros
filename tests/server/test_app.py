import json

import pytest
import scanpy as sc
from fastapi.testclient import TestClient

from spapros.server.app import create_app
from tests.server.conftest import FAKES
from tests.server.conftest import wait_for


@pytest.fixture(scope="module")
def h5ad(tmp_path_factory):
    adata = sc.read_h5ad("data/small_data_raw_counts.h5ad")
    path = tmp_path_factory.mktemp("data") / "small.h5ad"
    adata[:200].copy().write_h5ad(path)
    return path


@pytest.fixture()
def client(tmp_path):
    app = create_app(data_dir=tmp_path, pipeline=f"{FAKES}:succeed")
    with TestClient(app) as c:
        yield c


def post_job(client, h5ad, marker_list=None, filename="small.h5ad", **options):
    options = {"celltype_key": "celltype", "n": 10, **options}
    files = {"file": (filename, h5ad.read_bytes(), "application/octet-stream")}
    if marker_list is not None:
        files["marker_list"] = ("markers.csv", marker_list.read_bytes(), "text/csv")
    return client.post("/jobs", files=files, data={"options": json.dumps(options)})


def test_upload_runs_job_and_serves_results(client, h5ad):
    r = post_job(client, h5ad, n=12)
    assert r.status_code == 201, r.text
    job = r.json()
    assert job["input_filename"] == "small.h5ad"
    assert job["options"]["n"] == 12
    assert job["options"]["genes_key"] == "highly_variable"

    store = client.app.state.store
    wait_for(store, job["id"], "succeeded")
    assert (store.path(job["id"]) / "input.h5ad").stat().st_size == h5ad.stat().st_size

    results = client.get(f"/jobs/{job['id']}/results").json()
    assert results["genes"] == ["A", "B"]
    assert results["job"]["status"] == "succeeded"

    f = client.get(f"/jobs/{job['id']}/files/probeset.csv")
    assert f.status_code == 200 and f.text.startswith("gene,selection")
    assert client.get(f"/jobs/{job['id']}/files/../job.json").status_code == 404
    assert "stage: selecting" in client.get(f"/jobs/{job['id']}/log").text
    assert [j["id"] for j in client.get("/jobs").json()] == [job["id"]]


def test_marker_list_is_stored(client, h5ad):
    from pathlib import Path

    job = post_job(client, h5ad, marker_list=Path("data/small_data_marker_list.csv")).json()
    assert (client.app.state.store.path(job["id"]) / "marker_list.csv").exists()


@pytest.mark.parametrize(
    "options, message",
    [
        ({"celltype_key": "nope"}, "cell type column 'nope' not found"),
        ({"genes_key": "missing"}, "gene subset column 'missing' not found"),
        ({"genes_key": "gene_ids"}, "must be boolean"),
        ({"n": 100000, "panel_capacity": 100000}, "larger than the"),
        ({"preselected_genes": ["NOT_A_GENE"]}, "NOT_A_GENE"),
        ({"critical_celltypes": ["celltype_1", "nope"]}, "critical cell types not found"),
    ],
)
def test_upload_rejects_options_that_do_not_fit_the_data(client, h5ad, options, message):
    r = post_job(client, h5ad, **options)
    assert r.status_code == 422
    assert message in r.json()["detail"]
    # Rejected uploads leave nothing behind.
    assert client.get("/jobs").json() == []


def test_upload_rejects_bad_options_and_files(client, h5ad, tmp_path):
    assert post_job(client, h5ad, n=0).status_code == 422
    assert post_job(client, h5ad, unknown_option=1).status_code == 422
    assert post_job(client, h5ad, filename="data.csv").status_code == 422
    junk = tmp_path / "junk.h5ad"
    junk.write_text("not hdf5")
    r = post_job(client, junk, filename="junk.h5ad")
    assert r.status_code == 422 and "could not read" in r.json()["detail"]
    assert client.get("/jobs").json() == []


def test_unknown_job_and_unfinished_results(tmp_path, h5ad):
    app = create_app(data_dir=tmp_path, pipeline=f"{FAKES}:succeed", start_queue=False)
    with TestClient(app) as client:
        assert client.get("/jobs/doesnotexist").status_code == 404
        job = post_job(client, h5ad).json()
        assert client.get(f"/jobs/{job['id']}/results").status_code == 409
        assert client.post(f"/jobs/{job['id']}/retry").status_code == 409
        assert client.delete(f"/jobs/{job['id']}").status_code == 409
        assert client.post(f"/jobs/{job['id']}/cancel").json()["status"] == "cancelled"
        assert client.post(f"/jobs/{job['id']}/cancel").status_code == 409
        assert client.delete(f"/jobs/{job['id']}").status_code == 204
        assert client.get(f"/jobs/{job['id']}").status_code == 404


def test_failed_job_can_be_retried(tmp_path, h5ad):
    app = create_app(data_dir=tmp_path, pipeline=f"{FAKES}:fail_once")
    with TestClient(app) as client:
        job = post_job(client, h5ad).json()
        store = client.app.state.store
        assert wait_for(store, job["id"], "failed")["error"] == "RuntimeError: interrupted"
        assert client.post(f"/jobs/{job['id']}/retry").status_code == 200
        wait_for(store, job["id"], "succeeded")
        assert client.get(f"/jobs/{job['id']}/results").json()["resumed_from"] == "half done"


def test_bad_budget_is_rejected(client, h5ad):
    r = post_job(client, h5ad, n=200, panel_capacity=150)
    assert r.status_code == 422 and "gene budget" in r.text


def test_upload_then_run_from_upload(client, h5ad):
    r = client.post("/uploads", files={"file": ("small.h5ad", h5ad.read_bytes())})
    assert r.status_code == 201, r.text
    upload = r.json()
    assert (upload["n_cells"], upload["n_genes"]) == (200, 2000)
    assert upload["looks_like_counts"] is True
    celltype = next(c for c in upload["obs_columns"] if c["name"] == "celltype")
    assert {v["name"] for v in celltype["values"]} >= {"celltype_1", "celltype_6"}
    assert sum(v["count"] for v in celltype["values"]) == 200
    assert next(c for c in upload["obs_columns"] if c["name"] == "size_factors")["values"] is None
    assert {"name": "highly_variable", "n_true": 2000} in upload["var_bool_columns"]
    assert [u["id"] for u in client.get("/uploads").json()] == [upload["id"]]

    options = {"celltype_key": "celltype", "n": 10, "critical_celltypes": ["celltype_1"]}
    r = client.post("/jobs", data={"upload_id": upload["id"], "options": json.dumps(options)})
    assert r.status_code == 201, r.text
    job = r.json()
    assert job["input_filename"] == "small.h5ad"
    assert job["options"]["critical_celltypes"] == ["celltype_1"]
    assert job["options"]["panel_capacity"] == 300
    store = client.app.state.store
    wait_for(store, job["id"], "succeeded")
    # The job keeps its own copy, so deleting the upload does not affect it.
    assert client.delete(f"/uploads/{upload['id']}").status_code == 204
    assert (store.path(job["id"]) / "input.h5ad").stat().st_size == h5ad.stat().st_size
    assert client.get(f"/uploads/{upload['id']}").status_code == 404


def test_job_needs_exactly_one_dataset(client, h5ad):
    options = json.dumps({"celltype_key": "celltype", "n": 10})
    assert client.post("/jobs", data={"options": options}).status_code == 422
    assert client.post("/jobs", data={"options": options, "upload_id": "nothere"}).status_code == 404
    r = client.post(
        "/jobs",
        files={"file": ("small.h5ad", h5ad.read_bytes())},
        data={"options": options, "upload_id": "x"},
    )
    assert r.status_code == 422


def test_bad_uploads_leave_nothing_behind(client, tmp_path):
    junk = tmp_path / "junk.h5ad"
    junk.write_text("not hdf5")
    r = client.post("/uploads", files={"file": ("junk.h5ad", junk.read_bytes())})
    assert r.status_code == 422 and "could not read" in r.json()["detail"]
    assert client.post("/uploads", files={"file": ("x.csv", b"a,b")}).status_code == 422
    assert client.get("/uploads").json() == []


def test_gui_and_report(client, h5ad):
    page = client.get("/")
    assert page.status_code == 200 and "static/app.js" in page.text
    assert client.get("/static/app.js").status_code == 200
    assert client.get("/static/style.css").status_code == 200

    job = post_job(client, h5ad).json()
    store = client.app.state.store
    assert client.get(f"/jobs/{job['id']}/report").status_code in (404, 200)
    wait_for(store, job["id"], "succeeded")
    assert client.get(f"/jobs/{job['id']}/report").status_code == 404  # the fake pipeline writes no report
    results = store.path(job["id"]) / "results"
    (results / "report.html").write_text("<html>report</html>")
    r = client.get(f"/jobs/{job['id']}/report")
    assert r.status_code == 200 and r.headers["content-type"].startswith("text/html") and "report" in r.text
    r = client.get(f"/jobs/{job['id']}/report?download=1")
    assert "attachment" in r.headers["content-disposition"]


def test_password(tmp_path):
    app = create_app(data_dir=tmp_path, pipeline=f"{FAKES}:succeed", start_queue=False, password="s3cret")
    with TestClient(app) as c:
        assert c.get("/health").status_code == 200
        assert c.get("/").status_code == 401
        assert c.get("/jobs").headers["www-authenticate"].startswith("Basic")
        assert c.get("/jobs", auth=("anyone", "wrong")).status_code == 401
        assert c.get("/jobs", auth=("anyone", "s3cret")).status_code == 200
