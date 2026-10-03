"""End to end: upload a real dataset and let the server run spapros selection and evaluation on it."""

import json

import pandas as pd
import pytest
import scanpy as sc
from fastapi.testclient import TestClient

from spapros.server.app import create_app
from tests.server.conftest import wait_for


@pytest.fixture(scope="module")
def h5ad(tmp_path_factory):
    adata = sc.read_h5ad("data/small_data_raw_counts.h5ad")
    adata = adata[adata.obs["celltype"].isin(["celltype_1", "celltype_3", "celltype_6", "celltype_7"])]
    path = tmp_path_factory.mktemp("data") / "small.h5ad"
    adata.write_h5ad(path)
    return path


def test_real_pipeline(tmp_path, h5ad):
    options = {
        "celltype_key": "celltype",
        "n": 10,
        "forest_hparams": {"n_trees": 5, "subsample": 200, "test_subsample": 300},
        "n_jobs": 2,
    }
    with TestClient(create_app(data_dir=tmp_path)) as client:
        r = client.post(
            "/jobs",
            files={"file": ("small.h5ad", h5ad.read_bytes())},
            data={"options": json.dumps(options)},
        )
        assert r.status_code == 201, r.text
        job_id = r.json()["id"]
        job = wait_for(client.app.state.store, job_id, {"succeeded", "failed"}, timeout=600)
        assert job["status"] == "succeeded", client.get(f"/jobs/{job_id}/log").text

        results = client.get(f"/jobs/{job_id}/results").json()
        assert results["n_genes"] == len(results["genes"]) == 10
        assert set(results["summary"]) == {"spapros", "PCA", "DE", "HVG", "random_seed0"}
        assert results["skipped_reference_sets"] == {}
        assert 0 <= results["summary"]["spapros"]["forest_clfs accuracy"] <= 1
        assert {"probeset.csv", "evaluation_summary.csv", "gene_sets.csv", "confusion_matrix_spapros.csv"} <= set(
            results["files"]
        )

        confusion = pd.read_csv(tmp_path / "jobs" / job_id / "results" / "confusion_matrix_spapros.csv", index_col=0)
        assert set(confusion.index) == {"celltype_1", "celltype_3", "celltype_6", "celltype_7"}
        # Checkpoints the selector writes, which a rerun resumes from.
        assert (tmp_path / "jobs" / job_id / "selection" / "probeset.csv").exists()
