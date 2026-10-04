"""The work one server job does: select a probe set with spapros, evaluate it, and write the feasibility report.

This module runs inside the job's worker process (see :mod:`spapros.server.jobs`).
"""

import re
import shutil
import time
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional

from spapros.server.options import RunOptions

SPAPROS_SET_ID = "spapros"
INPUT_FILE = "input.h5ad"
MARKER_FILE = "marker_list.csv"
REFERENCE_METHODS = ["PCA", "DE", "HVG", "random"]
N_RANDOM_SEEDS = 3
CURVE_DIR = "panel_size_curve"


def _set_id(name: str) -> str:
    """Turn a reference set name like ``random (seed=0)`` into a file-safe id like ``random_seed0``."""
    return re.sub(r"[^A-Za-z0-9_-]+", "", name.replace(" ", "_").replace("=", ""))


def run_pipeline(job_dir: Path, options: Dict[str, Any], report_stage: Callable[[str], None]) -> Dict[str, Any]:
    """Run probe set selection, evaluation and the feasibility report for one job.

    Intermediate results go to ``job_dir/selection`` (the selector's ``save_dir``) and ``job_dir/evaluation`` (the
    evaluator's ``results_dir``). Both are reloaded by spapros when they exist, so rerunning an interrupted job resumes
    instead of starting over.

    Args:
        job_dir: The job's directory, holding ``input.h5ad`` and optionally ``marker_list.csv``.
        options: :class:`RunOptions` as a dict.
        report_stage: Called with the name of each stage as it starts.

    Returns:
        A JSON-serialisable summary: the selected genes, summary metrics per set, reference methods that could not run on
        this dataset, the feasibility verdict, and the result files.
    """
    import pandas as pd
    import scanpy as sc

    import spapros as sp

    started = time.time()
    opts = RunOptions(**options)
    job_dir = Path(job_dir)
    results_dir = job_dir / "results"
    results_dir.mkdir(exist_ok=True)
    marker_file: Optional[str] = str(job_dir / MARKER_FILE) if (job_dir / MARKER_FILE).exists() else None

    report_stage("loading")
    adata = sc.read_h5ad(job_dir / INPUT_FILE)
    if opts.normalize:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    report_stage("selecting")
    selector_kwargs: Dict[str, Any] = dict(
        celltype_key=opts.celltype_key,
        genes_key=opts.genes_key,
        n=opts.n,
        preselected_genes=opts.preselected_genes,
        prior_genes=opts.prior_genes,
        n_pca_genes=opts.n_pca_genes,
        n_min_markers=opts.n_min_markers,
        marker_list=marker_file,
        seed=opts.seed,
        save_dir=str(job_dir / "selection"),
        n_jobs=opts.n_jobs,
        verbosity=0,
    )
    if opts.forest_hparams:
        selector_kwargs["forest_hparams"] = opts.forest_hparams
    selector = sp.se.ProbesetSelector(adata, **selector_kwargs)
    selector.select_probeset()
    probeset = selector.probeset
    probeset.to_csv(results_dir / "probeset.csv")
    genes = probeset[probeset["selection"]].index.tolist()

    gene_sets = {SPAPROS_SET_ID: genes}
    skipped_reference_sets: Dict[str, str] = {}
    if opts.evaluate_reference_sets:
        report_stage("selecting_reference_sets")
        # Baselines are for comparison only, so a method that fails on this dataset is skipped, not fatal.
        for method in REFERENCE_METHODS:
            try:
                reference = sp.se.select_reference_probesets(
                    adata,
                    n=len(genes),
                    genes_key=opts.genes_key,
                    obs_key=opts.celltype_key,
                    methods=[method],
                    # Several random sets, so "better than random" is not judged on one lucky or unlucky draw.
                    seeds=[opts.seed + i for i in range(N_RANDOM_SEEDS if method == "random" else 1)],
                    verbosity=0,
                )
            except Exception as e:
                skipped_reference_sets[method] = f"{type(e).__name__}: {e}"
                print(f"reference set {method} skipped: {skipped_reference_sets[method]}", flush=True)
                continue
            for name, df in reference.items():
                gene_sets[_set_id(name)] = df[df["selection"]].index.tolist()
        pd.DataFrame({set_id: pd.Series(g) for set_id, g in gene_sets.items()}).to_csv(
            results_dir / "gene_sets.csv", index=False
        )

    report_stage("evaluating")
    evaluator = sp.ev.ProbesetEvaluator(
        adata,
        celltype_key=opts.celltype_key,
        results_dir=str(job_dir / "evaluation"),
        scheme=opts.evaluation_scheme,
        marker_list=marker_file,
        verbosity=0,
        n_jobs=opts.n_jobs,
    )
    evaluated = []
    for set_id, set_genes in gene_sets.items():
        try:
            evaluator.evaluate_probeset(set_genes, set_id=set_id, update_summary=False)
        except Exception as e:
            if set_id == SPAPROS_SET_ID:
                raise
            skipped_reference_sets[set_id] = f"evaluation failed: {type(e).__name__}: {e}"
            print(f"reference set {set_id} skipped: {skipped_reference_sets[set_id]}", flush=True)
            continue
        evaluated.append(set_id)
    evaluator.summary_statistics(set_ids=evaluated)

    summary = evaluator.summary_results
    summary.to_csv(results_dir / "evaluation_summary.csv")
    # Per cell type confusion matrices of the forest classifiers: the core input for a feasibility verdict.
    for set_id in evaluated:
        clf_file = job_dir / "evaluation" / "forest_clfs" / f"forest_clfs_{evaluator.ref_name}_{set_id}.csv"
        if clf_file.exists():
            shutil.copy(clf_file, results_dir / f"confusion_matrix_{set_id}.csv")
    adata.obs[opts.celltype_key].value_counts().rename("n_cells").rename_axis("celltype").to_csv(
        results_dir / "cell_counts.csv"
    )

    if opts.panel_size_curve:
        report_stage("panel_size_curve")
        evaluate_panel_size_curve(adata, probeset, opts, job_dir)

    report_stage("reporting")
    verdict = write_report(
        job_dir,
        opts,
        skipped_sets=skipped_reference_sets,
        run={"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars), "runtime_seconds": time.time() - started},
    )

    return {
        "genes": genes,
        "n_genes": len(genes),
        "summary": {
            set_id: {metric: (None if pd.isna(v) else float(v)) for metric, v in row.items()}
            for set_id, row in summary.iterrows()
        },
        "skipped_reference_sets": skipped_reference_sets,
        "verdict": {"verdict": verdict["verdict"], "reason": verdict["reason"]},
        "files": sorted(p.name for p in results_dir.iterdir() if p.is_file()),
    }


def ranked_genes(probeset) -> List[str]:
    """Genes of ``selector.probeset`` in selection order. The first ``n`` are the selected panel."""
    ranked = probeset[probeset["rank"].notna() | (probeset["pca_score"] > 0)]
    return ranked.sort_values("gene_nr").index.tolist()


def evaluate_panel_size_curve(adata, probeset, opts: RunOptions, job_dir: Path) -> None:
    """Classify cell types with the top ``k`` ranked genes for several ``k`` around the budget.

    Only the ``forest_clfs`` metric is computed, which is much cheaper than a full evaluation. The confusion matrices go
    to ``results/panel_size_curve/top_<k>.csv``; the one at the budget is the spapros panel's own.
    """
    import spapros as sp
    from spapros.server.feasibility import curve_sizes

    genes = ranked_genes(probeset)
    out_dir = job_dir / "results" / CURVE_DIR
    out_dir.mkdir(exist_ok=True)
    evaluator = sp.ev.ProbesetEvaluator(
        adata,
        celltype_key=opts.celltype_key,
        results_dir=str(job_dir / "evaluation_curve"),
        scheme="custom",
        metrics=["forest_clfs"],
        verbosity=0,
        n_jobs=opts.n_jobs,
    )
    spapros_cm = job_dir / "results" / f"confusion_matrix_{SPAPROS_SET_ID}.csv"
    for k in curve_sizes(opts.n, len(genes)):
        if k == opts.n and spapros_cm.exists():
            shutil.copy(spapros_cm, out_dir / f"top_{k}.csv")
            continue
        set_id = f"top_{k}"
        evaluator.evaluate_probeset(genes[:k], set_id=set_id, update_summary=False)
        clf_file = job_dir / "evaluation_curve" / "forest_clfs" / f"forest_clfs_{evaluator.ref_name}_{set_id}.csv"
        shutil.copy(clf_file, out_dir / f"top_{k}.csv")


def write_report(
    job_dir: Path, opts: RunOptions, skipped_sets: Optional[Dict[str, str]] = None, run: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """Judge feasibility from the job's result files and write ``results/verdict.json`` and ``results/report.html``.

    Everything is read back from ``job_dir/results``, so the report can be rebuilt without rerunning spapros.
    """
    import json

    import pandas as pd

    import spapros as sp
    from spapros.server.feasibility import assess
    from spapros.server.report import render_report

    results_dir = Path(job_dir) / "results"
    confusion = {
        p.stem[len("confusion_matrix_") :]: pd.read_csv(p, index_col=0)
        for p in sorted(results_dir.glob("confusion_matrix_*.csv"))
    }
    summary_file = results_dir / "evaluation_summary.csv"
    summary = pd.read_csv(summary_file, index_col=0) if summary_file.exists() else None
    counts = pd.read_csv(results_dir / "cell_counts.csv", index_col=0)["n_cells"]
    curve = {
        int(p.stem[len("top_") :]): pd.read_csv(p, index_col=0) for p in (results_dir / CURVE_DIR).glob("top_*.csv")
    }
    probeset = pd.read_csv(results_dir / "probeset.csv", index_col=0)
    n_required = int((probeset["pre_selected"].astype(bool) | probeset["required_marker"].astype(bool)).sum())

    verdict = assess(
        cell_counts={str(k): int(v) for k, v in counts.items() if v > 0},
        confusion=confusion,
        summary=summary,
        budget=opts.n,
        panel_capacity=opts.panel_capacity,
        reserved_slots=opts.reserved_slots,
        critical_celltypes=opts.critical_celltypes,
        curve=curve,
        n_required_genes=n_required,
        skipped_sets=skipped_sets,
        settings=opts.feasibility,
    )
    (results_dir / "verdict.json").write_text(json.dumps(verdict, indent=2))

    job_meta = {}
    if (Path(job_dir) / "job.json").exists():
        job_meta = json.loads((Path(job_dir) / "job.json").read_text())
    run_info = {
        "dataset": job_meta.get("input_filename", "dataset"),
        "created": job_meta.get("created_at"),
        "celltype_key": opts.celltype_key,
        "spapros_version": sp.__version__,
        "options": opts.model_dump(),
        **(run or {}),
    }
    files = sorted(p.name for p in results_dir.iterdir() if p.is_file() and p.name != "report.html")
    html = render_report(verdict, run_info, confusion[SPAPROS_SET_ID], probeset=probeset, files=files)
    (results_dir / "report.html").write_text(html)
    return verdict
