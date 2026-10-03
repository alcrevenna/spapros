"""The work one server job does: select a probe set with spapros, then evaluate it.

This module runs inside the job's worker process (see :mod:`spapros.server.jobs`).
"""

import re
import shutil
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Dict
from typing import Optional

from spapros.server.options import RunOptions

SPAPROS_SET_ID = "spapros"
INPUT_FILE = "input.h5ad"
MARKER_FILE = "marker_list.csv"
REFERENCE_METHODS = ["PCA", "DE", "HVG", "random"]


def _set_id(name: str) -> str:
    """Turn a reference set name like ``random (seed=0)`` into a file-safe id like ``random_seed0``."""
    return re.sub(r"[^A-Za-z0-9_-]+", "", name.replace(" ", "_").replace("=", ""))


def run_pipeline(job_dir: Path, options: Dict[str, Any], report_stage: Callable[[str], None]) -> Dict[str, Any]:
    """Run probe set selection and evaluation for one job.

    Intermediate results go to ``job_dir/selection`` (the selector's ``save_dir``) and ``job_dir/evaluation`` (the
    evaluator's ``results_dir``). Both are reloaded by spapros when they exist, so rerunning an interrupted job resumes
    instead of starting over.

    Args:
        job_dir: The job's directory, holding ``input.h5ad`` and optionally ``marker_list.csv``.
        options: :class:`RunOptions` as a dict.
        report_stage: Called with the name of each stage as it starts.

    Returns:
        A JSON-serialisable summary: the selected genes, summary metrics per set, reference methods that could not run on
        this dataset, and the result files.
    """
    import pandas as pd
    import scanpy as sc

    import spapros as sp

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
                    seeds=[opts.seed],
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

    return {
        "genes": genes,
        "n_genes": len(genes),
        "summary": {
            set_id: {metric: (None if pd.isna(v) else float(v)) for metric, v in row.items()}
            for set_id, row in summary.iterrows()
        },
        "skipped_reference_sets": skipped_reference_sets,
        "files": sorted(p.name for p in results_dir.iterdir() if p.is_file()),
    }
