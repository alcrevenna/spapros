"""Feasibility verdict for one spapros run.

spapros reports metrics, not a go/no-go. This module turns the metrics of a run into one of four verdicts, following
the project's feasibility rule (v1):

- ``FEASIBLE``: the panel resolves the cell types that matter and fits the gene budget with headroom.
- ``FEASIBLE WITH CAVEATS``: workable, but some cell types are marginal or unresolved, or the budget is tight.
- ``NOT FEASIBLE``: critical cell types cannot be resolved at this panel size, or required genes do not fit.
- ``INCONCLUSIVE``: the input cannot support a judgement (too few cells per type, too many types excluded).

Everything here works on plain pandas objects so the verdict can be recomputed from a job's saved result files.
"""

import math
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
import pandas as pd

from spapros.server.options import FeasibilitySettings

FEASIBLE = "FEASIBLE"
FEASIBLE_WITH_CAVEATS = "FEASIBLE WITH CAVEATS"
NOT_FEASIBLE = "NOT FEASIBLE"
INCONCLUSIVE = "INCONCLUSIVE"

GO = "GO"
CAVEAT = "CAVEAT"
FAIL = "FAIL"

RESOLVED = "resolved"
MARGINAL = "marginal"
UNRESOLVED = "unresolved"
NOT_ASSESSED = "not assessed"

SPAPROS_SET_ID = "spapros"
BASELINE_SET_IDS = ["PCA", "DE", "HVG"]
RANDOM_PREFIX = "random"


def curve_sizes(budget: int, n_ranked: int) -> List[int]:
    """Panel sizes evaluated for the panel size curve: 25/50/75/100% of the budget, plus 125/150% when available."""
    sizes = [math.ceil(f * budget) for f in (0.25, 0.5, 0.75, 1.0)]
    sizes += [math.ceil(f * budget) for f in (1.25, 1.5) if math.ceil(f * budget) <= n_ranked]
    return sorted({k for k in sizes if 1 <= k <= n_ranked})


def per_type_accuracy(confusion: pd.DataFrame) -> pd.Series:
    """Recall per cell type: the diagonal of a row-normalised confusion matrix (rows: true type, columns: predicted)."""
    cm = confusion.rename(index=str, columns=str)
    return pd.Series({ct: float(cm.loc[ct, ct]) if ct in cm.columns else 0.0 for ct in cm.index}, dtype=float)


def classify(acc: float, s: FeasibilitySettings) -> str:
    if acc >= s.acc_resolved:
        return RESOLVED
    if acc >= s.acc_marginal:
        return MARGINAL
    return UNRESOLVED


def _level(value: Optional[float], go: float, caveat: float) -> Optional[str]:
    if value is None or pd.isna(value):
        return None
    if value >= go:
        return GO
    if value >= caveat:
        return CAVEAT
    return FAIL


def _fmt(x: Optional[float]) -> str:
    return "n/a" if x is None or pd.isna(x) else f"{x:.2f}"


def _pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def _summary_value(summary: pd.DataFrame, set_id: str, prefix: str) -> Optional[float]:
    """First summary column starting with ``prefix`` for ``set_id``, or None."""
    if summary is None or set_id not in summary.index:
        return None
    for col in summary.columns:
        if str(col).startswith(prefix):
            v = summary.loc[set_id, col]
            return None if pd.isna(v) else float(v)
    return None


def _mean_acc(set_id: str, confusion: Dict[str, pd.DataFrame], summary: pd.DataFrame) -> Optional[float]:
    v = _summary_value(summary, set_id, "forest_clfs accuracy")
    if v is None and set_id in confusion:
        v = float(per_type_accuracy(confusion[set_id]).mean())
    return v


def _frac_resolved(acc: pd.Series, s: FeasibilitySettings) -> float:
    return float((acc >= s.acc_resolved).mean()) if len(acc) else 0.0


def assess(
    cell_counts: Dict[str, int],
    confusion: Dict[str, pd.DataFrame],
    summary: Optional[pd.DataFrame],
    budget: int,
    panel_capacity: int,
    reserved_slots: int = 0,
    critical_celltypes: Optional[List[str]] = None,
    curve: Optional[Dict[int, pd.DataFrame]] = None,
    n_required_genes: int = 0,
    skipped_sets: Optional[Dict[str, str]] = None,
    settings: Optional[FeasibilitySettings] = None,
) -> Dict[str, Any]:
    """Compute the feasibility verdict.

    Args:
        cell_counts: Cells per annotated cell type in the input data.
        confusion: Mean ``forest_clfs`` confusion matrix per evaluated gene set id. Must contain ``"spapros"``. Cell
            types missing from it were not assessed (spapros skips types with too few cells).
        summary: ``ProbesetEvaluator.summary_results``, indexed by set id.
        budget: Genes spapros was allowed to select (``n``).
        panel_capacity: Genes the platform's panel holds.
        reserved_slots: Panel slots kept for genes added outside spapros.
        critical_celltypes: Cell types that must be resolved.
        curve: Confusion matrix of the top ``k`` ranked genes, per ``k``.
        n_required_genes: Genes that must be on the panel (preselected genes).
        skipped_sets: Reference sets that could not be selected or evaluated, with the reason.
        settings: Thresholds; the defaults of :class:`FeasibilitySettings` when omitted.

    Returns:
        A JSON-serialisable dict, written to ``verdict.json`` by the pipeline.
    """
    s = settings or FeasibilitySettings()
    critical = [str(c) for c in (critical_celltypes or [])]
    cell_counts = {str(k): int(v) for k, v in cell_counts.items()}
    curve = curve or {}
    summary = summary if summary is not None else pd.DataFrame()

    acc = per_type_accuracy(confusion[SPAPROS_SET_ID])
    assessed = list(acc.index)
    random_ids = [i for i in confusion if i.startswith(RANDOM_PREFIX)]
    baseline_ids = [i for i in BASELINE_SET_IDS if i in confusion]
    baseline_acc = {i: per_type_accuracy(confusion[i]) for i in baseline_ids}

    # Per cell type table.
    cm = confusion[SPAPROS_SET_ID].rename(index=str, columns=str)
    celltypes = []
    for ct in sorted(set(cell_counts) | set(assessed)):
        row: Dict[str, Any] = {
            "celltype": ct,
            "n_cells": cell_counts.get(ct),
            "critical": ct in critical,
            "acc": None,
            "class": NOT_ASSESSED,
            "confused_with": None,
            "confusion_rate": None,
            "best_baseline": None,
            "best_baseline_acc": None,
        }
        if ct in acc.index:
            a = float(acc[ct])
            row["acc"] = a
            row["class"] = classify(a, s)
            if a < s.acc_resolved:
                off = cm.loc[ct].drop(ct, errors="ignore")
                if len(off) and off.max() > 0:
                    row["confused_with"] = str(off.idxmax())
                    row["confusion_rate"] = float(off.max())
            best = [(i, float(b[ct])) for i, b in baseline_acc.items() if ct in b.index]
            if best:
                row["best_baseline"], row["best_baseline_acc"] = max(best, key=lambda x: x[1])
        celltypes.append(row)
    order = {UNRESOLVED: 0, MARGINAL: 1, NOT_ASSESSED: 2, RESOLVED: 3}
    celltypes.sort(key=lambda r: (order[r["class"]], r["acc"] if r["acc"] is not None else 0, r["celltype"]))
    by_type = {r["celltype"]: r for r in celltypes}

    # Input quality (decides INCONCLUSIVE).
    excluded = [ct for ct in cell_counts if ct not in acc.index]
    n_cells_total = sum(cell_counts.values()) or 1
    excluded_cells = sum(cell_counts[ct] for ct in excluded) / n_cells_total
    excluded_types = len(excluded) / (len(cell_counts) or 1)
    unknown_critical = [ct for ct in critical if ct not in by_type]
    input_problems = []
    if len(assessed) < 2:
        input_problems.append(
            f"only {len(assessed)} cell type(s) have enough cells to be assessed; at least 2 are needed"
        )
    if excluded_types > s.max_excluded_types:
        input_problems.append(
            f"{len(excluded)} of {len(cell_counts)} cell types ({_pct(excluded_types)}) have too few cells to be "
            f"assessed ({', '.join(excluded)}); the limit is {_pct(s.max_excluded_types)}"
        )
    if excluded_cells > s.max_excluded_cells:
        input_problems.append(
            f"{_pct(excluded_cells)} of cells belong to cell types that cannot be assessed; the limit is "
            f"{_pct(s.max_excluded_cells)}"
        )
    critical_excluded = [ct for ct in critical if ct in excluded]
    if critical_excluded:
        input_problems.append(f"critical cell type(s) {', '.join(critical_excluded)} have too few cells to be assessed")
    if unknown_critical:
        input_problems.append(f"critical cell type(s) {', '.join(unknown_critical)} do not occur in the data")
    input_quality = {
        "ok": not input_problems,
        "problems": input_problems,
        "n_celltypes": len(cell_counts),
        "n_assessed": len(assessed),
        "not_assessed": excluded,
        "excluded_types_share": excluded_types,
        "excluded_cells_share": excluded_cells,
    }

    # Panel level resolution.
    frac = _frac_resolved(acc, s)
    mean_acc = _mean_acc(SPAPROS_SET_ID, confusion, summary)
    critical_classes = [by_type[ct]["class"] for ct in critical if ct in acc.index]
    critical_unresolved = [ct for ct in critical if ct in acc.index and by_type[ct]["class"] == UNRESOLVED]
    critical_marginal = [ct for ct in critical if ct in acc.index and by_type[ct]["class"] == MARGINAL]
    panel = {
        "frac_resolved": frac,
        "frac_resolved_level": _level(frac, s.frac_resolved_go, s.frac_resolved_caveat),
        "mean_acc": mean_acc,
        "mean_acc_level": _level(mean_acc, s.mean_acc_go, s.mean_acc_caveat),
        "smoothed_perct_acc": _summary_value(summary, SPAPROS_SET_ID, "forest_clfs perct acc"),
        "critical_level": (
            None if not critical_classes else FAIL if critical_unresolved else CAVEAT if critical_marginal else GO
        ),
        "n_resolved": int(sum(r["class"] == RESOLVED for r in celltypes)),
        "n_marginal": int(sum(r["class"] == MARGINAL for r in celltypes)),
        "n_unresolved": int(sum(r["class"] == UNRESOLVED for r in celltypes)),
    }
    levels = [panel["frac_resolved_level"], panel["mean_acc_level"]]
    panel_level = FAIL if FAIL in levels else CAVEAT if CAVEAT in levels else GO

    # Baselines.
    set_mean_acc = {i: _mean_acc(i, confusion, summary) for i in confusion}
    random_acc = float(np.mean([set_mean_acc[i] for i in random_ids])) if random_ids else None
    best_baseline = max(baseline_ids, key=lambda i: set_mean_acc[i] or 0) if baseline_ids else None
    beats_random = None if random_acc is None or mean_acc is None else mean_acc - random_acc >= s.beats_random_margin
    within_best = (
        None
        if best_baseline is None or mean_acc is None
        else mean_acc >= (set_mean_acc[best_baseline] or 0) - s.baseline_tolerance
    )
    rescues = [
        {"celltype": r["celltype"], "baseline": r["best_baseline"], "acc": r["best_baseline_acc"]}
        for r in celltypes
        if r["class"] == UNRESOLVED and (r["best_baseline_acc"] or 0) >= s.acc_resolved
    ]
    baselines = {
        "mean_acc": set_mean_acc,
        "random_mean_acc": random_acc,
        "random_sets": random_ids,
        "best_baseline": best_baseline,
        "beats_random": beats_random,
        "within_best_baseline": within_best,
        "rescues": rescues,
        "skipped": skipped_sets or {},
    }

    # Fit to the panel size.
    curve_points = []
    for k in sorted(curve):
        a = per_type_accuracy(curve[k])
        crit_ok = all(ct in a.index and a[ct] >= s.acc_resolved for ct in critical if ct in acc.index)
        f = _frac_resolved(a, s)
        curve_points.append(
            {
                "k": int(k),
                "frac_resolved": f,
                "mean_acc": float(a.mean()) if len(a) else None,
                "meets_target": bool(f >= s.frac_resolved_go and crit_ok),
            }
        )
    n_needed = next((p["k"] for p in curve_points if p["meets_target"]), None)
    if not curve_points:
        size_level = None
    elif n_needed is not None and n_needed <= (1 - s.headroom) * budget:
        size_level = GO
    else:
        # Not reaching the target within the budget is a caveat here; the resolution checks above decide failure.
        size_level = CAVEAT
    required_fits = n_required_genes <= budget
    size = {
        "panel_capacity": panel_capacity,
        "reserved_slots": reserved_slots,
        "budget": budget,
        "n_needed": n_needed,
        "headroom_target": int(math.floor((1 - s.headroom) * budget)),
        "level": size_level,
        "n_required_genes": n_required_genes,
        "required_fits": required_fits,
        "curve": curve_points,
    }

    # Secondary metrics: never fail on their own.
    secondary = []

    def gap_vs_pca(prefix: str, label: str) -> None:
        v = _summary_value(summary, SPAPROS_SET_ID, prefix)
        ref = _summary_value(summary, "PCA", prefix)
        if v is None:
            return
        flag = ref is not None and v < ref - s.secondary_gap_vs_pca
        secondary.append(
            {
                "metric": label,
                "value": v,
                "reference": ref,
                "flag": flag,
                "caveat": (
                    (
                        f"{label} of the panel ({_fmt(v)}) trails the PCA baseline ({_fmt(ref)}) by more than "
                        f"{_fmt(s.secondary_gap_vs_pca)}: the panel recovers the labels but loses unlabelled structure"
                    )
                    if flag
                    else None
                ),
            }
        )

    gap_vs_pca("cluster_similarity nmi_5_20", "Coarse clustering similarity (NMI 5-20)")
    gap_vs_pca("cluster_similarity nmi_21_60", "Fine clustering similarity (NMI 21-60)")
    gap_vs_pca("knn_overlap mean_overlap_AUC", "Neighbourhood overlap (kNN AUC)")
    gc = _summary_value(summary, SPAPROS_SET_ID, "gene_corr perct max")
    if gc is not None:
        flag = gc < s.gene_corr_min
        secondary.append(
            {
                "metric": "Non-redundant genes (max correlation < 0.8)",
                "value": gc,
                "reference": s.gene_corr_min,
                "flag": flag,
                "caveat": (
                    f"only {_pct(gc)} of genes are non-redundant; panel slots are spent on near-duplicates"
                    if flag
                    else None
                ),
            }
        )
    mc = _summary_value(summary, SPAPROS_SET_ID, "marker_corr per celltype")
    if mc is None:
        mc = _summary_value(summary, SPAPROS_SET_ID, "marker_corr")
    if mc is not None:
        flag = mc < s.marker_corr_min
        secondary.append(
            {
                "metric": "Marker correlation",
                "value": mc,
                "reference": s.marker_corr_min,
                "flag": flag,
                "caveat": (
                    f"the panel agrees poorly with the curated markers (correlation {_fmt(mc)})" if flag else None
                ),
            }
        )

    # Verdict logic.
    failures: List[str] = []
    caveats: List[str] = []
    notes: List[str] = []
    if panel["frac_resolved_level"] == FAIL:
        failures.append(
            f"only {_pct(frac)} of assessed cell types are resolved (needs at least {_pct(s.frac_resolved_caveat)})"
        )
    if panel["mean_acc_level"] == FAIL:
        failures.append(f"mean cell type recall is {_fmt(mean_acc)} (needs at least {_fmt(s.mean_acc_caveat)})")
    if not required_fits:
        failures.append(f"the {n_required_genes} required genes do not fit the budget of {budget}")
    if critical_unresolved:
        failures.append(f"critical cell type(s) unresolved: {', '.join(critical_unresolved)}")
    if beats_random is False and panel_level != GO:
        failures.append(
            f"the panel is not clearly better than random genes ({_fmt(mean_acc)} vs {_fmt(random_acc)}): the cell "
            "types are hard to separate with any gene set of this size"
        )

    if panel["frac_resolved_level"] == CAVEAT:
        caveats.append(f"{_pct(frac)} of assessed cell types are resolved (target {_pct(s.frac_resolved_go)})")
    if panel["mean_acc_level"] == CAVEAT:
        caveats.append(f"mean cell type recall is {_fmt(mean_acc)} (target {_fmt(s.mean_acc_go)})")
    if critical_marginal:
        caveats.append(f"critical cell type(s) only marginally resolved: {', '.join(critical_marginal)}")
    if size_level == CAVEAT:
        if n_needed is None or n_needed > budget:
            extra = f"; a panel of {n_needed} genes would" if n_needed else ""
            caveats.append(
                f"the resolution target ({_pct(s.frac_resolved_go)} of types and all critical types) is not reached "
                f"within {budget} genes{extra}"
            )
        else:
            caveats.append(
                f"about {n_needed} of {budget} genes are needed, leaving less than {_pct(s.headroom)} headroom"
            )
    for r in celltypes:
        if r["critical"] or r["class"] not in (MARGINAL, UNRESOLVED):
            continue
        why = f", mostly confused with {r['confused_with']} ({_fmt(r['confusion_rate'])})" if r["confused_with"] else ""
        caveats.append(f"{r['celltype']} is {r['class']} (recall {_fmt(r['acc'])}{why})")
    if within_best is False:
        caveats.append(
            f"spapros under-performed the {best_baseline} baseline ({_fmt(mean_acc)} vs "
            f"{_fmt(set_mean_acc[best_baseline])}); check the selection settings and rerun"
        )
    for resc in rescues:
        caveats.append(
            f"{resc['celltype']} is resolved by the {resc['baseline']} baseline (recall {_fmt(resc['acc'])}); adding "
            "its top genes for this type may close the gap"
        )
    caveats += [m["caveat"] for m in secondary if m["flag"]]
    if beats_random is False and panel_level == GO:
        notes.append("the task is easy: random gene sets of this size also separate the cell types well")
    if beats_random is None:
        notes.append("comparison with random gene sets was not run")
    if not curve_points:
        notes.append("the panel size curve was not computed, so the genes actually needed are unknown")

    if input_problems:
        verdict = INCONCLUSIVE
        reason = f"The input cannot support a judgement: {input_problems[0]}."
    elif failures:
        verdict = NOT_FEASIBLE
        reason = f"Not feasible with {budget} genes: {failures[0]}."
    elif caveats:
        verdict = FEASIBLE_WITH_CAVEATS
        reason = (
            f"Workable with {budget} genes ({panel['n_resolved']} of {len(assessed)} cell types resolved), "
            f"but read {'the caveat' if len(caveats) == 1 else f'the {len(caveats)} caveats'} first."
        )
    else:
        verdict = FEASIBLE
        needed = f"; about {n_needed} genes would suffice" if n_needed else ""
        reason = f"All {len(assessed)} assessed cell types are resolved with {budget} genes{needed}."

    return {
        "verdict": verdict,
        "reason": reason,
        "failures": failures,
        "caveats": caveats,
        "notes": notes,
        "critical_celltypes": critical,
        "input_quality": input_quality,
        "panel": panel,
        "size": size,
        "baselines": baselines,
        "secondary": secondary,
        "celltypes": celltypes,
        "thresholds": s.model_dump(),
    }
