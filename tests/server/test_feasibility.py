import numpy as np
import pandas as pd
import pytest

from spapros.server import feasibility as fz
from spapros.server.options import FeasibilitySettings
from spapros.server.options import RunOptions
from spapros.server.report import render_report

TYPES = ["A", "B", "C", "D", "E"]


def cm(diag, types=TYPES, partner=None):
    """Row-normalised confusion matrix with the given diagonal; the rest of each row goes to ``partner`` or the next type."""
    m = pd.DataFrame(0.0, index=types, columns=types)
    for i, (t, d) in enumerate(zip(types, diag)):
        m.loc[t, t] = d
        other = partner.get(t) if partner and t in partner else types[(i + 1) % len(types)]
        m.loc[t, other] += 1 - d
    return m


def summary(**mean_acc):
    return pd.DataFrame({"forest_clfs accuracy": mean_acc}).rename_axis("set_id")


COUNTS = {t: 500 for t in TYPES}


def run(diags, counts=COUNTS, curve=None, critical=None, budget=100, **kw):
    """Assess with spapros at ``diags`` and baselines a bit worse, random much worse."""
    conf = {"spapros": cm(diags)}
    conf.setdefault("PCA", cm([max(0, d - 0.05) for d in diags]))
    conf.setdefault("random_seed0", cm([0.3] * len(diags)))
    conf.update(kw.pop("confusion", {}))
    if curve is None:
        curve = {25: cm([0.5] * 5), 50: cm([0.95] * 5), 75: cm([0.95] * 5), 100: conf["spapros"]}
    return fz.assess(
        cell_counts=counts,
        confusion=conf,
        summary=kw.pop("summary", None),
        budget=budget,
        panel_capacity=budget,
        critical_celltypes=critical,
        curve=curve,
        **kw,
    )


def test_feasible():
    v = run([0.95, 0.9, 0.92, 0.97, 0.88])
    assert v["verdict"] == fz.FEASIBLE, v
    assert v["size"]["n_needed"] == 50
    assert v["size"]["level"] == fz.GO
    assert v["panel"]["n_resolved"] == 5
    assert "about 50 genes" in v["reason"]


def test_marginal_type_is_a_named_caveat():
    v = run([0.95, 0.9, 0.92, 0.97, 0.7])
    assert v["verdict"] == fz.FEASIBLE_WITH_CAVEATS
    assert any(c.startswith("E is marginal") and "confused with A" in c for c in v["caveats"])
    assert v["celltypes"][0]["celltype"] == "E"  # worst first


def test_unresolved_critical_type_fails():
    v = run([0.95, 0.9, 0.92, 0.97, 0.5], critical=["E"])
    assert v["verdict"] == fz.NOT_FEASIBLE
    assert "critical cell type(s) unresolved: E" in v["failures"]
    assert v["celltypes"][0]["critical"]


def test_marginal_critical_type_is_a_caveat():
    v = run([0.95, 0.9, 0.92, 0.97, 0.7], critical=["E"])
    assert v["verdict"] == fz.FEASIBLE_WITH_CAVEATS
    assert any("only marginally resolved: E" in c for c in v["caveats"])


def test_low_share_resolved_fails():
    v = run([0.5, 0.5, 0.9, 0.95, 0.5])
    assert v["verdict"] == fz.NOT_FEASIBLE
    assert v["panel"]["frac_resolved_level"] == fz.FAIL


def test_tight_budget_is_a_caveat():
    curve = {25: cm([0.5] * 5), 50: cm([0.7] * 5), 75: cm([0.7] * 5), 100: cm([0.95] * 5)}
    v = run([0.95] * 5, curve=curve)
    assert v["verdict"] == fz.FEASIBLE_WITH_CAVEATS
    assert v["size"]["n_needed"] == 100
    assert any("about 100 of 100 genes" in c for c in v["caveats"])


def test_target_only_reached_by_bigger_panel():
    curve = {50: cm([0.5] * 5), 100: cm([0.8, 0.8, 0.8, 0.8, 0.7]), 125: cm([0.9] * 5)}
    v = run([0.8, 0.8, 0.8, 0.8, 0.7], curve=curve)
    assert v["size"]["n_needed"] == 125
    assert any("a panel of 125 genes would" in c for c in v["caveats"])


def test_too_few_cells_is_inconclusive():
    counts = {**COUNTS, "F": 10, "G": 12}
    v = run([0.95] * 5, counts=counts)
    assert v["verdict"] == fz.INCONCLUSIVE
    assert v["input_quality"]["not_assessed"] == ["F", "G"]
    assert {r["celltype"]: r["class"] for r in v["celltypes"]}["F"] == fz.NOT_ASSESSED


def test_critical_type_not_assessed_is_inconclusive():
    counts = {**COUNTS, "F": 10}
    v = run([0.95] * 5, counts=counts, critical=["F"])
    assert v["verdict"] == fz.INCONCLUSIVE
    assert any("critical" in p for p in v["input_quality"]["problems"])


def test_no_better_than_random_fails_when_not_go():
    diags = [0.8, 0.8, 0.8, 0.65, 0.8]
    v = run(diags, confusion={"random_seed0": cm(diags), "random_seed1": cm(diags)})
    assert v["verdict"] == fz.NOT_FEASIBLE
    assert any("not clearly better than random" in f for f in v["failures"])


def test_easy_task_is_only_a_note():
    diags = [0.95] * 5
    v = run(diags, confusion={"random_seed0": cm(diags)})
    assert v["verdict"] == fz.FEASIBLE
    assert any("task is easy" in n for n in v["notes"])


def test_baseline_rescue_and_underperformance():
    diags = [0.95, 0.95, 0.95, 0.95, 0.5]
    v = run(diags, confusion={"DE": cm([0.99] * 5)})
    assert any("resolved by the DE baseline" in c for c in v["caveats"])
    assert any("under-performed the DE baseline" in c for c in v["caveats"])
    assert v["celltypes"][0]["best_baseline"] == "DE"


def test_required_genes_must_fit():
    v = run([0.95] * 5, n_required_genes=120)
    assert v["verdict"] == fz.NOT_FEASIBLE


def test_secondary_metrics_add_caveats():
    s = pd.DataFrame(
        {
            "forest_clfs accuracy": {"spapros": 0.93, "PCA": 0.88, "random_seed0": 0.3},
            "knn_overlap mean_overlap_AUC": {"spapros": 0.2, "PCA": 0.5, "random_seed0": 0.1},
            "gene_corr perct max < 0.8": {"spapros": 0.8, "PCA": 0.9, "random_seed0": 1.0},
        }
    )
    v = run([0.95, 0.9, 0.92, 0.97, 0.88], summary=s)
    assert v["verdict"] == fz.FEASIBLE_WITH_CAVEATS
    flagged = [m["metric"] for m in v["secondary"] if m["flag"]]
    assert flagged == ["Neighbourhood overlap (kNN AUC)", "Non-redundant genes (max correlation < 0.8)"]


def test_settings_change_the_verdict():
    diags = [0.95, 0.9, 0.92, 0.97, 0.78]
    assert run(diags)["verdict"] == fz.FEASIBLE_WITH_CAVEATS
    assert run(diags, settings=FeasibilitySettings(acc_resolved=0.75))["verdict"] == fz.FEASIBLE


def test_curve_sizes():
    assert fz.curve_sizes(300, 2000) == [75, 150, 225, 300, 375, 450]
    assert fz.curve_sizes(300, 320) == [75, 150, 225, 300]


def test_options_budget():
    assert RunOptions(celltype_key="ct").n == 300
    assert RunOptions(celltype_key="ct", panel_capacity=140, reserved_slots=20).n == 120
    assert RunOptions(celltype_key="ct", n=50).n == 50
    with pytest.raises(ValueError):
        RunOptions(celltype_key="ct", panel_capacity=100, n=150)
    with pytest.raises(ValueError):
        RunOptions(celltype_key="ct", panel_capacity=10, reserved_slots=10)
    with pytest.raises(ValueError):
        FeasibilitySettings(acc_resolved=0.5, acc_marginal=0.6)


def test_report_renders_and_escapes():
    v = run([0.95, 0.9, 0.92, 0.97, 0.7])
    probeset = pd.DataFrame({"selection": [True, True, False], "rank": [1, 2, 3]}, index=["G1", "<G2>", "G3"])
    html = render_report(
        v,
        {"dataset": "data<1>.h5ad", "n_cells": 2500, "options": {"celltype_key": "ct"}},
        cm([0.95, 0.9, 0.92, 0.97, 0.7]),
        probeset=probeset,
        files=["probeset.csv"],
    )
    assert "FEASIBLE WITH CAVEATS" in html
    assert "data&lt;1&gt;.h5ad" in html and "&lt;G2&gt;" in html and "<G2>" not in html
    assert html.count("<svg") == 3  # confusion matrix, panel size curve, baselines
    assert 'href="files/probeset.csv"' in html
    assert np.isfinite(v["panel"]["mean_acc"])
