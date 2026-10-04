"""Self-contained HTML feasibility report for one spapros run.

The report is a single HTML file with inline SVG figures, so it can be viewed in the GUI, downloaded, mailed and
printed to PDF without anything else.
"""

import html
import io
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

import pandas as pd

from spapros.server import feasibility as fz

VERDICT_STYLE = {
    fz.FEASIBLE: "go",
    fz.FEASIBLE_WITH_CAVEATS: "caveat",
    fz.NOT_FEASIBLE: "fail",
    fz.INCONCLUSIVE: "unknown",
}
CLASS_STYLE = {fz.RESOLVED: "go", fz.MARGINAL: "caveat", fz.UNRESOLVED: "fail", fz.NOT_ASSESSED: "unknown"}
LEVEL_STYLE = {fz.GO: "go", fz.CAVEAT: "caveat", fz.FAIL: "fail", None: "unknown"}

# Colours shared by the figures; chosen to read on both the light and the dark report background.
INK = "#6b7280"
ACCENT = "#2563eb"
SECOND = "#d97706"


def _e(x: Any) -> str:
    return html.escape("" if x is None else str(x))


def _num(x: Optional[float], digits: int = 2) -> str:
    return "–" if x is None or pd.isna(x) else f"{x:.{digits}f}"


def _pct(x: Optional[float]) -> str:
    return "–" if x is None or pd.isna(x) else f"{100 * x:.0f}%"


def _badge(text: str, style: str) -> str:
    return f'<span class="badge {style}">{_e(text)}</span>'


# Figures -------------------------------------------------------------------------------------------------------------


def _svg(fig) -> str:
    import matplotlib.pyplot as plt

    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    svg = buf.getvalue()
    return svg[svg.index("<svg") :]


def _style_axes(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK)
    ax.tick_params(colors=INK, labelcolor=INK)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def confusion_figure(confusion: pd.DataFrame, acc_resolved: float) -> str:
    """Heatmap of the spapros confusion matrix; diagonal cells below the resolved cut-off are outlined."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["svg.fonttype"] = "none"
    cm = confusion.rename(index=str, columns=str)
    n = len(cm)
    size = min(max(3.5, 0.45 * n + 1.5), 12)
    fig, ax = plt.subplots(figsize=(size + 1, size))
    im = ax.imshow(cm.values, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(cm.columns)), cm.columns, rotation=90)
    ax.set_yticks(range(n), cm.index)
    ax.set_xlabel("predicted cell type")
    ax.set_ylabel("true cell type")
    if n <= 25:
        for i in range(n):
            for j in range(len(cm.columns)):
                v = cm.values[i, j]
                if v >= 0.01:
                    ax.text(
                        j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white" if v > 0.6 else "black"
                    )
    for i, ct in enumerate(cm.index):
        if ct in cm.columns and cm.loc[ct, ct] < acc_resolved:
            j = list(cm.columns).index(ct)
            ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#dc2626", linewidth=2))
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(colors=INK, labelcolor=INK)
    _style_axes(ax)
    return _svg(fig)


def curve_figure(curve: List[Dict[str, Any]], budget: int, n_needed: Optional[int], frac_go: float) -> str:
    """Share of resolved cell types and mean recall against panel size."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["svg.fonttype"] = "none"
    ks = [p["k"] for p in curve]
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.plot(ks, [p["frac_resolved"] for p in curve], "o-", color=ACCENT, label="cell types resolved")
    ax.plot(ks, [p["mean_acc"] for p in curve], "s--", color=SECOND, label="mean recall")
    ax.axhline(frac_go, color=INK, linewidth=0.8, linestyle=":")
    ax.axvline(budget, color=INK, linewidth=1)
    ax.text(budget, 1.04, f" budget {budget}", color=INK, fontsize=8, ha="left", va="bottom")
    if n_needed is not None:
        ax.axvline(n_needed, color=ACCENT, linewidth=1, linestyle="--")
        ax.text(n_needed, -0.1, f"needed ≈{n_needed} ", color=ACCENT, fontsize=8, ha="right", va="top")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("genes in panel (top ranked)")
    ax.set_ylabel("score")
    ax.legend(frameon=False, loc="lower right", labelcolor=INK)
    _style_axes(ax)
    return _svg(fig)


def baseline_figure(mean_acc: Dict[str, Optional[float]]) -> str:
    """Mean recall of spapros against the baseline gene sets."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["svg.fonttype"] = "none"
    items = sorted(((k, v) for k, v in mean_acc.items() if v is not None), key=lambda kv: kv[0] != fz.SPAPROS_SET_ID)
    fig, ax = plt.subplots(figsize=(6.5, 0.4 * len(items) + 1))
    names = [k for k, _ in items][::-1]
    values = [v for _, v in items][::-1]
    colors = [ACCENT if k == fz.SPAPROS_SET_ID else "#9ca3af" for k in names]
    ax.barh(names, values, color=colors)
    for i, v in enumerate(values):
        ax.text(v + 0.01, i, f"{v:.2f}", va="center", fontsize=8, color=INK)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("mean cell type recall")
    _style_axes(ax)
    return _svg(fig)


# HTML ----------------------------------------------------------------------------------------------------------------

CSS = """
:root { --bg:#ffffff; --fg:#111827; --muted:#6b7280; --line:#e5e7eb; --panel:#f9fafb;
  --go:#15803d; --go-bg:#dcfce7; --caveat:#a16207; --caveat-bg:#fef3c7; --fail:#b91c1c; --fail-bg:#fee2e2;
  --unknown:#4b5563; --unknown-bg:#e5e7eb; --accent:#2563eb; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#0f1115; --fg:#e5e7eb; --muted:#9ca3af;
  --line:#2a2f3a; --panel:#161a22; --go:#4ade80; --go-bg:#12301d; --caveat:#facc15; --caveat-bg:#3a2f0b;
  --fail:#f87171; --fail-bg:#3b1414; --unknown:#d1d5db; --unknown-bg:#2a2f3a; --accent:#60a5fa; } }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--fg); font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:1000px; margin:0 auto; padding:24px 16px 64px; }
h1 { font-size:22px; margin:0 0 4px; } h2 { font-size:17px; margin:32px 0 8px; border-bottom:1px solid var(--line);
  padding-bottom:4px; }
.muted { color:var(--muted); } a { color:var(--accent); }
.verdict { border-radius:10px; padding:16px 18px; margin:16px 0; border:1px solid var(--line); }
.verdict .label { font-size:20px; font-weight:700; letter-spacing:.02em; }
.verdict.go { background:var(--go-bg); } .verdict.go .label { color:var(--go); }
.verdict.caveat { background:var(--caveat-bg); } .verdict.caveat .label { color:var(--caveat); }
.verdict.fail { background:var(--fail-bg); } .verdict.fail .label { color:var(--fail); }
.verdict.unknown { background:var(--unknown-bg); } .verdict.unknown .label { color:var(--unknown); }
.badge { display:inline-block; padding:1px 8px; border-radius:999px; font-size:12px; font-weight:600; white-space:nowrap; }
.badge.go { background:var(--go-bg); color:var(--go); } .badge.caveat { background:var(--caveat-bg); color:var(--caveat); }
.badge.fail { background:var(--fail-bg); color:var(--fail); }
.badge.unknown { background:var(--unknown-bg); color:var(--unknown); }
table { border-collapse:collapse; width:100%; font-size:14px; }
th, td { text-align:left; padding:6px 8px; border-bottom:1px solid var(--line); vertical-align:top; }
th { color:var(--muted); font-weight:600; } td.num, th.num { text-align:right; font-variant-numeric:tabular-nums; }
.scroll { overflow-x:auto; }
.grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr)); gap:8px 24px; }
.kv { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:10px 12px; }
.kv .k { color:var(--muted); font-size:12px; } .kv .v { font-size:16px; font-weight:600; }
figure { margin:12px 0; } figure svg { max-width:100%; height:auto; }
ul.items li { margin:4px 0; }
details summary { cursor:pointer; color:var(--accent); }
@media print { main { max-width:none; } details { display:block; } a { color:inherit; } }
"""


def _kv(label: str, value: Any) -> str:
    return f'<div class="kv"><div class="k">{_e(label)}</div><div class="v">{_e(value)}</div></div>'


def _list(items: List[str], style: str) -> str:
    if not items:
        return ""
    return '<ul class="items">' + "".join(f"<li>{_badge(style.upper(), style)} {_e(i)}</li>" for i in items) + "</ul>"


def render_report(
    verdict: Dict[str, Any],
    run: Dict[str, Any],
    confusion: pd.DataFrame,
    probeset: Optional[pd.DataFrame] = None,
    files: Optional[List[str]] = None,
) -> str:
    """Render the report.

    Args:
        verdict: Output of :func:`spapros.server.feasibility.assess`.
        run: Run description: ``dataset``, ``n_cells``, ``n_genes``, ``celltype_key``, ``spapros_version``,
            ``options``, ``runtime_seconds``, ``created``.
        confusion: Confusion matrix of the spapros gene set.
        probeset: ``selector.probeset``; the selected genes are listed with rank and marker cell types.
        files: Result files of the job, linked relative to the report's URL (``files/<name>``).
    """
    v = verdict
    th = v["thresholds"]
    size = v["size"]
    panel = v["panel"]
    base = v["baselines"]
    style = VERDICT_STYLE[v["verdict"]]
    opts = run.get("options", {})
    out: List[str] = []
    w = out.append

    w(
        "<!doctype html><html lang=en><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
    )
    w(f"<title>Feasibility report: {_e(run.get('dataset'))}</title><style>{CSS}</style></head><body><main>")
    w(f"<h1>Panel feasibility report</h1><div class=muted>{_e(run.get('dataset'))} · {_e(run.get('created'))}</div>")

    # 1. Verdict banner.
    needed = f", needed ≈ {size['n_needed']}" if size.get("n_needed") else ""
    w(f'<section class="verdict {style}"><div class=label>{_e(v["verdict"])}</div><div>{_e(v["reason"])}</div>')
    w(
        f"<div class=muted>Gene budget B = {size['budget']} of panel capacity C = {size['panel_capacity']}"
        f" ({size['reserved_slots']} reserved){_e(needed)}</div></section>"
    )
    if v["input_quality"]["problems"]:
        w("<h2>Input problems</h2>" + _list(v["input_quality"]["problems"], "fail"))
    if v["failures"] or v["caveats"] or v["notes"]:
        w("<h2>Why</h2>")
        w(_list(v["failures"], "fail") + _list(v["caveats"], "caveat"))
        if v["notes"]:
            w('<ul class="items">' + "".join(f"<li class=muted>{_e(n)}</li>" for n in v["notes"]) + "</ul>")

    # 2. Run summary.
    w("<h2>Run summary</h2><div class=grid>")
    w(_kv("Cells", run.get("n_cells")) + _kv("Genes in data", run.get("n_genes")))
    w(_kv("Cell types (assessed)", f"{v['input_quality']['n_celltypes']} ({v['input_quality']['n_assessed']})"))
    w(_kv("Cell type column", run.get("celltype_key")) + _kv("Genes selected", size["budget"]))
    w(_kv("Critical cell types", ", ".join(v["critical_celltypes"]) or "none"))
    w(
        _kv(
            "Resolved / marginal / unresolved",
            f"{panel['n_resolved']} / {panel['n_marginal']} / {panel['n_unresolved']}",
        )
    )
    w(_kv("Mean recall", _num(panel["mean_acc"])) + _kv("Share resolved", _pct(panel["frac_resolved"])))
    runtime = run.get("runtime_seconds")
    w(_kv("Runtime", f"{runtime / 60:.1f} min" if runtime else "–") + _kv("spapros", run.get("spapros_version")))
    w("</div>")

    # 3. Cell type table.
    w(
        "<h2>Cell types</h2><p class=muted>Recall is the share of a type's cells the classifier assigns correctly using "
        f"only the panel genes. Resolved ≥ {th['acc_resolved']:.2f}, marginal ≥ {th['acc_marginal']:.2f}. Sorted worst "
        "first.</p><div class=scroll><table><tr><th>Cell type</th><th class=num>Cells</th><th class=num>Recall</th>"
        "<th>Class</th><th>Mostly confused with</th><th>Best baseline</th></tr>"
    )
    for r in v["celltypes"]:
        crit = " " + _badge("critical", "unknown") if r["critical"] else ""
        conf = f"{_e(r['confused_with'])} ({_num(r['confusion_rate'])})" if r["confused_with"] else ""
        bb = f"{_e(r['best_baseline'])} ({_num(r['best_baseline_acc'])})" if r["best_baseline"] else ""
        w(
            f"<tr><td>{_e(r['celltype'])}{crit}</td><td class=num>{_e(r['n_cells'])}</td>"
            f"<td class=num>{_num(r['acc'])}</td><td>{_badge(r['class'], CLASS_STYLE[r['class']])}</td>"
            f"<td>{conf}</td><td>{bb}</td></tr>"
        )
    w("</table></div>")

    # 4. Confusion matrix.
    w(
        "<h2>Confusion matrix</h2><p class=muted>Rows are true cell types, columns the classifier's prediction. "
        "Outlined diagonal cells are below the resolved cut-off.</p>"
    )
    w(f"<figure>{confusion_figure(confusion, th['acc_resolved'])}</figure>")

    # 5. Panel size curve.
    w("<h2>Panel size</h2>")
    if size["curve"]:
        w(
            "<p class=muted>The top-ranked genes of the spapros list evaluated at several panel sizes. The target is "
            f"{_pct(th['frac_resolved_go'])} of cell types and every critical type resolved, with "
            f"{_pct(th['headroom'])} of the budget to spare (≤ {size['headroom_target']} genes).</p>"
        )
        w(f"<figure>{curve_figure(size['curve'], size['budget'], size['n_needed'], th['frac_resolved_go'])}</figure>")
        w(
            "<div class=scroll><table><tr><th class=num>Genes</th><th class=num>Resolved</th><th class=num>Mean recall"
            "</th><th>Meets target</th></tr>"
        )
        for p in size["curve"]:
            ok = _badge("yes", "go") if p["meets_target"] else _badge("no", "unknown")
            w(
                f"<tr><td class=num>{p['k']}</td><td class=num>{_pct(p['frac_resolved'])}</td>"
                f"<td class=num>{_num(p['mean_acc'])}</td><td>{ok}</td></tr>"
            )
        w("</table></div>")
    else:
        w("<p class=muted>Not computed for this run.</p>")
    req = "fit" if size["required_fits"] else "do not fit"
    w(
        f"<p>Genes the panel must contain (preselected genes and the markers spapros requires per cell type): "
        f"{size['n_required_genes']}, which {req} the budget.</p>"
    )

    # 6. Baselines.
    w("<h2>Comparison with simple gene sets</h2>")
    if len(base["mean_acc"]) > 1:
        w(
            "<p class=muted>PCA, DE (differential expression), HVG (highly variable) and random gene sets of the same "
            "size, evaluated the same way.</p>"
        )
        w(f"<figure>{baseline_figure(base['mean_acc'])}</figure><ul class=items>")
        if base["beats_random"] is not None:
            ok = base["beats_random"]
            w(
                f"<li>{_badge('pass' if ok else 'no', 'go' if ok else 'caveat')} Beats random genes by at least "
                f"{th['beats_random_margin']:.2f} (random mean {_num(base['random_mean_acc'])})</li>"
            )
        if base["within_best_baseline"] is not None:
            ok = base["within_best_baseline"]
            w(
                f"<li>{_badge('pass' if ok else 'no', 'go' if ok else 'caveat')} Within {th['baseline_tolerance']:.2f} "
                f"of the best simple baseline ({_e(base['best_baseline'])})</li>"
            )
        w("</ul>")
    else:
        w("<p class=muted>Baseline gene sets were not evaluated for this run.</p>")
    if base["skipped"]:
        w("<p class=muted>Skipped: " + "; ".join(f"{_e(k)} ({_e(m)})" for k, m in base["skipped"].items()) + "</p>")

    # 7. Secondary metrics.
    w("<h2>Secondary metrics</h2>")
    if v["secondary"]:
        w(
            "<p class=muted>These never fail a panel on their own but can add a caveat.</p><div class=scroll><table>"
            "<tr><th>Metric</th><th class=num>Panel</th><th class=num>Reference</th><th>Status</th></tr>"
        )
        for m in v["secondary"]:
            w(
                f"<tr><td>{_e(m['metric'])}</td><td class=num>{_num(m['value'])}</td>"
                f"<td class=num>{_num(m['reference'])}</td>"
                f"<td>{_badge('caveat', 'caveat') if m['flag'] else _badge('ok', 'go')}</td></tr>"
            )
        w("</table></div>")
    else:
        w("<p class=muted>None available for this run.</p>")
    w("<p class=muted>Expression fit to the platform's detection range is not checked.</p>")

    # 8. Gene panel.
    w("<h2>Gene panel</h2>")
    if files:
        w("<p>Downloads: " + " · ".join(f'<a href="files/{_e(f)}" download>{_e(f)}</a>' for f in files) + "</p>")
    if probeset is not None and "selection" in probeset.columns:
        sel = probeset[probeset["selection"].astype(bool)].copy()
        if "rank" in sel.columns:
            sel = sel.sort_values("rank")
        marker_col = next((c for c in ("marker_celltypes", "celltypes_marker") if c in sel.columns), None)
        w(
            f"<details><summary>{len(sel)} selected genes</summary><div class=scroll><table><tr><th>Gene</th>"
            "<th class=num>Rank</th><th>Marker for</th></tr>"
        )
        for gene, row in sel.iterrows():
            marker = row[marker_col] if marker_col and isinstance(row[marker_col], str) else ""
            rank = row["rank"] if "rank" in sel.columns else ""
            w(f"<tr><td>{_e(gene)}</td><td class=num>{_e(rank)}</td><td>{_e(marker)}</td></tr>")
        w("</table></div></details>")

    # 9. Methods.
    w(
        "<h2>Methods and limits</h2><p>Recalls come from spapros' XGBoost classifier (5-fold cross validation × 5 seeds, "
        "balanced class weights) trained on the panel genes only. Cell types with fewer than 40 cells are not assessed. "
        "Classification on single-cell RNA-seq is an upper bound for spatial data, where capture is lower and cell "
        "segmentation adds errors. Probe design is not checked.</p>"
    )
    w("<details><summary>Thresholds used</summary><div class=scroll><table>")
    for k, val in th.items():
        w(f"<tr><td>{_e(k)}</td><td class=num>{_e(val)}</td></tr>")
    w("</table></div></details>")
    w("<details><summary>Run options</summary><div class=scroll><table>")
    for k, val in opts.items():
        if k != "feasibility":
            w(f"<tr><td>{_e(k)}</td><td>{_e(val)}</td></tr>")
    w("</table></div></details>")
    w("</main></body></html>")
    return "".join(out)
