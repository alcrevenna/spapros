"""Run options accepted by the spapros server."""

from typing import Any
from typing import Dict
from typing import List
from typing import Literal
from typing import Optional

from pydantic import BaseModel
from pydantic import ConfigDict
from pydantic import Field
from pydantic import model_validator


class FeasibilitySettings(BaseModel):
    """Thresholds of the feasibility verdict written at the end of each run.

    The defaults follow the project's feasibility rule (v1). Accuracies are per cell type recalls, the diagonal of the
    ``forest_clfs`` confusion matrix.
    """

    model_config = ConfigDict(extra="forbid")

    acc_resolved: float = Field(0.80, ge=0, le=1, description="A cell type with recall at or above this is resolved.")
    acc_marginal: float = Field(0.60, ge=0, le=1, description="Below acc_resolved but at or above this is marginal.")
    frac_resolved_go: float = Field(0.90, ge=0, le=1, description="Share of resolved cell types for a clear go.")
    frac_resolved_caveat: float = Field(0.70, ge=0, le=1, description="Below this share the panel fails.")
    mean_acc_go: float = Field(0.85, ge=0, le=1, description="Mean recall over cell types for a clear go.")
    mean_acc_caveat: float = Field(0.75, ge=0, le=1, description="Below this mean recall the panel fails.")
    beats_random_margin: float = Field(0.05, ge=0, le=1, description="How much spapros must beat random gene sets by.")
    baseline_tolerance: float = Field(0.02, ge=0, le=1, description="How far spapros may trail the best baseline.")
    headroom: float = Field(0.20, ge=0, lt=1, description="Share of the gene budget that should be left unused.")
    secondary_gap_vs_pca: float = Field(
        0.10, ge=0, le=1, description="Caveat when clustering or kNN recovery trails the PCA baseline by more."
    )
    gene_corr_min: float = Field(0.90, ge=0, le=1, description="Caveat when fewer genes are non-redundant.")
    marker_corr_min: float = Field(0.50, ge=0, le=1, description="Caveat when marker correlation is lower.")
    max_excluded_types: float = Field(
        0.20, ge=0, le=1, description="Inconclusive above this share of unassessed types."
    )
    max_excluded_cells: float = Field(
        0.05, ge=0, le=1, description="Inconclusive above this share of unassessed cells."
    )

    @model_validator(mode="after")
    def _ordered(self) -> "FeasibilitySettings":
        for low, high in [
            ("acc_marginal", "acc_resolved"),
            ("frac_resolved_caveat", "frac_resolved_go"),
            ("mean_acc_caveat", "mean_acc_go"),
        ]:
            if getattr(self, low) > getattr(self, high):
                raise ValueError(f"{low} must not be larger than {high}")
        return self


class RunOptions(BaseModel):
    """Options for one probe set selection + evaluation run.

    The fields map onto :class:`spapros.se.ProbesetSelector` and :class:`spapros.ev.ProbesetEvaluator`. Anything not
    listed here keeps the spapros default.
    """

    model_config = ConfigDict(extra="forbid")

    celltype_key: str = Field(..., description="Column in adata.obs with the cell type annotation.")
    panel_capacity: int = Field(300, ge=1, description="Genes the spatial platform's panel holds.")
    reserved_slots: int = Field(
        0, ge=0, description="Panel slots kept for genes added outside spapros (controls, genes of interest)."
    )
    n: Optional[int] = Field(
        None,
        ge=1,
        description="Number of genes to select. Defaults to the gene budget panel_capacity - reserved_slots.",
    )
    critical_celltypes: List[str] = Field(
        default_factory=list,
        description="Cell types that must be resolved; any unresolved one makes the run infeasible.",
    )
    genes_key: Optional[str] = Field(
        "highly_variable",
        description="Boolean column in adata.var restricting the genes to select from. None uses all genes.",
    )
    normalize: bool = Field(
        True, description="Run sc.pp.normalize_total + sc.pp.log1p first. Disable if adata.X is already log-normalised."
    )
    preselected_genes: List[str] = Field(default_factory=list)
    prior_genes: List[str] = Field(default_factory=list)
    n_pca_genes: int = Field(100, ge=0)
    n_min_markers: int = Field(2, ge=0)
    forest_hparams: Optional[Dict[str, Any]] = Field(
        None, description="Overrides for the selector's forest hyperparameters (n_trees, subsample, test_subsample)."
    )
    evaluation_scheme: Literal["quick", "full"] = "quick"
    evaluate_reference_sets: bool = Field(
        True, description="Also select and evaluate PCA/DE/HVG/random baseline sets of the same size for comparison."
    )
    panel_size_curve: bool = Field(
        True, description="Evaluate smaller and larger prefixes of the ranked gene list to estimate the genes needed."
    )
    seed: int = 0
    n_jobs: int = Field(-1, description="CPUs per run; -1 uses all.")
    feasibility: FeasibilitySettings = Field(default_factory=FeasibilitySettings)

    @model_validator(mode="after")
    def _budget(self) -> "RunOptions":
        budget = self.panel_capacity - self.reserved_slots
        if budget < 1:
            raise ValueError("reserved_slots must leave at least one gene for spapros")
        if self.n is None:
            self.n = budget
        elif self.n > budget:
            raise ValueError(f"n={self.n} is larger than the gene budget panel_capacity - reserved_slots = {budget}")
        return self
