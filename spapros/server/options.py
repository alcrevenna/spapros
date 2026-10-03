"""Run options accepted by the spapros server."""

from typing import Any
from typing import Dict
from typing import List
from typing import Literal
from typing import Optional

from pydantic import BaseModel
from pydantic import ConfigDict
from pydantic import Field


class RunOptions(BaseModel):
    """Options for one probe set selection + evaluation run.

    The fields map onto :class:`spapros.se.ProbesetSelector` and :class:`spapros.ev.ProbesetEvaluator`. Anything not
    listed here keeps the spapros default.
    """

    model_config = ConfigDict(extra="forbid")

    celltype_key: str = Field(..., description="Column in adata.obs with the cell type annotation.")
    n: int = Field(50, ge=1, description="Number of genes to select.")
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
    seed: int = 0
    n_jobs: int = Field(-1, description="CPUs per run; -1 uses all.")
