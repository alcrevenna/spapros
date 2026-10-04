"""Uploaded datasets, kept apart from jobs so the GUI can inspect a file before a run and reuse it for several runs.

Every upload lives in ``<data_dir>/uploads/<upload_id>/`` with ``input.h5ad`` and ``meta.json`` (what
:func:`inspect_h5ad` found). A job takes a hard link to the file, so deleting an upload never breaks a job.
"""

import json
import os
import shutil
import uuid
from datetime import datetime
from datetime import timezone
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

UPLOAD_FILE = "input.h5ad"
MAX_CATEGORIES = 500


def inspect_h5ad(path: Path) -> Dict[str, Any]:
    """Summarise an ``.h5ad`` for the run form: size, candidate cell type columns with their counts, gene subsets.

    Raises:
        ValueError: With a message meant for the user when the file cannot be read.
    """
    import anndata
    import numpy as np
    import pandas as pd

    try:
        adata = anndata.read_h5ad(path, backed="r")
    except Exception as e:
        raise ValueError(f"could not read the file as .h5ad: {e}") from e
    try:
        obs_columns = []
        for name in adata.obs.columns:
            col = adata.obs[name]
            n_unique = int(col.nunique())
            entry: Dict[str, Any] = {"name": str(name), "n_unique": n_unique, "values": None}
            if 2 <= n_unique <= MAX_CATEGORIES and not pd.api.types.is_float_dtype(col):
                counts = col.value_counts()
                entry["values"] = [{"name": str(k), "count": int(v)} for k, v in counts.items() if v > 0]
            obs_columns.append(entry)
        var_bool_columns = [
            {"name": str(name), "n_true": int(adata.var[name].sum())}
            for name in adata.var.columns
            if adata.var[name].dtype == bool
        ]
        looks_like_counts: Optional[bool] = None
        try:
            sample = adata.X[: min(200, adata.n_obs)]
            values = (
                sample.data if hasattr(sample, "data") and not isinstance(sample, np.ndarray) else np.asarray(sample)
            )
            values = np.asarray(values).ravel()
            looks_like_counts = bool(values.size == 0 or (values.min() >= 0 and np.allclose(values, np.round(values))))
        except Exception:
            pass
        return {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "obs_columns": obs_columns,
            "var_bool_columns": var_bool_columns,
            "looks_like_counts": looks_like_counts,
        }
    finally:
        adata.file.close()


class UploadStore:
    def __init__(self, data_dir: os.PathLike):
        self.root = Path(data_dir) / "uploads"
        self.root.mkdir(parents=True, exist_ok=True)

    def path(self, upload_id: str) -> Path:
        return self.root / upload_id / UPLOAD_FILE

    def create(self, write: Any, filename: str) -> Dict[str, Any]:
        """Store a new upload. ``write(dest)`` writes the file; the upload is removed again if it is not a usable .h5ad."""
        upload_id = uuid.uuid4().hex[:12]
        d = self.root / upload_id
        d.mkdir()
        try:
            write(d / UPLOAD_FILE)
            meta = {
                "id": upload_id,
                "filename": filename,
                "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "size_bytes": (d / UPLOAD_FILE).stat().st_size,
                **inspect_h5ad(d / UPLOAD_FILE),
            }
        except BaseException:
            shutil.rmtree(d)
            raise
        (d / "meta.json").write_text(json.dumps(meta))
        return meta

    def get(self, upload_id: str) -> Optional[Dict[str, Any]]:
        if not upload_id.isalnum():
            return None
        try:
            return json.loads((self.root / upload_id / "meta.json").read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            return None

    def list(self) -> List[Dict[str, Any]]:
        metas = [self.get(p.name) for p in self.root.iterdir() if p.is_dir()]
        return sorted((m for m in metas if m is not None), key=lambda m: m["created_at"])

    def delete(self, upload_id: str) -> None:
        shutil.rmtree(self.root / upload_id)

    def link_into(self, upload_id: str, dest: Path) -> None:
        """Put the upload's file at ``dest``: a hard link when possible, otherwise a copy."""
        try:
            os.link(self.path(upload_id), dest)
        except OSError:
            shutil.copyfile(self.path(upload_id), dest)
