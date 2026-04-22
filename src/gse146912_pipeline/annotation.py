from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc


def score_major_celltypes(
    adata: sc.AnnData,
    marker_ens: dict[str, list[str]],
    use_raw: bool = True,
    min_genes_per_set: int = 2,
) -> dict[str, Any]:
    """对 adata 做各类型 score_genes；写入 obs cell_type_major 等。返回统计。"""
    applied = []
    skipped = []
    for ctype, ens_ids in marker_ens.items():
        in_data = [e for e in ens_ids if e in (adata.raw.var_names if use_raw and adata.raw else adata.var_names)]
        if len(in_data) < min_genes_per_set:
            skipped.append({"type": ctype, "reason": "too_few_genes_in_matrix", "n": len(in_data)})
            continue
        sc.tl.score_genes(adata, gene_list=in_data, score_name=f"score_{ctype}", use_raw=use_raw)
        applied.append(ctype)

    score_cols = [f"score_{k}" for k in applied if f"score_{k}" in adata.obs]
    if not score_cols:
        raise RuntimeError("无任何细胞类型成功打分，请检查 raw 与 marker Ensembl 是否匹配")

    score_mat = adata.obs[score_cols].copy()
    adata.obs["cell_type_major"] = score_mat.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["cell_type_max_score"] = score_mat.max(axis=1)
    adata.obs["cell_type_second_score"] = score_mat.apply(
        lambda r: float(np.sort(r.values)[-2]), axis=1
    )
    adata.obs["cell_type_score_margin"] = (
        adata.obs["cell_type_max_score"] - adata.obs["cell_type_second_score"]
    )
    return {
        "cell_types_scored": applied,
        "skipped": skipped,
        "score_columns": score_cols,
    }


def apply_low_confidence_label(
    adata: sc.AnnData, margin_threshold: float, out_col: str = "cell_type_major_filter"
) -> None:
    adata.obs[out_col] = adata.obs["cell_type_major"].copy()
    adata.obs.loc[adata.obs["cell_type_score_margin"] < margin_threshold, out_col] = "Low_confidence"
