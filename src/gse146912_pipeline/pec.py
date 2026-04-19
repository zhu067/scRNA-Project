from __future__ import annotations

from typing import Any

import numpy as np
import scanpy as sc


def pec_subset_from_adata(
    adata: sc.AnnData,
    pec_label: str = "PEC",
    cell_type_col: str = "cell_type_major",
) -> tuple[sc.AnnData, np.ndarray]:
    mask = (adata.obs[cell_type_col] == pec_label).values
    return adata[mask].copy(), mask


def run_pec_subclustering(
    adata_pec_fullgenes: sc.AnnData,
    cfg: dict[str, Any],
    random_seed: int,
) -> sc.AnnData:
    """在 PEC 对象上（已为全基因或 symbol 空间）做 HVG、scale、PCA、neighbors、UMAP、多分辨率 Leiden。"""
    pc = cfg.get("pec_subcluster", {})
    batch_key = cfg.get("processing", {}).get("batch_key", "batch")
    bk = batch_key if batch_key in adata_pec_fullgenes.obs else None

    sc.pp.highly_variable_genes(
        adata_pec_fullgenes,
        n_top_genes=pc.get("hvg_n_top_genes", 2500),
        batch_key=bk,
        flavor=cfg.get("processing", {}).get("hvg_flavor", "seurat"),
        subset=False,
    )
    adata_pec_fullgenes.raw = adata_pec_fullgenes.copy()
    adata_hvg = adata_pec_fullgenes[:, adata_pec_fullgenes.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=pc.get("scale_max_value", 10), zero_center=False)
    sc.tl.pca(adata_hvg, svd_solver=cfg.get("processing", {}).get("pca_solver", "arpack"))
    sc.pp.neighbors(
        adata_hvg,
        n_neighbors=pc.get("n_neighbors", 15),
        n_pcs=min(
            pc.get("n_pcs_max", 40),
            adata_hvg.obsm["X_pca"].shape[1],
        ),
    )
    sc.tl.umap(adata_hvg)
    ld = cfg.get("leiden", {})
    flavor = ld.get("flavor", "igraph")
    n_it = ld.get("n_iterations", 2)
    for key, res in (pc.get("leiden_resolutions") or {}).items():
        sc.tl.leiden(
            adata_hvg,
            resolution=res,
            key_added=f"leiden_pec_{key}",
            flavor=flavor,
            n_iterations=n_it,
        )
    default_key = pc.get("default_cluster_key", "leiden_pec_r08")
    return adata_hvg
