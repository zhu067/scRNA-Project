from __future__ import annotations

from typing import Any

import numpy as np
import scanpy as sc

from gse146912_pipeline.provenance import RunManifest


def annotate_mitochondrial_genes(
    adata: sc.AnnData,
    org: str = "mmusculus",
    attrname: str = "ensembl_gene_id",
) -> dict[str, Any]:
    """尝试用 scanpy Biomart 查询线粒体基因；失败则返回 status=skipped。"""
    try:
        mito = sc.queries.mitochondrial_genes(org, attrname=attrname)
        col = mito.columns[0]
        mito_ids = set(mito[col].astype(str))
        adata.var["mt"] = adata.var_names.astype(str).isin(mito_ids)
        n_mito = int(adata.var["mt"].sum())
        return {
            "status": "ok",
            "organism": org,
            "attrname": attrname,
            "n_mito_genes_in_var": n_mito,
            "n_mito_in_biomart_query": len(mito_ids),
        }
    except Exception as e:
        adata.var["mt"] = False
        return {
            "status": "skipped",
            "error": f"{type(e).__name__}: {e}",
            "hint": "Install pybiomart or run online; mito QC metrics will be zeros.",
        }


def compute_qc_metrics(adata: sc.AnnData, mito_key: str = "mt") -> None:
    # Scanpy 若干版本在 qc_vars=None 时会报错，故仅在存在线粒体标记时传入
    if mito_key in adata.var.columns and adata.var[mito_key].any():
        sc.pp.calculate_qc_metrics(adata, qc_vars=[mito_key])
    else:
        sc.pp.calculate_qc_metrics(adata)


def apply_qc_filters(
    adata: sc.AnnData,
    min_genes: int | None,
    max_genes: int | None,
    min_counts: int | None,
    max_mito_pct: float | None,
) -> tuple[sc.AnnData, dict[str, Any]]:
    """按阈值过滤细胞；无阈值则原样返回。"""
    n0 = adata.n_obs
    report: dict[str, Any] = {"n_before": n0, "filters_applied": []}
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
        report["filters_applied"].append({"min_genes": min_genes})
    if max_genes is not None:
        adata = adata[adata.obs["n_genes_by_counts"] <= max_genes]
        report["filters_applied"].append({"max_genes": max_genes})
    if min_counts is not None:
        adata = adata[adata.obs["total_counts"] >= min_counts]
        report["filters_applied"].append({"min_counts": min_counts})
    if max_mito_pct is not None and "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] <= max_mito_pct]
        report["filters_applied"].append({"max_mito_pct": max_mito_pct})
    report["n_after"] = adata.n_obs
    report["cells_removed"] = n0 - adata.n_obs
    return adata, report


def batch_balance_summary(
    adata: sc.AnnData, batch_key: str, min_fraction_warn: float
) -> dict[str, Any]:
    vc = adata.obs[batch_key].value_counts(normalize=True)
    small = vc[vc < min_fraction_warn]
    return {
        "batch_key": batch_key,
        "n_batches": int(vc.size),
        "counts": adata.obs[batch_key].value_counts().to_dict(),
        "fractions": vc.to_dict(),
        "batches_below_min_fraction_warn": small.to_dict(),
    }


def run_scrublet_optional(
    adata: sc.AnnData,
    expected_doublet_rate: float,
    random_seed: int,
    manifest: RunManifest | None,
    assume_log_normalized: bool,
) -> dict[str, Any]:
    """可选双细胞检测。若数据已 log 归一化，结果仅作参考。"""
    out: dict[str, Any] = {"ran": False}
    if assume_log_normalized:
        out["skipped"] = "data treated as log-normalized (xmax heuristic); scrublet unreliable"
        return out
    try:
        import scrublet as scr
    except ImportError:
        out["skipped"] = "scrublet not installed"
        return out

    counts = adata.X
    if hasattr(counts, "toarray"):
        counts = counts.toarray()
    else:
        counts = np.asarray(counts)
    scrub = scr.Scrublet(counts, expected_doublet_rate=expected_doublet_rate, random_state=random_seed)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs["scrublet_score"] = doublet_scores
    adata.obs["scrublet_predicted_doublet"] = predicted_doublets
    out["ran"] = True
    out["predicted_doublet_rate"] = float(np.mean(predicted_doublets))
    if manifest:
        manifest.add_step(
            "scrublet",
            n_cells_before=adata.n_obs,
            n_cells_after=adata.n_obs,
            parameters={"expected_doublet_rate": expected_doublet_rate, "random_seed": random_seed},
            extra=out,
        )
    return out
