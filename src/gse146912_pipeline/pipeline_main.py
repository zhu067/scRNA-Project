from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc

from gse146912_pipeline.annotation import apply_low_confidence_label, score_major_celltypes
from gse146912_pipeline.config import load_config, merge_config, project_root, resolve_paths
from gse146912_pipeline.event_flow import run_event_flow_prototype
from gse146912_pipeline.flowsig_stage import run_flowsig_stage
from gse146912_pipeline.injury_comm import injury_comm_stage
from gse146912_pipeline.markers import map_marker_symbols_to_ensembl
from gse146912_pipeline.pec import run_pec_subclustering
from gse146912_pipeline.provenance import RunManifest, hash_input_h5ad
from gse146912_pipeline.qc import (
    annotate_mitochondrial_genes,
    apply_qc_filters,
    batch_balance_summary,
    compute_qc_metrics,
    run_scrublet_optional,
)


def _xmax_approx(adata: ad.AnnData, step: int) -> float:
    xmax = 0.0
    for i in range(0, adata.n_obs, step):
        chunk = adata.X[i : i + step]
        if hasattr(chunk, "data") and chunk.data.size:
            xmax = max(xmax, float(chunk.data.max()))
    return xmax


def _run_pec_save(
    adata: ad.AnnData,
    cfg: dict[str, Any],
    rs: int,
    manifest: RunManifest,
    out_dir: Path,
) -> None:
    if not cfg.get("pec_subcluster", {}).get("enabled", True):
        manifest.add_step("pec_subcluster_skipped", extra={"reason": "disabled_in_config"})
        return
    pec_mask = (adata.obs["cell_type_major"] == "PEC").values
    n_pec = int(pec_mask.sum())
    if n_pec < 10:
        manifest.add_step("pec_subcluster_skipped", extra={"reason": "too_few_pec", "n_pec": n_pec})
        return
    Xpec = adata.raw[pec_mask].X
    obs_pec = adata.obs.loc[pec_mask].copy()
    var_full = adata.raw.var.copy()
    adata_pec = ad.AnnData(Xpec, obs=obs_pec, var=var_full)
    adata_pec.var_names_make_unique()
    adata_pec = run_pec_subclustering(adata_pec, cfg, rs)
    out_pec = out_dir / cfg.get("outputs", {}).get("pec_h5ad", "pec.h5ad")
    adata_pec.write_h5ad(out_pec, compression="gzip")
    manifest.add_step(
        "pec_subcluster",
        n_cells_before=adata.n_obs,
        n_cells_after=n_pec,
        extra={"path": str(out_pec)},
    )


def run_pipeline(
    config_path: Path,
    local_config_path: Path | None,
    stages: list[str],
) -> None:
    cfg = load_config(config_path)
    if local_config_path and local_config_path.is_file():
        cfg = merge_config(cfg, load_config(local_config_path))

    root = project_root(cfg)
    cfg = resolve_paths(cfg, root)
    rs = int(cfg.get("project", {}).get("random_seed", 0))
    np.random.seed(rs)
    sc.settings.seed = rs

    R = cfg["_resolved"]
    out_dir = Path(R["output_run_dir"])
    fig_dir = Path(R["figures_dir"])
    tab_dir = Path(R["tables_dir"])
    prov_dir = Path(R["provenance_dir"])
    for d in (out_dir, fig_dir, tab_dir, prov_dir):
        d.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(fig_dir)

    manifest = RunManifest(prov_dir)
    manifest.set_config_snapshot(cfg)
    prov_cfg = cfg.get("provenance", {})

    if "all" in stages:
        stages = ["main", "pec_subcluster", "injury_comm", "event_flow", "flowsig"]
    stages_set = set(stages)

    manifest.add_step(
        "init",
        parameters={"random_seed": rs, "stages": list(stages_set), "config_path": str(config_path)},
    )

    out_main = out_dir / cfg.get("outputs", {}).get("analyzed_h5ad", "GSE146912_merged_analyzed.h5ad")

    # 仅下游：按需加载 analyzed；纯 event_flow 仅需 R 产出的 CellChat 差异 CSV
    if "main" not in stages_set:
        needs_h5ad = bool(stages_set & {"pec_subcluster", "injury_comm"})
        adata = None
        if needs_h5ad:
            if not out_main.is_file():
                raise FileNotFoundError(
                    f"未包含 main 阶段且缺少 {out_main}，请先运行 --stages main"
                )
            manifest.set_input_fingerprint(
                hash_input_h5ad(out_main, bool(prov_cfg.get("hash_input_full_file", False)))
            )
            adata = sc.read_h5ad(out_main)
            manifest.add_step("load_analyzed_h5ad", n_cells_after=adata.n_obs, n_vars=adata.n_vars)
        else:
            manifest.add_step("load_analyzed_skipped", extra={"reason": "stages_do_not_need_h5ad"})
        if adata is not None and "pec_subcluster" in stages_set:
            _run_pec_save(adata, cfg, rs, manifest, out_dir)
        if (
            adata is not None
            and "injury_comm" in stages_set
            and cfg.get("injury_comm", {}).get("enabled", False)
        ):
            try:
                res = injury_comm_stage(adata, cfg, tab_dir, manifest)
                manifest.add_step("injury_comm", extra=res)
            except ImportError as e:
                manifest.add_step("injury_comm_skipped", extra={"error": str(e)})
        if "event_flow" in stages_set and cfg.get("event_flow", {}).get("enabled", False):
            try:
                ef_res = run_event_flow_prototype(cfg, tab_dir, manifest)
                manifest.add_step("event_flow", extra=ef_res)
            except FileNotFoundError as e:
                manifest.add_step("event_flow_skipped", extra={"error": str(e)})
        if "flowsig" in stages_set and cfg.get("flowsig", {}).get("enabled", False):
            try:
                fs_res = run_flowsig_stage(cfg, tab_dir, manifest)
                manifest.add_step("flowsig", extra=fs_res)
            except (FileNotFoundError, ImportError, ValueError) as e:
                manifest.add_step("flowsig_skipped", extra={"error": str(e)})
        manifest._flush()
        print("完成（下游阶段）。输出目录:", out_dir)
        print("过程记录:", manifest.path)
        return

    raw_path = Path(R["raw_h5ad"])
    if not raw_path.is_file():
        raise FileNotFoundError(f"输入 h5ad 不存在: {raw_path}")
    manifest.set_input_fingerprint(
        hash_input_h5ad(raw_path, bool(prov_cfg.get("hash_input_full_file", False)))
    )

    adata = sc.read_h5ad(raw_path)
    n_load = adata.n_obs
    manifest.add_step("load_h5ad", n_cells_after=n_load, n_vars=adata.n_vars)

    norm_cfg = cfg.get("normalization", {})
    xmax = _xmax_approx(adata, int(norm_cfg.get("xmax_sample_step", 8000)))
    skip_norm = xmax < float(norm_cfg.get("xmax_skip_normalize", 100))
    manifest.write_step_file(
        "normalization_heuristic",
        {"xmax_approx": xmax, "skip_normalize_total_log1p": skip_norm},
    )
    if skip_norm:
        manifest.write_step_file(
            "paper_replication_warning",
            {
                "level": "warning",
                "message": (
                    "Skipped normalize_total + log1p (input X appears already log/scaled per heuristic). "
                    "For Liu WB et al. Kidney Int 2023 parity, verify raw-count vs preprocessed "
                    "input against paper Methods; see docs/METHODS_PAPER_ALIGNMENT.md. "
                    "SCENIC/pySCENIC inputs must follow that tool's requirements separately."
                ),
                "xmax_approx": xmax,
            },
        )

    qc_cfg = cfg.get("qc", {})
    if qc_cfg.get("annotate_mito", True):
        mres = annotate_mitochondrial_genes(
            adata,
            org=qc_cfg.get("mito_biomart_org", "mmusculus"),
            attrname=qc_cfg.get("mito_attrname", "ensembl_gene_id"),
        )
        manifest.add_step("mitochondrial_annotation", extra=mres)
    compute_qc_metrics(adata)
    manifest.add_step("calculate_qc_metrics", n_cells_after=adata.n_obs)

    adata, frep = apply_qc_filters(
        adata,
        qc_cfg.get("filter_min_genes"),
        qc_cfg.get("filter_max_genes"),
        qc_cfg.get("filter_min_counts"),
        qc_cfg.get("filter_max_mito_pct"),
    )
    manifest.add_step(
        "qc_filter",
        n_cells_before=frep["n_before"],
        n_cells_after=frep["n_after"],
        extra=frep,
    )

    bb = batch_balance_summary(
        adata,
        cfg.get("processing", {}).get("batch_key", "batch"),
        float(cfg.get("batch_balance", {}).get("min_fraction_warn", 0.02)),
    )
    manifest.write_step_file("batch_balance", bb)

    sb = qc_cfg.get("scrublet", {})
    if sb.get("enabled", False):
        run_scrublet_optional(
            adata,
            float(sb.get("expected_doublet_rate", 0.06)),
            rs,
            manifest,
            assume_log_normalized=skip_norm,
        )

    if not skip_norm:
        sc.pp.normalize_total(adata, target_sum=float(norm_cfg.get("target_sum", 1e4)))
        sc.pp.log1p(adata)
        manifest.add_step("normalize_total_log1p", parameters={"target_sum": norm_cfg.get("target_sum")})
    adata.raw = adata.copy()
    manifest.add_step("set_raw", n_cells_after=adata.n_obs, n_vars=adata.raw.shape[1])

    proc = cfg.get("processing", {})
    batch_key = proc.get("batch_key", "batch")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=int(proc.get("hvg_n_top_genes", 3000)),
        batch_key=batch_key if batch_key in adata.obs else None,
        flavor=proc.get("hvg_flavor", "seurat"),
        subset=bool(proc.get("hvg_subset", True)),
    )
    sc.pp.scale(
        adata,
        max_value=float(proc.get("scale_max_value", 10)),
        zero_center=bool(proc.get("scale_zero_center", False)),
    )
    sc.tl.pca(adata, svd_solver=proc.get("pca_solver", "arpack"))
    sc.pp.neighbors(
        adata,
        n_neighbors=int(proc.get("n_neighbors", 15)),
        n_pcs=int(proc.get("n_pcs_neighbors", 40)),
    )
    if proc.get("umap", True):
        sc.tl.umap(adata)
    ld = cfg.get("leiden", {})
    sc.tl.leiden(
        adata,
        resolution=float(ld.get("resolution", 0.5)),
        key_added=ld.get("key_added", "leiden"),
        flavor=ld.get("flavor", "igraph"),
        n_iterations=int(ld.get("n_iterations", 2)),
    )
    manifest.add_step(
        "embedding_clustering",
        n_cells_after=adata.n_obs,
        parameters={
            "hvg_n_top_genes": proc.get("hvg_n_top_genes"),
            "leiden_resolution": ld.get("resolution"),
            "n_pcs_neighbors": proc.get("n_pcs_neighbors"),
        },
        extra={"n_clusters": int(adata.obs[ld.get("key_added", "leiden")].nunique())},
    )

    ann = cfg.get("annotation", {})
    marker_sets = ann.get("marker_sets", {})
    marker_ens, mstats = map_marker_symbols_to_ensembl(
        marker_sets, species=ann.get("mygene_species", "mouse")
    )
    manifest.write_step_file("marker_symbol_to_ensembl", mstats)

    astats = score_major_celltypes(adata, marker_ens, use_raw=True)
    manifest.write_step_file("cell_type_scoring", astats)
    apply_low_confidence_label(adata, float(ann.get("score_margin_lowconf", 0.03)))
    manifest.add_step(
        "cell_type_major",
        n_cells_after=adata.n_obs,
        extra={"value_counts": adata.obs["cell_type_major"].value_counts().to_dict()},
    )

    adata.write_h5ad(out_main, compression="gzip")
    manifest.add_step("save_analyzed", extra={"path": str(out_main)})

    if "pec_subcluster" in stages_set:
        _run_pec_save(adata, cfg, rs, manifest, out_dir)

    if "injury_comm" in stages_set and cfg.get("injury_comm", {}).get("enabled", False):
        try:
            res = injury_comm_stage(adata, cfg, tab_dir, manifest)
            manifest.add_step("injury_comm", extra=res)
        except ImportError as e:
            manifest.add_step("injury_comm_skipped", extra={"error": str(e)})

    if "event_flow" in stages_set and cfg.get("event_flow", {}).get("enabled", False):
        try:
            ef_res = run_event_flow_prototype(cfg, tab_dir, manifest)
            manifest.add_step("event_flow", extra=ef_res)
        except FileNotFoundError as e:
            manifest.add_step("event_flow_skipped", extra={"error": str(e)})
    if "flowsig" in stages_set and cfg.get("flowsig", {}).get("enabled", False):
        try:
            fs_res = run_flowsig_stage(cfg, tab_dir, manifest)
            manifest.add_step("flowsig", extra=fs_res)
        except (FileNotFoundError, ImportError, ValueError) as e:
            manifest.add_step("flowsig_skipped", extra={"error": str(e)})

    manifest._flush()
    print("完成。输出目录:", out_dir)
    print("过程记录:", manifest.path)


def main_cli() -> None:
    ap = argparse.ArgumentParser(description="GSE146912 可复现主流程")
    ap.add_argument(
        "--config",
        type=Path,
        default=Path("config/default.yaml"),
        help="主配置文件路径",
    )
    ap.add_argument(
        "--local-config",
        type=Path,
        default=None,
        help="可选本地覆盖（如 config/local.yaml）",
    )
    ap.add_argument(
        "--stages",
        nargs="+",
        default=["main"],
        help="main | pec_subcluster | injury_comm | event_flow | flowsig | all",
    )
    args = ap.parse_args()
    run_pipeline(args.config, args.local_config, args.stages)


if __name__ == "__main__":
    main_cli()
