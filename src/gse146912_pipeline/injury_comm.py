from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

import pandas as pd
import scanpy as sc

from gse146912_pipeline.genes import collapse_ensembl_to_symbol
from gse146912_pipeline.injury import assign_injury_group


def prepare_injury_comparison_adata(
    adata: sc.AnnData,
    control_token: str,
    nephritis_token: str,
) -> tuple[sc.AnnData, dict[str, Any]]:
    adata = adata.copy()
    adata.obs["injury_group"] = assign_injury_group(
        adata.obs["Sample"], control_token, nephritis_token
    )
    mask = adata.obs["injury_group"].isin(["control", "nephritis"])
    sub = adata[mask].copy()
    stats = {
        "n_cells_total": adata.n_obs,
        "n_cells_injury_comparison": sub.n_obs,
        "injury_group_counts": sub.obs["injury_group"].value_counts().to_dict(),
    }
    return sub, stats


def run_liana_rank_aggregate(
    adata_sym: sc.AnnData,
    groupby: str,
    resource_name: str,
    expr_prop: float,
    key_added: str = "liana_global",
) -> pd.DataFrame:
    import liana as li

    li.mt.rank_aggregate(
        adata_sym,
        groupby=groupby,
        resource_name=resource_name,
        expr_prop=expr_prop,
        use_raw=False,
        key_added=key_added,
        verbose=False,
    )
    return adata_sym.uns[key_added].copy()


def injury_comm_stage(
    adata: sc.AnnData,
    cfg: dict[str, Any],
    tables_dir: Path,
    provenance: Any,
) -> dict[str, Any]:
    """肾炎 vs 对照：symbol 矩阵导出供 CellChat（R）；可选 LIANA（非论文主流程）。"""
    ic = cfg.get("injury_comm", {})
    sub, st = prepare_injury_comparison_adata(
        adata,
        ic.get("control_token", "control"),
        ic.get("nephritis_token", "nephritis"),
    )
    adata_sym, rep = collapse_ensembl_to_symbol(sub)
    provenance.write_step_file(
        "injury_symbol_collapse",
        {"collapse": asdict(rep), "injury_stats": st},
    )
    h5 = Path(cfg["_resolved"]["output_run_dir"]) / cfg.get("outputs", {}).get(
        "cellchat_export_h5ad", "exp.h5ad"
    )
    h5.parent.mkdir(parents=True, exist_ok=True)
    adata_sym.write_h5ad(h5, compression="gzip")
    comm_dir = tables_dir / "PEC_injury_comm"
    comm_dir.mkdir(parents=True, exist_ok=True)
    ctab = pd.crosstab(
        adata_sym.obs["injury_group"],
        adata_sym.obs["cell_type_major"],
    )
    counts_csv = comm_dir / "injury_group_cell_type_major_counts.csv"
    ctab.to_csv(counts_csv)
    out: dict[str, Any] = {
        "export_h5ad": str(h5),
        "injury_celltype_counts_csv": str(counts_csv),
        "paper_primary": "Kidney International 2023 (GSE146912): CellChat (R) + SCENIC/pySCENIC; see docs/METHODS_PAPER_ALIGNMENT.md",
    }

    li_cfg = ic.get("liana", {})
    if not bool(li_cfg.get("enabled", False)):
        provenance.write_step_file(
            "liana_skipped",
            {
                "reason": "injury_comm.liana.enabled is false (default aligns with paper CellChat-first)",
            },
        )
        return out

    df = run_liana_rank_aggregate(
        adata_sym,
        groupby="cell_type_major",
        resource_name=li_cfg.get("resource_name", "mouseconsensus"),
        expr_prop=float(li_cfg.get("expr_prop", 0.1)),
    )
    out_csv = tables_dir / "PEC_injury_comm" / "LIANA_global.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    out["liana_rows"] = len(df)
    out["liana_csv"] = str(out_csv)
    return out
