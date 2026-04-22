from __future__ import annotations

from typing import Any

import mygene


def map_marker_symbols_to_ensembl(
    marker_sets: dict[str, list[str]],
    species: str = "mouse",
    batch_size: int = 200,
) -> tuple[dict[str, list[str]], dict[str, Any]]:
    """将各细胞类型的 marker symbol 转为 Ensembl ID；返回 (marker_ens, 统计信息)。"""
    flat = sorted({g for genes in marker_sets.values() for g in genes})
    mg = mygene.MyGeneInfo()
    sym_to_ens: dict[str, str] = {}
    for i in range(0, len(flat), batch_size):
        chunk = flat[i : i + batch_size]
        res = mg.querymany(
            chunk, scopes="symbol", fields="ensembl.gene", species=species, verbose=False
        )
        for r in res:
            q = r.get("query")
            eg = r.get("ensembl")
            if isinstance(eg, list) and eg:
                eg = eg[0]
            if isinstance(eg, dict):
                gid = eg.get("gene")
                if gid and q:
                    sym_to_ens[q] = gid

    missing = sorted(set(flat) - set(sym_to_ens))
    marker_ens = {
        k: [sym_to_ens[s] for s in v if s in sym_to_ens] for k, v in marker_sets.items()
    }
    per_type = {
        k: {
            "n_genes_requested": len(v),
            "n_genes_mapped": len(marker_ens[k]),
            "coverage": len(marker_ens[k]) / max(len(v), 1),
        }
        for k, v in marker_sets.items()
    }
    stats: dict[str, Any] = {
        "species": species,
        "n_symbols_requested": len(flat),
        "n_symbols_mapped_to_ensembl": len(sym_to_ens),
        "symbol_mapping_coverage": len(sym_to_ens) / max(len(flat), 1),
        "missing_symbols": missing,
        "per_cell_type": per_type,
    }
    return marker_ens, stats
