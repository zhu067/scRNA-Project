from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import anndata as ad
import mygene
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix, csr_matrix
import scipy.sparse as sp


@dataclass
class CollapseReport:
    n_genes_x: int = 0
    n_symbols_unique: int = 0
    n_raw_genes: int = 0
    used_separate_raw_matrix: bool = False
    extra: dict[str, Any] = field(default_factory=dict)


def collapse_ensembl_to_symbol(
    adata_obj: ad.AnnData,
    species: str = "mouse",
    mygene_batch: int = 2500,
) -> tuple[ad.AnnData, CollapseReport]:
    """Ensembl → symbol，同 symbol 列对表达取平均；raw 与 X 列数不同时单独构造 Mr。"""
    rep = CollapseReport()
    ids_x = adata_obj.var_names.astype(str).tolist()
    ids_r = None
    if adata_obj.raw is not None and adata_obj.raw.shape[1] != adata_obj.shape[1]:
        ids_r = adata_obj.raw.var_names.astype(str).tolist()
        rep.used_separate_raw_matrix = True

    all_ids = list(dict.fromkeys(ids_x + (ids_r or [])))
    mg = mygene.MyGeneInfo()
    id_to_sym: dict[str, str] = {}
    for i in range(0, len(all_ids), mygene_batch):
        chunk = all_ids[i : i + mygene_batch]
        res = mg.querymany(
            chunk, scopes="ensembl.gene", fields="symbol", species=species, verbose=False
        )
        for r in res:
            id_to_sym[r["query"]] = r.get("symbol") or r["query"]

    sym_x = np.array([id_to_sym.get(i, i) for i in ids_x], dtype=str)
    u = np.unique(sym_x)
    sym_r_arr: np.ndarray | None = None
    if ids_r is not None:
        sym_r_arr = np.array([id_to_sym.get(i, i) for i in ids_r], dtype=str)
        u = np.unique(np.concatenate([sym_x, sym_r_arr]))
    sym_to_idx = {s: j for j, s in enumerate(u)}
    n_new = len(u)
    inv_x = np.array([sym_to_idx[s] for s in sym_x], dtype=int)

    X = adata_obj.X
    if not sp.issparse(X):
        X = csr_matrix(X)
    if X.shape[1] != len(inv_x):
        raise ValueError(f"X columns {X.shape[1]} != inv_x {len(inv_x)}")

    cx = np.bincount(inv_x, minlength=n_new)
    wx = 1.0 / np.maximum(cx[inv_x], 1)
    Mx = coo_matrix(
        (wx, (np.arange(len(inv_x), dtype=int), inv_x)), shape=(X.shape[1], n_new)
    ).tocsc()
    Xn = X @ Mx
    out = ad.AnnData(X=Xn, obs=adata_obj.obs.copy(), var=pd.DataFrame(index=u))

    if adata_obj.raw is not None:
        Xr = adata_obj.raw.X
        if not sp.issparse(Xr):
            Xr = csr_matrix(Xr)
        if Xr.shape[1] == X.shape[1]:
            out.raw = ad.AnnData(X=Xr @ Mx, obs=out.obs, var=out.var.copy())
        else:
            assert sym_r_arr is not None
            inv_r = np.array([sym_to_idx[s] for s in sym_r_arr], dtype=int)
            cr = np.bincount(inv_r, minlength=n_new)
            wr = 1.0 / np.maximum(cr[inv_r], 1)
            Mr = coo_matrix(
                (wr, (np.arange(len(inv_r), dtype=int), inv_r)), shape=(Xr.shape[1], n_new)
            ).tocsc()
            out.raw = ad.AnnData(X=Xr @ Mr, obs=out.obs, var=out.var.copy())

    rep.n_genes_x = len(ids_x)
    rep.n_raw_genes = len(ids_r) if ids_r else adata_obj.raw.shape[1] if adata_obj.raw else 0
    rep.n_symbols_unique = n_new
    rep.extra = {
        "n_ensembl_ids_queried": len(all_ids),
        "n_ids_with_symbol": len(id_to_sym),
        "id_to_symbol_coverage": len(id_to_sym) / max(len(all_ids), 1),
    }
    return out, rep
