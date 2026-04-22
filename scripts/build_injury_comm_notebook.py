"""Emit sc-RNAGES146912_PEC_injury_cellcomm.ipynb. Run from repo root: python scripts/build_injury_comm_notebook.py"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "sc-RNAGES146912_PEC_injury_cellcomm.ipynb"

cells = []


def md(s):
    lines = [ln + "\n" for ln in s.strip().split("\n")]
    cells.append({"cell_type": "markdown", "metadata": {}, "source": lines})


def code(s):
    lines = s.strip().split("\n")
    src = [ln + "\n" for ln in lines[:-1]] + ([lines[-1] + "\n"] if lines else [])
    cells.append({"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": src})


md(
    r"""
## 新月体 / 增生性肾炎：PEC 激活与健康对照比较 + 细胞通讯（GSE146912）

**数据说明**：本数据集（Chung 等，肾小球 scRNA）中 `*_nephritis` 为 **肾毒性血清肾炎（NTS）** 模型，常用于新月体/免疫介导损伤；**与临床抗 GBM（抗 IV 型胶原 NC1）并非同一造模**。若需严格「正常 vs 抗 GBM」，请换用含该模型的数据；此处按 **`control` vs `nephritis`（NTS）** 实现，分组逻辑可直接迁移。

**与刘文斌等 Kidney International 2023（GSE146912；壁层上皮细胞亚群与新月体形成等）一致的主结论路径**：① **TF / regulon：SCENIC**（本仓库用 **`docs/METHODS_PAPER_ALIGNMENT.md`** + **`scripts/pyscenic_mouse_GSE146912.sh`**，勿以 DoRothEA+ULM 替代）；② **细胞通讯：CellChat（R）**（**`scripts/cellchat_GSE146912_PEC.R`** + 本笔记导出的 `h5ad`）。③ 本笔记内 **PAGA + DPT**、**PEC DEG** 等仍为 Scanpy 探索流程。

**可选（非论文必需）**：将下方 `EXPLORATORY_PYTHON_MULTIMETHOD = True` 时，才运行 **DoRothEA+ULM** 与 **LIANA**（需 `pip install decoupler liana`），仅作方法比较，**不能**声称与原文 SCENIC/CellChat Figure 一一对应。

**依赖**：`scanpy`、`mygene`（首格自动 pip）；`decoupler`/`liana` 仅探索模式需要。
"""
)

code(
    r"""
import os
import sys
import subprocess

def pip_install(pkgs):
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *pkgs])

for pkg in ("mygene",):
    try:
        __import__(pkg)
    except ImportError:
        pip_install([pkg])

import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT = Path(".").resolve()
OUTPUT_ROOT = Path(os.environ.get("SCRNA_OUTPUT_DIR", PROJECT / "results"))
OUT_DIR = OUTPUT_ROOT / "GSE146912_run"
FIG_DIR = OUT_DIR / "figures" / "PEC_injury_comm"
TAB_DIR = OUT_DIR / "tables" / "PEC_injury_comm"
for d in (FIG_DIR, TAB_DIR):
    d.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIG_DIR)
sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)

# 论文主路径：pySCENIC + R CellChat。仅在为 True 时安装并运行 DoRothEA / LIANA。
EXPLORATORY_PYTHON_MULTIMETHOD = False

ANALYZED = OUT_DIR / "GSE146912_merged_analyzed.h5ad"

if not ANALYZED.is_file():
    raise FileNotFoundError(
        f"未找到 {ANALYZED}，请先运行 sc-RNAGES146912.ipynb 完成注释并保存。"
    )

adata = sc.read_h5ad(ANALYZED)
if "cell_type_major" not in adata.obs.columns:
    raise KeyError("adata 缺少 cell_type_major，请先运行主笔记本注释单元。")

print(adata)
print(adata.obs["cell_type_major"].value_counts())
"""
)

code(
    r"""
def injury_group_from_sample(s: str) -> str:
    s = str(s).lower()
    if "control" in s:
        return "control"
    if "nephritis" in s:
        return "nephritis"
    return "other"


adata.obs["injury_group"] = adata.obs["Sample"].map(injury_group_from_sample)
keep = adata.obs["injury_group"].isin(["control", "nephritis"])
adata_cmp = adata[keep].copy()
print("对照/肾炎子集:", adata_cmp.shape)
print(adata_cmp.obs["injury_group"].value_counts())
print(adata_cmp.obs.groupby(["injury_group", "cell_type_major"]).size().unstack(fill_value=0))
"""
)

code(
    r"""
def collapse_ensembl_to_symbol(adata_obj: ad.AnnData, species: str = "mouse") -> ad.AnnData:
    # Ensembl 转 symbol，同 symbol 列取平均；raw 与 X 基因数不同时单独构造折叠矩阵。
    import mygene
    import scipy.sparse as sp
    from scipy.sparse import coo_matrix

    ids_x = adata_obj.var_names.astype(str).tolist()
    ids_r = None
    if adata_obj.raw is not None and adata_obj.raw.shape[1] != adata_obj.shape[1]:
        ids_r = adata_obj.raw.var_names.astype(str).tolist()

    all_ids = list(dict.fromkeys(ids_x + (ids_r or [])))
    mg = mygene.MyGeneInfo()
    id_to_sym = {}
    step = 2500
    for i in range(0, len(all_ids), step):
        chunk = all_ids[i : i + step]
        res = mg.querymany(
            chunk, scopes="ensembl.gene", fields="symbol", species=species, verbose=False
        )
        for r in res:
            id_to_sym[r["query"]] = r.get("symbol") or r["query"]

    sym_x = np.array([id_to_sym.get(i, i) for i in ids_x], dtype=str)
    u = np.unique(sym_x)
    if ids_r is not None:
        sym_r = np.array([id_to_sym.get(i, i) for i in ids_r], dtype=str)
        u = np.unique(np.concatenate([sym_x, sym_r]))
    sym_to_idx = {s: j for j, s in enumerate(u)}
    n_new = len(u)
    inv_x = np.array([sym_to_idx[s] for s in sym_x], dtype=int)

    X = adata_obj.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    assert X.shape[1] == len(inv_x), (X.shape, len(inv_x))
    cx = np.bincount(inv_x, minlength=n_new)
    wx = 1.0 / np.maximum(cx[inv_x], 1)
    Mx = coo_matrix((wx, (np.arange(len(inv_x), dtype=int), inv_x)), shape=(X.shape[1], n_new)).tocsc()
    Xn = X @ Mx
    out = ad.AnnData(X=Xn, obs=adata_obj.obs.copy(), var=pd.DataFrame(index=u))

    if adata_obj.raw is not None:
        Xr = adata_obj.raw.X
        if not sp.issparse(Xr):
            Xr = sp.csr_matrix(Xr)
        if Xr.shape[1] == X.shape[1]:
            out.raw = ad.AnnData(X=Xr @ Mx, obs=out.obs, var=out.var.copy())
        else:
            inv_r = np.array([sym_to_idx[s] for s in sym_r], dtype=int)
            assert Xr.shape[1] == len(inv_r), (Xr.shape, len(inv_r))
            cr = np.bincount(inv_r, minlength=n_new)
            wr = 1.0 / np.maximum(cr[inv_r], 1)
            Mr = coo_matrix(
                (wr, (np.arange(len(inv_r), dtype=int), inv_r)), shape=(Xr.shape[1], n_new)
            ).tocsc()
            out.raw = ad.AnnData(X=Xr @ Mr, obs=out.obs, var=out.var.copy())
    return out


adata_sym = collapse_ensembl_to_symbol(adata_cmp)
print("symbol 矩阵:", adata_sym.shape)
"""
)

code(
    r"""
# --- A. PEC：肾炎 vs 对照 差异基因 + 增殖模块分 ---
pec = adata_cmp[adata_cmp.obs["cell_type_major"] == "PEC"].copy()
if pec.n_obs < 30:
    raise RuntimeError("PEC 细胞过少。")

sc.tl.rank_genes_groups(
    pec,
    groupby="injury_group",
    reference="control",
    method="wilcoxon",
    use_raw=True,
    key_added="deg_nephritis_vs_control",
)
deg = sc.get.rank_genes_groups_df(pec, group="nephritis", key="deg_nephritis_vs_control")
deg_path = TAB_DIR / "PEC_DEG_nephritis_vs_control.csv"
deg.to_csv(deg_path, index=False)
print("DEG:", deg_path)

import mygene

mg = mygene.MyGeneInfo()
prolif = ["Mki67", "Top2a", "Pcna", "Cdk1", "Ccnb1"]
pr = mg.querymany(prolif, scopes="symbol", fields="ensembl.gene", species="mouse", verbose=False)
ens = []
for r in pr:
    eg = r.get("ensembl")
    if isinstance(eg, list) and eg:
        eg = eg[0]
    if isinstance(eg, dict) and eg.get("gene"):
        e = eg["gene"]
        if e in pec.raw.var_names:
            ens.append(e)
if len(ens) >= 2:
    sc.tl.score_genes(pec, gene_list=ens, score_name="proliferation_score", use_raw=True)
    sc.pl.violin(
        pec,
        keys=["proliferation_score"],
        groupby="injury_group",
        save="_PEC_proliferation.pdf",
        show=False,
    )
    plt.show()
else:
    print("增殖标志基因匹配不足，跳过。")
"""
)

code(
    r"""
# --- B. TF / regulon：论文为 SCENIC；此处默认不跑 DoRothEA（与 SCENIC 不等价）---
if not EXPLORATORY_PYTHON_MULTIMETHOD:
    print(
        "跳过 DoRothEA+ULM。与 Kidney Int 2023 一致的 TF 结论请运行：\n"
        "  docs/METHODS_PAPER_ALIGNMENT.md\n"
        "  scripts/pyscenic_mouse_GSE146912.sh"
    )
else:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "decoupler"])
    import mygene
    import decoupler as dc
    import requests
    from urllib.error import HTTPError, URLError

    pec_sym = adata_sym[adata_sym.obs["cell_type_major"] == "PEC"].copy()

    def dorothea_mouse_from_human_homologene(levels=("A", "B")):
        net_h = dc.op.dorothea(organism="human", levels=list(levels))
        sym_u = list(dict.fromkeys(net_h["source"].tolist() + net_h["target"].tolist()))
        mg = mygene.MyGeneInfo()
        res = mg.querymany(sym_u, scopes="symbol", fields="homologene", species="human", verbose=False)
        mouse_ids = []
        human_keys = []
        for r in res:
            hg = r.get("homologene")
            if not hg or "genes" not in hg:
                continue
            for g in hg["genes"]:
                if g[0] == 10090:
                    human_keys.append(r["query"])
                    mouse_ids.append(str(g[1]))
                    break
        if not mouse_ids:
            raise RuntimeError("HomoloGene 中未解析到小鼠同源基因")
        m2 = mg.querymany(mouse_ids, scopes="entrezgene", fields="symbol", species="mouse", verbose=False)
        id2sym = {str(x["query"]): x.get("symbol") for x in m2}
        h2m = {h: id2sym.get(mid) for h, mid in zip(human_keys, mouse_ids) if id2sym.get(mid)}
        net = net_h.copy()
        net["source"] = net["source"].map(h2m)
        net["target"] = net["target"].map(h2m)
        net = net.dropna(subset=["source", "target"])
        gset = set(pec_sym.var_names.astype(str))
        net = net[net["target"].isin(gset)]
        net = net.drop_duplicates(subset=["source", "target", "weight"])
        return net

    try:
        net = dc.op.dorothea(organism="mouse", levels=["A", "B"])
    except (HTTPError, URLError, OSError, requests.exceptions.HTTPError, requests.exceptions.RequestException) as exc:
        print("mouse DoRothEA（HCOP）不可用，已改用 HomoloGene 映射:", type(exc).__name__, exc)
        net = dorothea_mouse_from_human_homologene()

    if len(net) < 100:
        print("警告: TF 网络行数较少:", len(net))

    dc.mt.ulm(pec_sym, net, raw=False)
    scores_adata = dc.pp.get_obsm(pec_sym, "score_ulm")
    tf_df = pd.DataFrame(scores_adata.X, index=pec_sym.obs_names, columns=scores_adata.var_names)
    summ = tf_df.join(pec_sym.obs[["injury_group"]]).groupby("injury_group").mean().T
    summ["delta_nephritis_minus_control"] = summ["nephritis"] - summ["control"]
    summ = summ.sort_values("delta_nephritis_minus_control", ascending=False)
    tf_out = TAB_DIR / "PEC_TF_DoRothEA_mean_by_condition_EXPLORATORY.csv"
    summ.to_csv(tf_out)
    print("探索性输出（非论文 SCENIC）:", tf_out)
    print(summ.head(20).to_string())
"""
)

code(
    r"""
# --- C. 轨迹：POD + PEC + TEC，PAGA + DPT ---
line = adata_cmp[adata_cmp.obs["cell_type_major"].isin(["POD", "PEC", "TEC"])].copy()
line_sym = collapse_ensembl_to_symbol(line)
sc.pp.highly_variable_genes(line_sym, n_top_genes=2000, flavor="seurat", subset=False)
line_sym.raw = line_sym.copy()
line_sym = line_sym[:, line_sym.var["highly_variable"]].copy()
sc.pp.scale(line_sym, max_value=10, zero_center=False)
sc.tl.pca(line_sym, svd_solver="arpack")
sc.pp.neighbors(line_sym, n_neighbors=20, n_pcs=30)
sc.tl.umap(line_sym)
sc.tl.paga(line_sym, groups="cell_type_major")
sc.pl.paga(line_sym, save="_PAGA_POD_PEC_TEC.pdf", show=False)
plt.show()

root_idx = np.flatnonzero((line_sym.obs["cell_type_major"] == "POD").values)
if root_idx.size > 0:
    line_sym.uns["iroot"] = int(root_idx[0])
    sc.tl.dpt(line_sym)
    sc.pl.umap(
        line_sym,
        color=["cell_type_major", "dpt_pseudotime"],
        legend_loc="on data",
        save="_dpt_lineage.pdf",
        show=False,
    )
    plt.show()
else:
    print("无 POD，跳过 DPT")
"""
)

md(
    r"""
### D. 细胞通讯

- **论文主路径**：**CellChat（R）** — 使用本笔记末尾导出的 `h5ad` 与 **`scripts/cellchat_GSE146912_PEC.R`**；对照/肾炎分层比较请分别子集后各建一个 CellChat 对象。
- **可选**：仅当 `EXPLORATORY_PYTHON_MULTIMETHOD = True` 时运行下方 **LIANA**（多方法共识，**不等价**于原文 CellChat 推断与图示）。
"""
)

code(
    r"""
glob_lr = None
pec_lr = None
_ord = None

if not EXPLORATORY_PYTHON_MULTIMETHOD:
    print("跳过 LIANA。请使用 scripts/cellchat_GSE146912_PEC.R。")
else:
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "liana"])
    import liana as li

    def _liana_order_col(df: pd.DataFrame):
        for c in ("magnitude_rank", "specificity_rank", "aggregate_rank"):
            if c in df.columns:
                return c
        cand = [c for c in df.columns if "rank" in c.lower()]
        if cand:
            return cand[0]
        raise ValueError("无法在 LIANA 结果中找到排序列: " + str(df.columns.tolist()))

    li.mt.rank_aggregate(
        adata_sym,
        groupby="cell_type_major",
        resource_name="mouseconsensus",
        expr_prop=0.1,
        use_raw=False,
        key_added="liana_global",
        verbose=True,
    )
    glob_lr = adata_sym.uns["liana_global"].copy()
    _ord = _liana_order_col(glob_lr)
    pec_lr = glob_lr[(glob_lr["source"] == "PEC") | (glob_lr["target"] == "PEC")].sort_values(
        _ord, ascending=True
    )
    pec_lr.to_csv(TAB_DIR / "LIANA_PEC_involved_global_EXPLORATORY.csv", index=False)
    print("PEC 相关 LR（全局，探索性）, 排序列=", _ord, pec_lr.shape)
    print(pec_lr.head(25).to_string())
"""
)

md(
    r"""
#### 气泡图 1（仅探索模式）：PEC 相关互作（LIANA `dotplot`）
"""
)

code(
    r"""
if EXPLORATORY_PYTHON_MULTIMETHOD and pec_lr is not None and _ord is not None:
    import liana.plotting as lipl

    def _liana_dotplot_cols(df: pd.DataFrame):
        for a, b in (
            ("magnitude_rank", "specificity_rank"),
            ("specificity_rank", "magnitude_rank"),
        ):
            if a in df.columns and b in df.columns:
                return a, b
        num = [
            c
            for c in df.columns
            if df[c].dtype.kind in "fiu"
            and c not in ("source", "target")
            and "ligand" not in c.lower()
            and "receptor" not in c.lower()
        ]
        if len(num) >= 2:
            return num[0], num[1]
        raise ValueError("找不到用于气泡图的两个数值列: " + str(df.columns.tolist()))

    _top_n_bubble = 45
    pec_top = pec_lr.nsmallest(_top_n_bubble, _ord).copy()
    _col_c, _col_s = _liana_dotplot_cols(pec_top)
    p = lipl.dotplot(
        liana_res=pec_top,
        colour=_col_c,
        size=_col_s,
        inverse_colour=False,
        inverse_size=False,
        top_n=None,
        figure_size=(14, 16),
        cmap="viridis",
        return_fig=True,
    )
    bubble_path = FIG_DIR / "LIANA_PEC_bubble_dotplot_global_EXPLORATORY.pdf"
    p.save(str(bubble_path), limitsize=False)
    print("气泡图已保存:", bubble_path)
    p
else:
    print("跳过 LIANA dotplot（探索模式关闭或无结果）")
"""
)

code(
    r"""
merged_pec = None

if EXPLORATORY_PYTHON_MULTIMETHOD:
    import liana as li

    def run_liana_subset(ad_sub, key_tag):
        ad_sub = ad_sub.copy()
        li.mt.rank_aggregate(
            ad_sub,
            groupby="cell_type_major",
            resource_name="mouseconsensus",
            expr_prop=0.1,
            use_raw=False,
            key_added=f"liana_{key_tag}",
            verbose=False,
        )
        return ad_sub.uns[f"liana_{key_tag}"].copy()

    def _liana_order_col2(df: pd.DataFrame):
        for c in ("magnitude_rank", "specificity_rank", "aggregate_rank"):
            if c in df.columns:
                return c
        cand = [c for c in df.columns if "rank" in c.lower()]
        if cand:
            return cand[0]
        raise ValueError("LIANA 列: " + str(df.columns.tolist()))

    ctrl = adata_sym[adata_sym.obs["injury_group"] == "control"].copy()
    neph = adata_sym[adata_sym.obs["injury_group"] == "nephritis"].copy()

    df_c = run_liana_subset(ctrl, "control")
    df_n = run_liana_subset(neph, "nephritis")
    rank_col = _liana_order_col2(df_c)
    if rank_col not in df_n.columns:
        rank_col = _liana_order_col2(df_n)

    base_keys = ["source", "target", "ligand", "receptor"]
    extra = [c for c in df_c.columns if c.endswith("_complex") and c in df_n.columns]
    merge_keys = [k for k in base_keys + extra if k in df_c.columns and k in df_n.columns]

    merged = df_c[merge_keys + [rank_col]].merge(
        df_n[merge_keys + [rank_col]],
        on=merge_keys,
        suffixes=("_ctrl", "_neph"),
        how="outer",
    )
    merged["delta_rank_ctrl_minus_neph"] = merged[f"{rank_col}_ctrl"] - merged[f"{rank_col}_neph"]
    merged_pec = merged[(merged["source"] == "PEC") | (merged["target"] == "PEC")].sort_values(
        "delta_rank_ctrl_minus_neph", ascending=False
    )
    merged_pec.to_csv(TAB_DIR / "LIANA_PEC_stratified_rankDiff_EXPLORATORY.csv", index=False)
    print("分层 PEC 通讯（探索性）, rank列=", rank_col, merged_pec.shape)
    print(merged_pec.head(30).to_string())
else:
    print("跳过分层 LIANA。")
"""
)

md(
    r"""
#### 气泡图 2：对照 vs 肾炎 — `delta_rank_ctrl_minus_neph`（越大表示肾炎中通讯越强）

横轴为 rank 差；点大小表示绝对变化。仅展示非缺失且差值较大的互作对。
"""
)

code(
    r"""
import matplotlib.pyplot as plt
import seaborn as sns

if not EXPLORATORY_PYTHON_MULTIMETHOD or merged_pec is None:
    print("跳过分层气泡图。")
elif merged_pec.dropna(subset=["delta_rank_ctrl_minus_neph"]).empty:
    print("无足够数据绘制分层气泡图")
else:
    plot_df = merged_pec.dropna(subset=["delta_rank_ctrl_minus_neph"]).copy()
    lig_col = "ligand_complex" if "ligand_complex" in plot_df.columns else "ligand"
    rec_col = "receptor_complex" if "receptor_complex" in plot_df.columns else "receptor"
    plot_df["LR_pair"] = (
        plot_df[lig_col].astype(str) + " → " + plot_df[rec_col].astype(str)
    )
    plot_df["cell_pair"] = plot_df["source"].astype(str) + " → " + plot_df["target"].astype(str)
    topk = plot_df.nlargest(40, "delta_rank_ctrl_minus_neph").copy()
    topk["PEC_side"] = np.where(topk["source"].astype(str) == "PEC", "PEC 为发送细胞", "PEC 为接收细胞")
    fig, ax = plt.subplots(figsize=(10, max(6, len(topk) * 0.22)))
    sns.scatterplot(
        data=topk,
        x="delta_rank_ctrl_minus_neph",
        y="LR_pair",
        size=np.abs(topk["delta_rank_ctrl_minus_neph"].values),
        hue="PEC_side",
        sizes=(80, 400),
        alpha=0.85,
        ax=ax,
    )
    ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("Δ rank (control − nephritis)；正值 → 肾炎中更强")
    ax.set_ylabel("Ligand → Receptor")
    ax.set_title("PEC 相关互作：肾炎 vs 对照")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    outp = FIG_DIR / "LIANA_PEC_bubble_stratified_delta_rank_EXPLORATORY.pdf"
    fig.savefig(outp, dpi=200, bbox_inches="tight")
    plt.close()
    print("分层气泡图:", outp)
"""
)

code(
    r"""
out_r = OUT_DIR / "GSE146912_for_CellChat_symbol_injury.h5ad"
adata_sym.write_h5ad(out_r, compression="gzip")
print("导出供 R CellChat:", out_r)
"""
)

md(
    r"""
### 论文主流程小结（勿与探索模式混淆）

- **SCENIC / pySCENIC**：见 **`docs/METHODS_PAPER_ALIGNMENT.md`** 与 **`scripts/pyscenic_mouse_GSE146912.sh`**（需 cisTarget 数据库）。
- **CellChat（R）**：见 **`scripts/cellchat_GSE146912_PEC.R`**；输入为本笔记导出的 `h5ad`。

**引用（节选）**：Kidney International 2023（GSE146912）；SCENIC（Aibar *et al.*, *Nat Methods* 2017）；CellChat（Jin *et al.*, *Nat Commun* 2021）；若使用探索模式中的 LIANA / decoupler，另引对应文献。
"""
)

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12.0"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
print("Wrote", NB)
