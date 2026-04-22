from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _read_square_matrix_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def _path_graph_weights(
    delta: pd.DataFrame,
    positive_only: bool,
) -> tuple[list[str], dict[tuple[str, str], float]]:
    nodes = list(delta.index.astype(str))
    idx = set(nodes)
    w: dict[tuple[str, str], float] = {}
    for s in nodes:
        for t in nodes:
            if s == t:
                continue
            val = float(delta.loc[s, t])
            if positive_only:
                val = max(val, 0.0)
            else:
                val = abs(val)
            if val > 0:
                w[(s, t)] = val
    return nodes, w


def _top_k_paths(
    nodes: list[str],
    edge_w: dict[tuple[str, str], float],
    max_len: int,
    top_k: int,
) -> list[tuple[list[tuple[str, str]], float]]:
    """简单路径枚举（|V| 较小时可用）：路径得分 = 边权之和。"""
    adj: dict[str, list[tuple[str, float]]] = {n: [] for n in nodes}
    for (s, t), wt in edge_w.items():
        adj[s].append((t, wt))

    scored: list[tuple[list[tuple[str, str]], float]] = []

    def dfs(cur: str, visited: frozenset[str], edges: list[tuple[str, str]], score: float) -> None:
        if edges:
            scored.append((list(edges), score))
        if len(edges) >= max_len:
            return
        for nbr, wt in adj.get(cur, []):
            if nbr in visited:
                continue
            dfs(nbr, visited | {nbr}, edges + [(cur, nbr)], score + wt)

    for s in nodes:
        dfs(s, frozenset({s}), [], 0.0)

    scored.sort(key=lambda x: -x[1])
    out: list[tuple[list[tuple[str, str]], float]] = []
    seen_path: set[tuple[tuple[str, str], ...]] = set()
    for e_list, sc in scored:
        key = tuple(e_list)
        if key in seen_path:
            continue
        seen_path.add(key)
        out.append((e_list, sc))
        if len(out) >= top_k:
            break
    return out


def run_event_flow_prototype(
    cfg: dict[str, Any],
    tables_dir: Path,
    provenance: Any | None = None,
) -> dict[str, Any]:
    """以 CellChat nephritis−control 差异权重矩阵构造有向加权图，输出边/节点/路径表与简易图。"""
    ef = cfg.get("event_flow", {})
    rel = Path(ef.get("cellchat_diff_dir", "tables/PEC_injury_comm/cellchat_R_diff/diff"))
    wname = str(ef.get("weight_matrix_csv", "cellchat_diff_weight_matrix.csv"))
    top_k = int(ef.get("top_k_paths", 15))
    max_len = int(ef.get("max_path_length", 4))
    pos_only = bool(ef.get("path_use_positive_delta_only", True))

    mat_path = tables_dir / rel / wname
    if not mat_path.is_file():
        raise FileNotFoundError(
            f"未找到差异权重矩阵 {mat_path}。请先运行 R：scripts/cellchat_GSE146912_PEC.R"
        )

    delta_w = _read_square_matrix_csv(mat_path)
    nodes = list(delta_w.index.astype(str))
    for c in nodes:
        if c not in delta_w.columns:
            raise ValueError(f"差异矩阵行列不一致，缺少列 {c}")

    rows = []
    for s in nodes:
        for t in nodes:
            if s == t:
                continue
            dv = float(delta_w.loc[s, t])
            rows.append({"source": s, "target": t, "delta_weight": dv})
    edges_df = pd.DataFrame(rows)

    outflow = edges_df.groupby("source")["delta_weight"].sum().reindex(nodes, fill_value=0.0)
    inflow = edges_df.groupby("target")["delta_weight"].sum().reindex(nodes, fill_value=0.0)
    netflow = outflow - inflow
    nodes_df = pd.DataFrame(
        {
            "cell_type_major": nodes,
            "inflow": [float(inflow.get(n, 0.0)) for n in nodes],
            "outflow": [float(outflow.get(n, 0.0)) for n in nodes],
            "netflow": [float(netflow.get(n, 0.0)) for n in nodes],
        }
    )

    act_csv = ef.get("pec_node_activity_csv")
    if act_csv:
        p = Path(act_csv)
        if not p.is_file():
            p = tables_dir / act_csv
        if p.is_file():
            act = pd.read_csv(p)
            if "cell_type_major" in act.columns:
                key = (
                    "deg_direction"
                    if "deg_direction" in act.columns
                    else ("direction" if "direction" in act.columns else None)
                )
                if key:
                    sub = act[["cell_type_major", key]].drop_duplicates("cell_type_major")
                    sub = sub.rename(columns={key: "node_activity"})
                    nodes_df = nodes_df.merge(sub, on="cell_type_major", how="left")

    _, edge_w = _path_graph_weights(delta_w, positive_only=pos_only)
    paths = _top_k_paths(nodes, edge_w, max_len=max_len, top_k=top_k)
    path_rows = []
    for rank, (elist, sc) in enumerate(paths, start=1):
        seq = [elist[0][0]] + [e[1] for e in elist] if elist else []
        path_rows.append(
            {
                "rank": rank,
                "path": " -> ".join(seq),
                "score": sc,
                "hop": len(elist),
            }
        )
    paths_df = pd.DataFrame(path_rows)

    out_dir = tables_dir / "PEC_injury_comm" / "event_flow"
    out_dir.mkdir(parents=True, exist_ok=True)
    edges_df.to_csv(out_dir / "edges_delta.csv", index=False)
    nodes_df.to_csv(out_dir / "nodes_flow.csv", index=False)
    paths_df.to_csv(out_dir / "paths_topk.csv", index=False)

    # --- 可视化：净流量条形图 + 简易网络（边宽 ~ |delta|，颜色区分正负）---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axb = axes[0]
    order = nodes_df.sort_values("netflow", ascending=True)["cell_type_major"].tolist()
    y = np.arange(len(order))
    vals = nodes_df.set_index("cell_type_major").loc[order, "netflow"].values
    axb.barh(y, vals, color=["#c0392b" if v < 0 else "#27ae60" for v in vals])
    axb.set_yticks(y)
    axb.set_yticklabels(order, fontsize=8)
    axb.axvline(0, color="k", lw=0.6)
    axb.set_xlabel("netflow (outflow - inflow)")
    axb.set_title("Net signal: right = net sender, left = net receiver")

    axn = axes[1]
    n = len(nodes)
    ang = np.linspace(0, 2 * math.pi, n, endpoint=False)
    pos = {nodes[i]: (math.cos(ang[i]), math.sin(ang[i])) for i in range(n)}
    thr = float(np.percentile(np.abs(edges_df["delta_weight"].values), 85))
    thr = max(thr, 1e-9)
    for _, r in edges_df.iterrows():
        s, t, dw = r["source"], r["target"], float(r["delta_weight"])
        if abs(dw) < thr * 0.25:
            continue
        x0, y0 = pos[s]
        x1, y1 = pos[t]
        axn.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>",
                lw=min(3.0, 0.4 + 2.5 * abs(dw) / thr),
                color="#2980b9" if dw >= 0 else "#8e44ad",
                alpha=0.55,
                shrinkA=14,
                shrinkB=14,
                connectionstyle="arc3,rad=0.12",
            ),
        )
    for nn, (x, y) in pos.items():
        nf = float(nodes_df.set_index("cell_type_major").at[nn, "netflow"])
        axn.scatter([x], [y], s=120 + 40 * min(abs(nf), 3), c=["#e74c3c" if nf < 0 else "#2ecc71"], zorder=3)
        axn.text(x * 1.18, y * 1.18, nn, fontsize=7, ha="center", va="center")
    axn.set_aspect("equal")
    axn.axis("off")
    axn.set_title("Diff comm. network (weak edges filtered)")

    fig.tight_layout()
    fig_path = out_dir / "event_flow_overview.png"
    fig.savefig(fig_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "out_dir": str(out_dir),
        "edges_csv": str(out_dir / "edges_delta.csv"),
        "nodes_csv": str(out_dir / "nodes_flow.csv"),
        "paths_csv": str(out_dir / "paths_topk.csv"),
        "figure_png": str(fig_path),
        "weight_matrix_used": str(mat_path),
    }
    if provenance is not None:
        provenance.write_step_file("event_flow_prototype", summary)
    return summary
