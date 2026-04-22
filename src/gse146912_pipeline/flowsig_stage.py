from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


def _read_square_matrix_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    if list(df.index) != list(df.columns):
        raise ValueError("CellChat 差异矩阵行列需一致且顺序对应。")
    return df


def _top_k_paths(g: nx.DiGraph, top_k: int, max_hops: int) -> pd.DataFrame:
    scored: list[tuple[str, float, int]] = []
    nodes = list(g.nodes())
    for s in nodes:
        for t in nodes:
            if s == t:
                continue
            for path in nx.all_simple_paths(g, source=s, target=t, cutoff=max_hops):
                if len(path) < 2:
                    continue
                score = 0.0
                for i in range(len(path) - 1):
                    score += float(g[path[i]][path[i + 1]].get("weight", 0.0))
                scored.append((" -> ".join(path), score, len(path) - 1))
    scored.sort(key=lambda x: -x[1])
    rows = []
    for i, (p, sc, hops) in enumerate(scored[:top_k], start=1):
        rows.append({"rank": i, "path": p, "score": sc, "hop": hops})
    return pd.DataFrame(rows)


def run_flowsig_stage(
    cfg: dict[str, Any],
    tables_dir: Path,
    provenance: Any | None = None,
) -> dict[str, Any]:
    """
    使用 FlowSig 构造通讯有向网络（最小接入版本）：
    - 输入: CellChat diff weight 矩阵 (nephritis - control)
    - 输出: edges / nodes / paths / overview 图
    """
    import flowsig

    fg = cfg.get("flowsig", {})
    rel = Path(fg.get("cellchat_diff_dir", "PEC_injury_comm/cellchat_R_diff/diff"))
    weight_csv = str(fg.get("weight_matrix_csv", "cellchat_diff_weight_matrix.csv"))
    top_k = int(fg.get("top_k_paths", 15))
    max_hops = int(fg.get("max_path_length", 4))
    positive_only = bool(fg.get("positive_delta_only", True))

    mat_path = tables_dir / rel / weight_csv
    if not mat_path.is_file():
        raise FileNotFoundError(f"未找到 FlowSig 输入矩阵: {mat_path}")

    delta = _read_square_matrix_csv(mat_path)
    nodes = list(delta.index)

    # FlowSig 的 CPDAG 邻接矩阵：这里使用差异权重阈值后的有向边（最小接入）
    adj = np.zeros((len(nodes), len(nodes)), dtype=float)
    edge_rows = []
    for i, s in enumerate(nodes):
        for j, t in enumerate(nodes):
            if i == j:
                continue
            d = float(delta.loc[s, t])
            if positive_only and d <= 0:
                continue
            w = d if positive_only else abs(d)
            if w > 0:
                adj[i, j] = w
                edge_rows.append({"source": s, "target": t, "delta_weight": d, "flowsig_weight": w})

    # 构造最小 AnnData + FlowSig network 结构，直接复用 FlowSig 网络工具函数
    a = ad.AnnData(X=np.zeros((1, len(nodes)), dtype=float))
    flow_var_info = pd.DataFrame(index=pd.Index(nodes, name="variable"))
    flow_var_info["Type"] = "module"
    a.uns["flowsig_network"] = {
        "flow_var_info": flow_var_info,
        "network": {"adjacency": adj},
    }

    flowsig.tl.apply_biological_flow(a, flowsig_network_key="flowsig_network")
    g = flowsig.tl.construct_intercellular_flow_network(
        a,
        flowsig_network_key="flowsig_network",
        adjacency_key="adjacency_validated",
    )

    out_dir = tables_dir / "PEC_injury_comm" / "flowsig"
    out_dir.mkdir(parents=True, exist_ok=True)

    edges_df = pd.DataFrame(edge_rows)
    if not edges_df.empty:
        edges_df.to_csv(out_dir / "flowsig_edges.csv", index=False)
    else:
        pd.DataFrame(columns=["source", "target", "delta_weight", "flowsig_weight"]).to_csv(
            out_dir / "flowsig_edges.csv", index=False
        )

    node_in = {n: 0.0 for n in nodes}
    node_out = {n: 0.0 for n in nodes}
    for u, v, d in g.edges(data=True):
        w = float(d.get("weight", 0.0))
        node_out[u] += w
        node_in[v] += w
    nodes_df = pd.DataFrame(
        {
            "cell_type_major": nodes,
            "inflow": [node_in[n] for n in nodes],
            "outflow": [node_out[n] for n in nodes],
            "netflow": [node_out[n] - node_in[n] for n in nodes],
        }
    )
    nodes_df.to_csv(out_dir / "flowsig_nodes.csv", index=False)

    paths_df = _top_k_paths(g, top_k=top_k, max_hops=max_hops)
    paths_df.to_csv(out_dir / "flowsig_paths_topk.csv", index=False)

    # 简图：环形网络
    fig, ax = plt.subplots(figsize=(6.5, 6.0))
    n = len(nodes)
    ang = np.linspace(0, 2 * math.pi, n, endpoint=False)
    pos = {nodes[i]: (math.cos(ang[i]), math.sin(ang[i])) for i in range(n)}
    for u, v, d in g.edges(data=True):
        w = float(d.get("weight", 0.0))
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>",
                lw=min(3.0, 0.6 + 2.0 * w),
                color="#2c7fb8",
                alpha=0.6,
                shrinkA=12,
                shrinkB=12,
                connectionstyle="arc3,rad=0.1",
            ),
        )
    for nn, (x, y) in pos.items():
        nf = float(nodes_df.set_index("cell_type_major").at[nn, "netflow"])
        ax.scatter([x], [y], s=160, c=["#d7301f" if nf < 0 else "#31a354"], zorder=3)
        ax.text(x * 1.18, y * 1.18, nn, fontsize=8, ha="center", va="center")
    ax.set_title("FlowSig network from CellChat delta")
    ax.axis("off")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig_path = out_dir / "flowsig_overview.png"
    fig.savefig(fig_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "flowsig_version": getattr(flowsig, "__version__", "unknown"),
        "weight_matrix_used": str(mat_path),
        "out_dir": str(out_dir),
        "edges_csv": str(out_dir / "flowsig_edges.csv"),
        "nodes_csv": str(out_dir / "flowsig_nodes.csv"),
        "paths_csv": str(out_dir / "flowsig_paths_topk.csv"),
        "figure_png": str(fig_path),
    }
    if provenance is not None:
        provenance.write_step_file("flowsig_stage", summary)
    return summary

