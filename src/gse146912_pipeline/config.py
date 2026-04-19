from __future__ import annotations

import os
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml


def project_root(config: dict[str, Any] | None = None) -> Path:
    env = os.environ.get("SCRNA_PROJECT_ROOT")
    if env:
        return Path(env).resolve()
    if config and config.get("project", {}).get("root"):
        return Path(config["project"]["root"]).resolve()
    return Path.cwd().resolve()


def load_config(path: Path | str) -> dict[str, Any]:
    path = Path(path)
    with path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return data


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(base)
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = merge_config(out[k], v)
        else:
            out[k] = deepcopy(v)
    return out


def resolve_paths(cfg: dict[str, Any], root: Path) -> dict[str, Any]:
    """将 paths.* 与 input 解析为绝对路径字符串（写入 provenance 用）。"""
    p = cfg.setdefault("_resolved", {})
    paths = cfg.get("paths", {})
    raw = root / paths.get("raw_h5ad", "")
    out_dir = root / paths.get("output_run_dir", "results/run")
    p["project_root"] = str(root)
    p["raw_h5ad"] = str(raw)
    p["output_run_dir"] = str(out_dir)
    p["figures_dir"] = str(out_dir / paths.get("figures_subdir", "figures"))
    p["tables_dir"] = str(out_dir / paths.get("tables_subdir", "tables"))
    p["provenance_dir"] = str(out_dir / paths.get("provenance_subdir", "provenance"))
    return cfg
