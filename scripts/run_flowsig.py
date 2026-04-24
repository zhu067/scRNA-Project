#!/usr/bin/env python3
"""Run FlowSig from the repository root."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gse146912_pipeline.config import load_config, merge_config, project_root, resolve_paths
from gse146912_pipeline.flowsig_analysis import run_flowsig_stage


def main() -> None:
    ap = argparse.ArgumentParser(description="Run FlowSig using the injury CellChat exports")
    ap.add_argument("--config", type=Path, default=Path("config/default.yaml"))
    ap.add_argument("--local-config", type=Path, default=Path("config/local.yaml"))
    args = ap.parse_args()

    cfg = load_config(args.config)
    if args.local_config and args.local_config.is_file():
        cfg = merge_config(cfg, load_config(args.local_config))
    root = project_root(cfg)
    cfg = resolve_paths(cfg, root)

    res = run_flowsig_stage(cfg)
    print("FlowSig outputs:")
    for key, value in res.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main()
