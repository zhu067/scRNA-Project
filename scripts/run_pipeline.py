#!/usr/bin/env python3
"""复现入口：从仓库根目录运行（或设置 SCRNA_PROJECT_ROOT）。"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gse146912_pipeline.pipeline_main import main_cli

if __name__ == "__main__":
    main_cli()
