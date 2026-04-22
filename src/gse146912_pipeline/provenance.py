from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def hash_input_h5ad(path: Path, full_file: bool) -> dict[str, Any]:
    st = path.stat()
    out: dict[str, Any] = {
        "path": str(path.resolve()),
        "size_bytes": st.st_size,
        "mtime_utc": datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(),
    }
    if full_file:
        out["sha256"] = _sha256_file(path)
    else:
        out["sha256"] = None
        out["sha256_note"] = "full_file hashing disabled; set provenance.hash_input_full_file: true"
    return out


class RunManifest:
    def __init__(self, provenance_dir: Path):
        self.provenance_dir = Path(provenance_dir)
        self.provenance_dir.mkdir(parents=True, exist_ok=True)
        self.path = self.provenance_dir / "run_manifest.json"
        self.steps: list[dict[str, Any]] = []
        self.meta: dict[str, Any] = {
            "started_at_utc": datetime.now(timezone.utc).isoformat(),
        }

    def set_input_fingerprint(self, info: dict[str, Any]) -> None:
        self.meta["input_h5ad"] = info

    def set_config_snapshot(self, cfg: dict[str, Any]) -> None:
        # 不写入 _resolved 循环引用
        snap = {k: v for k, v in cfg.items() if not k.startswith("_")}
        self.meta["config_snapshot"] = snap

    def add_step(
        self,
        name: str,
        n_cells_before: int | None = None,
        n_cells_after: int | None = None,
        n_vars: int | None = None,
        parameters: dict[str, Any] | None = None,
        extra: dict[str, Any] | None = None,
    ) -> None:
        rec: dict[str, Any] = {
            "step": name,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "n_cells_before": n_cells_before,
            "n_cells_after": n_cells_after,
            "n_vars": n_vars,
            "parameters": parameters or {},
        }
        if extra:
            rec.update(extra)
        self.steps.append(rec)
        self._flush()

    def _flush(self) -> None:
        payload = {"meta": self.meta, "steps": self.steps}
        with self.path.open("w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)

    def write_step_file(self, name: str, data: dict[str, Any]) -> None:
        p = self.provenance_dir / f"step_{name}.json"
        with p.open("w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
