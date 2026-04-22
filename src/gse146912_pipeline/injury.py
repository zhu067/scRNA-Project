from __future__ import annotations

import pandas as pd


def assign_injury_group(
    sample: pd.Series,
    control_token: str = "control",
    nephritis_token: str = "nephritis",
) -> pd.Series:
    """根据 Sample 字符串粗分：control / nephritis / other（大小写不敏感）。"""

    def _one(s: str) -> str:
        sl = str(s).lower()
        if control_token.lower() in sl:
            return "control"
        if nephritis_token.lower() in sl:
            return "nephritis"
        return "other"

    return sample.map(_one)
