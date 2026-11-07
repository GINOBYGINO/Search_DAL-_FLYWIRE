"""讀取 Codex `data/download/` 表，統一 root_id 為 int64。"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DOWNLOAD = ROOT / "data" / "download"


def _as_root_id(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce").astype("Int64")


def _read_csv(name: str, **kwargs) -> pd.DataFrame:
    path = DOWNLOAD / name
    if not path.exists():
        raise FileNotFoundError(f"缺少資料檔：{path}")
    return pd.read_csv(path, **kwargs)


def load_connections(
    *,
    no_threshold: bool = False,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """連線表。預設 `connections_princeton.csv.gz`。"""
    name = (
        "connections_princeton_no_threshold.csv.gz"
        if no_threshold
        else "connections_princeton.csv.gz"
    )
    usecols = columns or [
        "pre_root_id",
        "post_root_id",
        "neuropil",
        "syn_count",
        "nt_type",
    ]
    df = _read_csv(name, usecols=usecols)
    df["pre_root_id"] = _as_root_id(df["pre_root_id"])
    df["post_root_id"] = _as_root_id(df["post_root_id"])
    df["syn_count"] = pd.to_numeric(df["syn_count"], errors="coerce").fillna(0).astype(int)
    return df.dropna(subset=["pre_root_id", "post_root_id"])


def load_neurons() -> pd.DataFrame:
    df = _read_csv("neurons.csv.gz")
    df["root_id"] = _as_root_id(df["root_id"])
    return df.dropna(subset=["root_id"])


def load_classification() -> pd.DataFrame:
    df = _read_csv("classification.csv.gz")
    df["root_id"] = _as_root_id(df["root_id"])
    return df.dropna(subset=["root_id"])


def load_cell_types() -> pd.DataFrame:
    df = _read_csv("consolidated_cell_types.csv.gz")
    df["root_id"] = _as_root_id(df["root_id"])
    return df.dropna(subset=["root_id"])


def load_coordinates() -> pd.DataFrame:
    df = _read_csv("coordinates.csv.gz")
    if "root_id" in df.columns:
        df["root_id"] = _as_root_id(df["root_id"])
        df = df.dropna(subset=["root_id"])
    return df


def roots_by_primary_type(cell_types: pd.DataFrame, type_name: str) -> list[int]:
    """精確比對 primary_type。"""
    m = cell_types["primary_type"].astype(str) == type_name
    return cell_types.loc[m, "root_id"].astype(int).tolist()


def roots_matching_type(
    cell_types: pd.DataFrame,
    pattern: str,
    *,
    regex: bool = False,
) -> pd.DataFrame:
    """primary_type 或 additional_type(s) 含 pattern。"""
    prim = cell_types["primary_type"].astype(str)
    add = (
        cell_types["additional_type(s)"].astype(str)
        if "additional_type(s)" in cell_types.columns
        else pd.Series("", index=cell_types.index)
    )
    if regex:
        mask = prim.str.contains(pattern, case=False, na=False, regex=True) | add.str.contains(
            pattern, case=False, na=False, regex=True
        )
    else:
        mask = prim.str.contains(pattern, case=False, na=False, regex=False) | add.str.contains(
            pattern, case=False, na=False, regex=False
        )
    return cell_types.loc[mask].copy()
