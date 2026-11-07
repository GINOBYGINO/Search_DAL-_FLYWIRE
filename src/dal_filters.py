"""DAL 結構候選篩選：各條件函式，門檻皆為參數。"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from .connectivity import CALYX_NEUROPILS, filter_edges, synapse_weight_to_set
from .io import (
    load_cell_types,
    load_classification,
    load_connections,
    load_coordinates,
    load_neurons,
    roots_by_primary_type,
)

# Wu 2017 GRASP：MBON–DAL 接觸在 SMP/SIP/CRE
MBON_DAL_NEUROPILS = ("CRE_L", "CRE_R", "SMP_L", "SMP_R", "SIP_L", "SIP_R")


@dataclass
class ScreenParams:
    min_syn_edge: int = 2
    min_syn_to_kc: int = 2
    min_syn_from_mbon: int = 2
    min_ca_syn: int = 2
    # 連線表：預設改用無門檻檔（含 1–4 突觸弱邊）
    no_threshold: bool = True
    # KC 集合：'kcabp' = 僅 KCab-p；'kcab' = 所有 primary_type 以 KCab 開頭
    kc_mode: str = "kcabp"
    mbon_types: tuple[str, ...] = ("MBON14", "MBON03")
    # 空 = 不限 neuropil；預設限制在 CRE/SMP/SIP
    mbon_neuropils: tuple[str, ...] = MBON_DAL_NEUROPILS
    nt_whitelist: tuple[str, ...] = ()
    exclude_classes: tuple[str, ...] = ("Kenyon_Cell", "MBON", "DAN", "ALPN", "LHCENT")
    exclude_super_classes: tuple[str, ...] = ("visual_projection",)
    exclude_primary_types: tuple[str, ...] = ("APL", "DPM", "LHMB1")
    require_side: str | None = None
    # 細胞本體粗框：預設關閉（避免誤殺）
    soma_y_min: float | None = None
    soma_y_max: float | None = None
    soma_z_min: float | None = None
    soma_z_max: float | None = None
    soma_x_min: float | None = None
    soma_x_max: float | None = None
    soma_midline_x: float | None = None
    soma_min_lateral_offset: float | None = None
    bidirectional_kc: bool = True
    calyx_binary_in_total: bool = True


def parse_position(series: pd.Series) -> pd.DataFrame:
    """將 `[x y z]` 字串拆成 x,y,z。"""
    cleaned = (
        series.astype(str)
        .str.replace("[", "", regex=False)
        .str.replace("]", "", regex=False)
        .str.replace(",", " ", regex=False)
        .str.strip()
        .str.split()
    )
    xyz = cleaned.apply(lambda xs: pd.Series(xs[:3] if isinstance(xs, list) else [np.nan] * 3))
    xyz.columns = ["x", "y", "z"]
    return xyz.apply(pd.to_numeric, errors="coerce")


def resolve_kc_ids(ct: pd.DataFrame, kc_mode: str) -> list[int]:
    mode = kc_mode.lower().strip()
    if mode in ("kcabp", "kcab-p", "p"):
        return roots_by_primary_type(ct, "KCab-p")
    if mode in ("kcab", "ab", "alpha_beta", "αβ"):
        m = ct["primary_type"].astype(str).str.startswith("KCab", na=False)
        return ct.loc[m, "root_id"].astype(int).tolist()
    raise ValueError(f"unknown kc_mode: {kc_mode}")


def score_output_to_set(
    edges: pd.DataFrame,
    target_ids: list[int],
    *,
    min_syn: int,
) -> pd.Series:
    """每個 pre 對 target 集合的輸出突觸總數；未達門檻為 0。"""
    w = synapse_weight_to_set(edges, post_ids=target_ids)
    s = w.copy()
    s[s < min_syn] = 0
    return s


def score_input_from_set(
    edges: pd.DataFrame,
    source_ids: list[int],
    *,
    min_syn: int,
    neuropils: tuple[str, ...] | None = None,
) -> pd.Series:
    """每個 post 來自 source 集合的輸入突觸總數；可限 neuropil。"""
    use = edges
    if neuropils:
        use = edges.loc[edges["neuropil"].isin(neuropils)]
    if use.empty:
        return pd.Series(dtype=float)
    w = synapse_weight_to_set(use, pre_ids=source_ids)
    s = w.copy()
    s[s < min_syn] = 0
    return s


def score_calyx_involvement(
    edges: pd.DataFrame,
    *,
    min_syn: int,
    neuropils: tuple[str, ...] = CALYX_NEUROPILS,
) -> pd.Series:
    """細胞作為 pre 或 post、落在萼部 neuropil 的突觸總數。"""
    sub = edges.loc[edges["neuropil"].isin(neuropils)]
    if sub.empty:
        return pd.Series(dtype=float)
    pre = sub.groupby("pre_root_id")["syn_count"].sum()
    post = sub.groupby("post_root_id")["syn_count"].sum()
    total = pre.add(post, fill_value=0)
    total[total < min_syn] = 0
    return total


def _apply_soma_mask(out: pd.DataFrame, params: ScreenParams) -> pd.Series:
    mask = pd.Series(True, index=out.index)
    has_xyz = out[["x", "y", "z"]].notna().all(axis=1)
    mask &= has_xyz
    if params.soma_x_min is not None:
        mask &= out["x"] >= params.soma_x_min
    if params.soma_x_max is not None:
        mask &= out["x"] <= params.soma_x_max
    if params.soma_y_min is not None:
        mask &= out["y"] >= params.soma_y_min
    if params.soma_y_max is not None:
        mask &= out["y"] <= params.soma_y_max
    if params.soma_z_min is not None:
        mask &= out["z"] >= params.soma_z_min
    if params.soma_z_max is not None:
        mask &= out["z"] <= params.soma_z_max
    if params.soma_midline_x is not None and params.soma_min_lateral_offset is not None:
        mask &= (out["x"] - params.soma_midline_x).abs() >= params.soma_min_lateral_offset
    return mask


def disable_soma(params: ScreenParams) -> ScreenParams:
    params.soma_y_min = None
    params.soma_y_max = None
    params.soma_z_min = None
    params.soma_z_max = None
    params.soma_x_min = None
    params.soma_x_max = None
    params.soma_midline_x = None
    params.soma_min_lateral_offset = None
    return params


def enable_soma_dal_box(params: ScreenParams) -> ScreenParams:
    """啟用先前粗估的背前外側框。"""
    params.soma_y_min = 110_000.0
    params.soma_y_max = 240_000.0
    params.soma_z_min = 100_000.0
    params.soma_midline_x = 512_000.0
    params.soma_min_lateral_offset = 60_000.0
    return params


def run_screen(
    params: ScreenParams | None = None,
    *,
    connections: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """回傳候選表（含各條件分數與總分）。可傳入已載入的 connections 以免重複讀檔。"""
    params = params or ScreenParams()
    conn = connections if connections is not None else load_connections(no_threshold=params.no_threshold)
    edges = filter_edges(conn, min_syn=params.min_syn_edge)
    ct = load_cell_types()
    cl = load_classification()
    neu = load_neurons()

    kc_ids = resolve_kc_ids(ct, params.kc_mode)
    mbon_ids: list[int] = []
    for t in params.mbon_types:
        mbon_ids.extend(roots_by_primary_type(ct, t))

    s_to_kc = score_output_to_set(edges, kc_ids, min_syn=params.min_syn_to_kc)
    s_from_kc = score_input_from_set(edges, kc_ids, min_syn=params.min_syn_to_kc)
    s_mbon = score_input_from_set(
        edges,
        mbon_ids,
        min_syn=params.min_syn_from_mbon,
        neuropils=params.mbon_neuropils or None,
    )
    s_ca = score_calyx_involvement(edges, min_syn=params.min_ca_syn)

    universe = sorted(
        set(s_to_kc[s_to_kc > 0].index.astype(int))
        | set(s_from_kc[s_from_kc > 0].index.astype(int))
        | set(s_mbon[s_mbon > 0].index.astype(int))
        | set(s_ca[s_ca > 0].index.astype(int))
    )
    out = pd.DataFrame({"root_id": universe})
    # 欄名維持 score_*_kcabp 以便與舊 CSV／日誌相容
    out["score_to_kcabp"] = out["root_id"].map(s_to_kc).fillna(0).astype(float)
    out["score_from_kcabp"] = out["root_id"].map(s_from_kc).fillna(0).astype(float)
    out["score_kcabp"] = out[["score_to_kcabp", "score_from_kcabp"]].max(axis=1)
    out["score_from_mbon"] = out["root_id"].map(s_mbon).fillna(0).astype(float)
    out["score_calyx"] = out["root_id"].map(s_ca).fillna(0).astype(float)
    out["n_kc_targets"] = len(kc_ids)

    if params.bidirectional_kc:
        out["pass_kcabp"] = (out["score_kcabp"] >= params.min_syn_to_kc).astype(int)
        out["pass_to_kcabp"] = (out["score_to_kcabp"] >= params.min_syn_to_kc).astype(int)
        out["pass_from_kcabp"] = (out["score_from_kcabp"] >= params.min_syn_to_kc).astype(int)
    else:
        out["pass_kcabp"] = (out["score_to_kcabp"] >= params.min_syn_to_kc).astype(int)
        out["pass_to_kcabp"] = out["pass_kcabp"]
        out["pass_from_kcabp"] = 0

    out["pass_mbon"] = (out["score_from_mbon"] >= params.min_syn_from_mbon).astype(int)
    out["pass_calyx"] = (out["score_calyx"] >= params.min_ca_syn).astype(int)

    out = out.merge(ct, on="root_id", how="left")
    out = out.merge(
        cl[["root_id", "flow", "super_class", "class", "sub_class", "side"]],
        on="root_id",
        how="left",
    )
    out = out.merge(neu[["root_id", "nt_type", "nt_type_score"]], on="root_id", how="left")

    if params.exclude_classes:
        out = out.loc[~out["class"].astype(str).isin(params.exclude_classes)]
    if params.exclude_super_classes:
        out = out.loc[~out["super_class"].astype(str).isin(params.exclude_super_classes)]
    if params.exclude_primary_types:
        out = out.loc[~out["primary_type"].astype(str).isin(params.exclude_primary_types)]

    if params.nt_whitelist:
        out["pass_nt"] = out["nt_type"].astype(str).isin(params.nt_whitelist).astype(int)
    else:
        out["pass_nt"] = 1

    if params.require_side:
        out = out.loc[out["side"].astype(str).str.lower() == params.require_side.lower()]

    need_soma = any(
        v is not None
        for v in (
            params.soma_x_min,
            params.soma_x_max,
            params.soma_y_min,
            params.soma_y_max,
            params.soma_z_min,
            params.soma_z_max,
            params.soma_midline_x if params.soma_min_lateral_offset else None,
        )
    )
    if need_soma:
        coords = load_coordinates()
        xyz = parse_position(coords["position"])
        coords = pd.concat([coords[["root_id"]], xyz], axis=1)
        coords = coords.dropna(subset=["x", "y", "z"]).drop_duplicates(
            subset=["root_id"], keep="first"
        )
        out = out.merge(coords, on="root_id", how="left")
        mask = _apply_soma_mask(out, params)
        out["pass_soma"] = mask.astype(int)
        out = out.loc[out["pass_soma"] == 1]
    else:
        out["pass_soma"] = 1

    out["score_structural"] = out["pass_kcabp"] + out["pass_mbon"] + out["pass_calyx"]
    out["score_pass_sum"] = out["score_structural"] + (
        out["pass_nt"] if params.nt_whitelist else 0
    )
    for col, new in [
        ("score_kcabp", "norm_kcabp"),
        ("score_from_mbon", "norm_mbon"),
    ]:
        mx = out[col].max()
        out[new] = out[col] / mx if mx and mx > 0 else 0.0
    if params.calyx_binary_in_total:
        out["norm_calyx"] = out["pass_calyx"].astype(float)
    else:
        mx = out["score_calyx"].max()
        out["norm_calyx"] = out["score_calyx"] / mx if mx and mx > 0 else 0.0
    out["score_total"] = (
        out["score_structural"] + out["norm_kcabp"] + out["norm_mbon"] + out["norm_calyx"]
    )

    kc_label = "KCab-p" if params.kc_mode.lower().startswith("kcabp") else "KCab*"

    def _reason(r: pd.Series) -> str:
        parts: list[str] = []
        if r["pass_to_kcabp"]:
            parts.append(f"to_{kc_label}")
        if r["pass_from_kcabp"]:
            parts.append(f"from_{kc_label}")
        if r["pass_mbon"]:
            parts.append("from_MBON")
        if r["pass_calyx"]:
            parts.append("calyx")
        if r["pass_nt"] and params.nt_whitelist:
            parts.append("nt")
        return ",".join(parts)

    out["reason"] = out.apply(_reason, axis=1)

    out = out.sort_values(
        ["score_pass_sum", "score_total", "score_kcabp", "score_from_mbon"],
        ascending=False,
    ).reset_index(drop=True)
    out.attrs["params"] = asdict(params)
    return out


def params_to_markdown(params: ScreenParams) -> str:
    d = asdict(params)
    lines = ["| 參數 | 值 |", "|---|---|"]
    for k, v in d.items():
        lines.append(f"| `{k}` | {v} |")
    return "\n".join(lines)
