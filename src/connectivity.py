"""依突觸門檻彙總輸入／輸出與 neuropil 分佈。"""

from __future__ import annotations

import pandas as pd

# FlyWire Codex：蕈狀體萼部
CALYX_NEUROPILS = ("MB_CA_L", "MB_CA_R")


def filter_edges(connections: pd.DataFrame, min_syn: int = 5) -> pd.DataFrame:
    return connections.loc[connections["syn_count"] >= min_syn].copy()


def aggregate_neuron_pairs(edges: pd.DataFrame) -> pd.DataFrame:
    """跨 neuropil 加總同一對 pre→post 的 syn_count。"""
    return (
        edges.groupby(["pre_root_id", "post_root_id"], as_index=False)["syn_count"]
        .sum()
        .sort_values("syn_count", ascending=False)
    )


def outputs_from(edges: pd.DataFrame, pre_ids) -> pd.DataFrame:
    pre_ids = set(int(x) for x in pre_ids)
    return edges.loc[edges["pre_root_id"].isin(pre_ids)].copy()


def inputs_to(edges: pd.DataFrame, post_ids) -> pd.DataFrame:
    post_ids = set(int(x) for x in post_ids)
    return edges.loc[edges["post_root_id"].isin(post_ids)].copy()


def partner_summary(
    edges: pd.DataFrame,
    *,
    direction: str = "downstream",
    seed_ids,
    meta: pd.DataFrame | None = None,
    top_n: int = 30,
) -> pd.DataFrame:
    """
    direction=downstream：seed 為 pre，列 post 夥伴。
    direction=upstream：seed 為 post，列 pre 夥伴。
    meta 需含 root_id，可含 primary_type / class 等。
    """
    seed_ids = [int(x) for x in seed_ids]
    if direction == "downstream":
        sub = outputs_from(edges, seed_ids)
        partner_col = "post_root_id"
    elif direction == "upstream":
        sub = inputs_to(edges, seed_ids)
        partner_col = "pre_root_id"
    else:
        raise ValueError("direction must be 'downstream' or 'upstream'")

    pairs = aggregate_neuron_pairs(sub)
    # rename partner
    if direction == "downstream":
        out = pairs.rename(columns={"post_root_id": "partner_id", "pre_root_id": "seed_id"})
    else:
        out = pairs.rename(columns={"pre_root_id": "partner_id", "post_root_id": "seed_id"})

    by_partner = (
        out.groupby("partner_id", as_index=False)["syn_count"]
        .sum()
        .sort_values("syn_count", ascending=False)
    )
    if meta is not None and "root_id" in meta.columns:
        by_partner = by_partner.merge(
            meta,
            left_on="partner_id",
            right_on="root_id",
            how="left",
        )
    return by_partner.head(top_n)


def neuropil_input_distribution(
    edges: pd.DataFrame,
    post_ids,
    *,
    top_n: int = 20,
) -> pd.DataFrame:
    """後突觸細胞的輸入突觸在各 neuropil 的總數與占比。"""
    sub = inputs_to(edges, post_ids)
    if sub.empty:
        return pd.DataFrame(columns=["neuropil", "syn_count", "fraction"])
    g = sub.groupby("neuropil", as_index=False)["syn_count"].sum()
    total = g["syn_count"].sum()
    g["fraction"] = g["syn_count"] / total if total else 0.0
    return g.sort_values("syn_count", ascending=False).head(top_n)


def calyx_input_fraction(edges: pd.DataFrame, post_ids) -> float:
    """輸入突觸落在 MB_CA_* 的比例。"""
    sub = inputs_to(edges, post_ids)
    if sub.empty:
        return float("nan")
    total = sub["syn_count"].sum()
    ca = sub.loc[sub["neuropil"].isin(CALYX_NEUROPILS), "syn_count"].sum()
    return float(ca) / float(total) if total else float("nan")


def synapse_weight_to_set(
    edges: pd.DataFrame,
    *,
    pre_ids=None,
    post_ids=None,
) -> pd.Series:
    """
    回傳 index=對方 root_id、value=突觸總數的 Series。
    指定 pre_ids：對每個 pre，加總到各 post（多 seed 則合併）。
    指定 post_ids：對每個 post，加總來自各 pre。
    """
    if pre_ids is not None:
        sub = outputs_from(edges, pre_ids)
        return sub.groupby("post_root_id")["syn_count"].sum()
    if post_ids is not None:
        sub = inputs_to(edges, post_ids)
        return sub.groupby("pre_root_id")["syn_count"].sum()
    raise ValueError("provide pre_ids or post_ids")
