"""DAL 候選延伸分析（左腦／時鐘／巧合／MB 迴路／睡眠路徑）。

分塊讀取 no_threshold 連線，只保留與候選／MBON／時鐘／FB 相關的邊。
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.dal_filters import MBON_DAL_NEUROPILS, parse_position  # noqa: E402
from src.io import (  # noqa: E402
    DOWNLOAD,
    load_cell_types,
    load_classification,
    load_coordinates,
    load_neurons,
    roots_by_primary_type,
)

TODAY = date.today().strftime("%Y%m%d")
RES = ROOT / "results"
FIG = RES / "extension_figures"
FIG.mkdir(exist_ok=True)
CONN_PATH = DOWNLOAD / "connections_princeton_no_threshold.csv.gz"

SEED_RIGHT = {
    "CB3198": 720575940633320785,
    "CB2122": 720575940628347654,
    "CB0339": 720575940632389028,
}
CLOCK_GROUPS = {
    "s-LNv": ["s-LNv"],
    "l-LNv": ["l-LNv"],
    "LNd": ["LNd_CRY+", "LNd_CRY+_ITP+", "LNd_CRY-"],
    "DN1": ["DN1a", "DN1pA", "DN1pB", "DN1pC", "DN1pD", "DN1pE", "DN1-l"],
    "DN2": ["DN2"],
    "DN3": [
        "s-CPDN3A",
        "s-CPDN3B",
        "s-CPDN3C",
        "s-CPDN3D",
        "s-CPDN3E",
        "l-CPDN3",
        "APDN3",
    ],
}
CLOCK_EXTRA = {"5th-LNv": ["5th-LNv"]}
VFB_PREFIXES = ("FB1", "FB2")


def log_append(path: Path, lines: list[str]):
    with path.open("a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def meta_join() -> pd.DataFrame:
    ct = load_cell_types()[["root_id", "primary_type"]]
    cl = load_classification()[["root_id", "class", "super_class", "sub_class", "side"]]
    nt = load_neurons()[["root_id", "nt_type", "nt_type_score"]]
    return ct.merge(cl, on="root_id", how="left").merge(nt, on="root_id", how="left")


def first_coord(coords: pd.DataFrame) -> pd.DataFrame:
    if "position" in coords.columns:
        xyz = parse_position(coords["position"])
        out = pd.concat([coords[["root_id"]].reset_index(drop=True), xyz], axis=1)
    else:
        out = coords[[c for c in ["root_id", "x", "y", "z"] if c in coords.columns]].copy()
    return out.groupby("root_id", as_index=False).first()


def roots_by_types(ct: pd.DataFrame, types: list[str]) -> list[int]:
    return ct.loc[ct["primary_type"].astype(str).isin(types), "root_id"].astype(int).tolist()


def roots_by_prefix(ct: pd.DataFrame, prefixes: tuple[str, ...]) -> list[int]:
    s = ct["primary_type"].astype(str)
    m = False
    for p in prefixes:
        m = m | s.str.startswith(p, na=False)
    return ct.loc[m, "root_id"].astype(int).tolist()


def _cast(df: pd.DataFrame) -> pd.DataFrame:
    df = df.dropna(subset=["pre_root_id", "post_root_id"]).copy()
    df["pre_root_id"] = pd.to_numeric(df["pre_root_id"], errors="coerce").astype("int64")
    df["post_root_id"] = pd.to_numeric(df["post_root_id"], errors="coerce").astype("int64")
    df["syn_count"] = pd.to_numeric(df["syn_count"], errors="coerce").fillna(0).astype("int32")
    return df


def load_relevant_connections(keep_ids: set[int], clock_ids: set[int], mbon_ids: set[int], fb_ids: set[int]) -> pd.DataFrame:
    """兩趟分塊讀取。"""
    usecols = ["pre_root_id", "post_root_id", "neuropil", "syn_count"]
    parts1 = []
    print("pass1 chunked load...", flush=True)
    for i, chunk in enumerate(pd.read_csv(CONN_PATH, usecols=usecols, chunksize=2_000_000)):
        chunk = _cast(chunk)
        m = (
            chunk["pre_root_id"].isin(keep_ids)
            | chunk["post_root_id"].isin(keep_ids)
            | chunk["pre_root_id"].isin(mbon_ids)
            | chunk["post_root_id"].isin(mbon_ids)
            | chunk["pre_root_id"].isin(fb_ids)
            | chunk["post_root_id"].isin(fb_ids)
            | chunk["neuropil"].isin(("MB_CA_L", "MB_CA_R"))
            & (chunk["pre_root_id"].isin(keep_ids) | chunk["post_root_id"].isin(keep_ids))
        )
        # simplify calyx: already covered if keep_ids involved; keep all calyx touching keep_ids
        m = (
            chunk["pre_root_id"].isin(keep_ids)
            | chunk["post_root_id"].isin(keep_ids)
            | (
                (chunk["pre_root_id"].isin(mbon_ids) | chunk["post_root_id"].isin(mbon_ids))
                & (
                    chunk["pre_root_id"].isin(keep_ids | fb_ids | mbon_ids)
                    | chunk["post_root_id"].isin(keep_ids | fb_ids | mbon_ids)
                )
            )
            | (chunk["pre_root_id"].isin(fb_ids) & chunk["post_root_id"].isin(mbon_ids | keep_ids))
            | (chunk["pre_root_id"].isin(mbon_ids) & chunk["post_root_id"].isin(fb_ids | keep_ids))
        )
        hit = chunk.loc[m]
        if len(hit):
            parts1.append(hit)
        if (i + 1) % 5 == 0:
            print(f"  pass1 chunk {i+1}", flush=True)
    conn1 = pd.concat(parts1, ignore_index=True) if parts1 else pd.DataFrame(columns=usecols)
    print(f"pass1 rows={len(conn1)}", flush=True)

    # partners of keep_ids become mids for clock 2-hop
    partners = set(conn1.loc[conn1["post_root_id"].isin(keep_ids), "pre_root_id"].astype(int))
    partners |= set(conn1.loc[conn1["pre_root_id"].isin(keep_ids), "post_root_id"].astype(int))
    print(f"partners of candidates={len(partners)}", flush=True)

    parts2 = []
    print("pass2 clock↔partner edges...", flush=True)
    for i, chunk in enumerate(pd.read_csv(CONN_PATH, usecols=usecols, chunksize=2_000_000)):
        chunk = _cast(chunk)
        m = (
            (chunk["pre_root_id"].isin(clock_ids) & chunk["post_root_id"].isin(partners | keep_ids))
            | (chunk["pre_root_id"].isin(partners | keep_ids) & chunk["post_root_id"].isin(clock_ids))
        )
        hit = chunk.loc[m]
        if len(hit):
            parts2.append(hit)
        if (i + 1) % 5 == 0:
            print(f"  pass2 chunk {i+1}", flush=True)
    conn2 = pd.concat(parts2, ignore_index=True) if parts2 else pd.DataFrame(columns=usecols)
    print(f"pass2 rows={len(conn2)}", flush=True)

    conn = pd.concat([conn1, conn2], ignore_index=True)
    conn = (
        conn.groupby(["pre_root_id", "post_root_id", "neuropil"], as_index=False)["syn_count"].sum()
    )
    print(f"merged relevant rows={len(conn)} mem_mb={conn.memory_usage(deep=True).sum()/1e6:.1f}", flush=True)
    return conn


def agg_pairs(sub: pd.DataFrame) -> pd.DataFrame:
    if sub.empty:
        return pd.DataFrame(columns=["pre_root_id", "post_root_id", "syn_count"])
    return sub.groupby(["pre_root_id", "post_root_id"], as_index=False)["syn_count"].sum()


def syn_between(conn, pre_ids, post_ids, min_syn, neuropils=None) -> int:
    pre_ids = set(int(x) for x in pre_ids)
    post_ids = set(int(x) for x in post_ids)
    m = (
        conn["pre_root_id"].isin(pre_ids)
        & conn["post_root_id"].isin(post_ids)
        & (conn["syn_count"] >= min_syn)
    )
    if neuropils is not None:
        m = m & conn["neuropil"].isin(neuropils)
    return int(conn.loc[m, "syn_count"].sum()) if m.any() else 0


def syn_calyx(conn, rid, min_syn=1) -> int:
    m = (
        conn["neuropil"].isin(("MB_CA_L", "MB_CA_R"))
        & (conn["syn_count"] >= min_syn)
        & ((conn["pre_root_id"] == rid) | (conn["post_root_id"] == rid))
    )
    return int(conn.loc[m, "syn_count"].sum()) if m.any() else 0


def twohop_stats(conn, sources, targets, min_syn, detail_top=0):
    sources = [int(x) for x in sources]
    targets = [int(x) for x in targets]
    down = agg_pairs(
        conn.loc[conn["pre_root_id"].isin(sources) & (conn["syn_count"] >= min_syn)]
    )
    up = agg_pairs(
        conn.loc[conn["post_root_id"].isin(targets) & (conn["syn_count"] >= min_syn)]
    )
    if down.empty or up.empty:
        return 0, 0, pd.DataFrame()
    mids = set(down["post_root_id"].astype(int)) & set(up["pre_root_id"].astype(int))
    mids -= set(sources)
    mids -= set(targets)
    if not mids:
        return 0, 0, pd.DataFrame()
    d1 = down[down["post_root_id"].isin(mids)].rename(
        columns={"pre_root_id": "source", "post_root_id": "mid", "syn_count": "syn_1"}
    )
    d2 = up[up["pre_root_id"].isin(mids)].rename(
        columns={"pre_root_id": "mid", "post_root_id": "target", "syn_count": "syn_2"}
    )
    two = d1.merge(d2, on="mid", how="inner")
    if two.empty:
        return 0, 0, pd.DataFrame()
    two["syn_path"] = two["syn_1"] + two["syn_2"]
    two = two.sort_values("syn_path", ascending=False)
    return int(len(two)), int(two["syn_path"].iloc[0]), two.head(detail_top) if detail_top else two.head(0)


def savefig(path: Path):
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()


def pick_left_counterparts(meta, coords, conn, mbon_ids, kc_ids) -> pd.DataFrame:
    rows = []
    for typ, right_id in SEED_RIGHT.items():
        same = meta.loc[meta["primary_type"] == typ].copy()
        lefts = same.loc[same["side"].astype(str).str.lower() == "left"]
        rights = same.loc[same["side"].astype(str).str.lower() == "right"]
        r_meta = same.loc[same["root_id"] == right_id]
        r_ca = syn_calyx(conn, right_id, 1)
        r_mbon = syn_between(conn, mbon_ids, [right_id], 1, MBON_DAL_NEUROPILS)
        r_kc = syn_between(conn, [right_id], kc_ids, 1) + syn_between(conn, kc_ids, [right_id], 1)
        r_xyz = coords.loc[coords["root_id"] == right_id]
        rx = float(r_xyz["x"].iloc[0]) if len(r_xyz) and pd.notna(r_xyz["x"].iloc[0]) else np.nan
        ry = float(r_xyz["y"].iloc[0]) if len(r_xyz) and pd.notna(r_xyz["y"].iloc[0]) else np.nan
        rz = float(r_xyz["z"].iloc[0]) if len(r_xyz) and pd.notna(r_xyz["z"].iloc[0]) else np.nan
        cand_scores = []
        for lid in lefts["root_id"].astype(int).tolist():
            l_ca = syn_calyx(conn, lid, 1)
            l_mbon = syn_between(conn, mbon_ids, [lid], 1, MBON_DAL_NEUROPILS)
            l_kc = syn_between(conn, [lid], kc_ids, 1) + syn_between(conn, kc_ids, [lid], 1)
            l_xyz = coords.loc[coords["root_id"] == lid]
            lx = float(l_xyz["x"].iloc[0]) if len(l_xyz) and pd.notna(l_xyz["x"].iloc[0]) else np.nan
            ly = float(l_xyz["y"].iloc[0]) if len(l_xyz) and pd.notna(l_xyz["y"].iloc[0]) else np.nan
            lz = float(l_xyz["z"].iloc[0]) if len(l_xyz) and pd.notna(l_xyz["z"].iloc[0]) else np.nan
            mid = 512_000.0
            pos_pen = 0.0
            if np.isfinite(rx) and np.isfinite(lx):
                pos_pen += abs(abs(rx - mid) - abs(lx - mid)) / 1000.0
            if np.isfinite(ry) and np.isfinite(ly):
                pos_pen += abs(ry - ly) / 1000.0
            if np.isfinite(rz) and np.isfinite(lz):
                pos_pen += abs(rz - lz) / 1000.0
            nt = meta.loc[meta["root_id"] == lid, "nt_type"]
            cand_scores.append(
                {
                    "type": typ,
                    "role": "left_candidate",
                    "root_id": lid,
                    "side": "left",
                    "mb_ca_syn": l_ca,
                    "mbon_in_syn": l_mbon,
                    "kc_bi_syn": l_kc,
                    "x": lx,
                    "y": ly,
                    "z": lz,
                    "match_score_lower_better": abs(l_ca - r_ca)
                    + 5 * abs(l_mbon - r_mbon)
                    + abs(l_kc - r_kc)
                    + pos_pen,
                    "paired_right_root_id": right_id,
                    "nt_type": nt.iloc[0] if len(nt) else None,
                }
            )
        if cand_scores:
            best = sorted(cand_scores, key=lambda d: d["match_score_lower_better"])[0]
            for d in cand_scores:
                d["is_best_left_match"] = int(d["root_id"] == best["root_id"])
            rows.extend(cand_scores)
        rows.append(
            {
                "type": typ,
                "role": "right_seed",
                "root_id": right_id,
                "side": "right",
                "mb_ca_syn": r_ca,
                "mbon_in_syn": r_mbon,
                "kc_bi_syn": r_kc,
                "x": rx,
                "y": ry,
                "z": rz,
                "match_score_lower_better": 0.0,
                "paired_right_root_id": right_id,
                "is_best_left_match": 0,
                "nt_type": r_meta["nt_type"].iloc[0] if len(r_meta) else None,
                "n_type_left": int(len(lefts)),
                "n_type_right": int(len(rights)),
            }
        )
        for rid in rights["root_id"].astype(int).tolist():
            if rid == right_id:
                continue
            rxyz = coords.loc[coords["root_id"] == rid]
            nt = meta.loc[meta["root_id"] == rid, "nt_type"]
            rows.append(
                {
                    "type": typ,
                    "role": "right_other_same_type",
                    "root_id": rid,
                    "side": "right",
                    "mb_ca_syn": syn_calyx(conn, rid, 1),
                    "mbon_in_syn": syn_between(conn, mbon_ids, [rid], 1, MBON_DAL_NEUROPILS),
                    "kc_bi_syn": syn_between(conn, [rid], kc_ids, 1)
                    + syn_between(conn, kc_ids, [rid], 1),
                    "x": float(rxyz["x"].iloc[0]) if len(rxyz) else np.nan,
                    "y": float(rxyz["y"].iloc[0]) if len(rxyz) else np.nan,
                    "z": float(rxyz["z"].iloc[0]) if len(rxyz) else np.nan,
                    "match_score_lower_better": np.nan,
                    "paired_right_root_id": right_id,
                    "is_best_left_match": 0,
                    "nt_type": nt.iloc[0] if len(nt) else None,
                }
            )
    return pd.DataFrame(rows)


def main():
    progress = RES / f"_extension_progress_{TODAY}.txt"
    progress.write_text("start\n", encoding="utf-8")
    log_path = RES / f"log_{TODAY}_extension.md"
    log_path.write_text(
        f"# 延伸分析日誌 {TODAY}\n\n- 資料：v783 `connections_princeton_no_threshold`（分塊過濾）\n"
        "- 聲明：只寫結構候選；nt_type 為預測值\n\n",
        encoding="utf-8",
    )

    print("meta...", flush=True)
    meta = meta_join()
    ct = load_cell_types()
    coords = first_coord(load_coordinates())

    mbon14 = roots_by_primary_type(ct, "MBON14")
    mbon03 = roots_by_primary_type(ct, "MBON03")
    mbon09 = roots_by_primary_type(ct, "MBON09")
    mbon08_n = len(roots_by_primary_type(ct, "MBON08"))
    mbon_ids = set(mbon14 + mbon03 + mbon09)

    cb_ids = set(
        meta.loc[meta["primary_type"].isin(SEED_RIGHT.keys()), "root_id"].astype(int)
    )
    clock_ids = set()
    for types in list(CLOCK_GROUPS.values()) + list(CLOCK_EXTRA.values()):
        clock_ids |= set(roots_by_types(ct, types))
    vfb_ids = set(roots_by_prefix(ct, VFB_PREFIXES))
    exr_ids = set(roots_by_prefix(ct, ("ExR",)))
    fb_ids = vfb_ids | exr_ids

    kc = set(meta.loc[meta["class"].astype(str) == "Kenyon_Cell", "root_id"].astype(int))
    if not kc:
        kc = set(
            meta.loc[meta["primary_type"].astype(str).str.startswith("KC", na=False), "root_id"].astype(int)
        )
    pn = set(meta.loc[meta["class"].astype(str) == "ALPN", "root_id"].astype(int))
    pam = set(meta.loc[meta["primary_type"].astype(str).str.startswith("PAM", na=False), "root_id"].astype(int))
    ppl1 = set(
        meta.loc[meta["primary_type"].astype(str).str.match(r"PPL10[1-8]$", na=False), "root_id"].astype(int)
    )

    conn = load_relevant_connections(cb_ids, clock_ids, mbon_ids, fb_ids)
    progress.write_text("loaded conn\n", encoding="utf-8")

    # ----- 1 left -----
    print("1 left...", flush=True)
    left_tbl = pick_left_counterparts(meta, coords, conn, list(mbon_ids), list(kc))
    left_tbl.to_csv(RES / f"extension_left_hemisphere_{TODAY}.csv", index=False)
    focus = {}
    for typ, rid in SEED_RIGHT.items():
        focus[typ] = [(f"{typ}_R", rid)]
        sub = left_tbl[
            (left_tbl["type"] == typ)
            & (left_tbl["role"] == "left_candidate")
            & (left_tbl["is_best_left_match"] == 1)
        ]
        if len(sub):
            focus[typ].append((f"{typ}_L", int(sub["root_id"].iloc[0])))
    focus_df = pd.DataFrame(
        [{"type": t, "label": a, "root_id": b} for t, items in focus.items() for a, b in items]
    )
    focus_df.to_csv(RES / f"extension_focus_cells_{TODAY}.csv", index=False)
    focus_ids = [int(x) for x in focus_df.root_id]
    label_of = dict(zip(focus_df.root_id.astype(int), focus_df.label))

    plot_df = left_tbl[
        (left_tbl["role"] == "right_seed")
        | ((left_tbl["role"] == "left_candidate") & (left_tbl["is_best_left_match"] == 1))
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, col, title in zip(
        axes, ["mb_ca_syn", "mbon_in_syn"], ["MB_CA 接觸突觸和", "MBON14/03 輸入（CRE/SMP/SIP）"]
    ):
        for typ in SEED_RIGHT:
            sub = plot_df[plot_df["type"] == typ].sort_values("side")
            if sub.empty:
                continue
            xs = [0 if s == "left" else 1 for s in sub["side"]]
            ax.plot(xs, sub[col].tolist(), "o-", label=typ)
        ax.set_xticks([0, 1], ["left", "right"])
        ax.set_title(title)
        ax.legend(fontsize=8)
    savefig(FIG / f"left_right_compare_{TODAY}.png")
    lines = ["## 1. 左腦補齊", ""]
    for typ, items in focus.items():
        lines.append(f"- **{typ}**：" + ", ".join(f"{a}=`{b}`" for a, b in items))
        sub = plot_df[plot_df["type"] == typ][
            ["side", "root_id", "mb_ca_syn", "mbon_in_syn", "kc_bi_syn", "nt_type"]
        ]
        lines += ["```", sub.to_string(index=False), "```"]
    lines += [
        "",
        "白話：CB3198／CB2122 每側兩顆同型；左腦是最像右種子的候選，不是鑑定。",
        "",
    ]
    log_append(log_path, lines)
    progress.write_text("done1\n", encoding="utf-8")

    # ----- 2 clock -----
    print("2 clock...", flush=True)
    clock_rows, edge_rows = [], []
    for thr in (1, 5):
        for grp, types in {**CLOCK_GROUPS, **CLOCK_EXTRA}.items():
            cids = roots_by_types(ct, types)
            for rid in focus_ids:
                d_in = syn_between(conn, cids, [rid], thr)
                d_out = syn_between(conn, [rid], cids, thr)
                n2a, b2a, _ = twohop_stats(conn, cids, [rid], thr, 0)
                n2b, b2b, _ = twohop_stats(conn, [rid], cids, thr, 0)
                clock_rows.append(
                    {
                        "min_syn": thr,
                        "focus_label": label_of[rid],
                        "root_id": rid,
                        "clock_group": grp,
                        "n_clock_cells": len(cids),
                        "direct_clock_to_focus": d_in,
                        "direct_focus_to_clock": d_out,
                        "n_twohop_clock_to_focus": n2a,
                        "n_twohop_focus_to_clock": n2b,
                        "best_twohop_syn_path": max(b2a, b2b),
                        "has_direct": int(d_in + d_out > 0),
                        "has_within_2hop": int((d_in + d_out > 0) or n2a or n2b),
                    }
                )
                if d_in + d_out:
                    sub = agg_pairs(
                        conn.loc[
                            (
                                (conn["pre_root_id"].isin(cids) & (conn["post_root_id"] == rid))
                                | ((conn["pre_root_id"] == rid) & conn["post_root_id"].isin(cids))
                            )
                            & (conn["syn_count"] >= thr)
                        ]
                    )
                    for r in sub.itertuples(index=False):
                        edge_rows.append(
                            {
                                "min_syn": thr,
                                "focus_label": label_of[rid],
                                "clock_group": grp,
                                "direction": "clock→focus"
                                if int(r.pre_root_id) in set(cids)
                                else "focus→clock",
                                "pre": int(r.pre_root_id),
                                "post": int(r.post_root_id),
                                "syn_count": int(r.syn_count),
                            }
                        )
    clock_df = pd.DataFrame(clock_rows)
    clock_df.to_csv(RES / f"extension_clock_summary_{TODAY}.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(RES / f"extension_clock_direct_edges_{TODAY}.csv", index=False)
    heat = (
        clock_df[clock_df.min_syn == 1]
        .pivot_table(index="focus_label", columns="clock_group", values="has_within_2hop", aggfunc="max")
        .reindex(columns=list(CLOCK_GROUPS) + list(CLOCK_EXTRA))
    )
    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(heat.fillna(0).values, aspect="auto", cmap="Greens", vmin=0, vmax=1)
    ax.set_xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
    ax.set_yticks(range(heat.shape[0]), heat.index)
    ax.set_title("時鐘↔候選 2跳內（≥1）")
    fig.colorbar(im, ax=ax, fraction=0.046)
    savefig(FIG / f"clock_twohop_heatmap_{TODAY}.png")
    lines = [
        "## 2. 時鐘神經元連結（Yin 2025 動機）",
        "",
        "- Codex 已確認 s-LNv／l-LNv／LNd*／DN1*／DN2／DN3*；附 5th-LNv。",
    ]
    for thr in (1, 5):
        sub = clock_df[clock_df.min_syn == thr]
        lines.append(
            f"- ≥{thr}：2跳內組合={int(sub.has_within_2hop.sum())}；直接邊組合={int(sub.has_direct.sum())}"
        )
    dh = clock_df[clock_df.has_direct == 1]
    if len(dh):
        lines += [
            "```",
            dh[
                [
                    "min_syn",
                    "focus_label",
                    "clock_group",
                    "direct_clock_to_focus",
                    "direct_focus_to_clock",
                ]
            ].to_string(index=False),
            "```",
        ]
    else:
        lines.append("- **直接邊：空的**。")
    lines += ["白話：有邊≠功能時鐘輸入。", ""]
    log_append(log_path, lines)
    progress.write_text("done2\n", encoding="utf-8")

    # ----- 3 coincidence -----
    print("3 coincidence...", flush=True)
    coin_rows = []
    for thr in (1, 5):
        ups = agg_pairs(conn.loc[conn["post_root_id"].isin(focus_ids) & (conn["syn_count"] >= thr)])
        for rid in focus_ids:
            u = ups[ups.post_root_id == rid]
            total = int(u.syn_count.sum()) if len(u) else 0
            syn_pn = int(u.loc[u.pre_root_id.isin(pn), "syn_count"].sum()) if len(u) else 0
            syn_kc = int(u.loc[u.pre_root_id.isin(kc), "syn_count"].sum()) if len(u) else 0
            syn_pam = int(u.loc[u.pre_root_id.isin(pam), "syn_count"].sum()) if len(u) else 0
            syn_ppl = int(u.loc[u.pre_root_id.isin(ppl1), "syn_count"].sum()) if len(u) else 0
            odor, reinf = syn_pn + syn_kc, syn_pam + syn_ppl
            coin_rows.append(
                {
                    "min_syn": thr,
                    "focus_label": label_of[rid],
                    "root_id": rid,
                    "upstream_syn_total": total,
                    "syn_from_PN_ALPN": syn_pn,
                    "syn_from_KC": syn_kc,
                    "syn_odor_PN_plus_KC": odor,
                    "syn_from_PAM": syn_pam,
                    "syn_from_PPL1": syn_ppl,
                    "syn_reinf_PAM_plus_PPL1": reinf,
                    "frac_odor": odor / total if total else 0.0,
                    "frac_reinf": reinf / total if total else 0.0,
                    "has_both_odor_and_reinf": int(odor > 0 and reinf > 0),
                }
            )
    coin_df = pd.DataFrame(coin_rows)
    coin_df.to_csv(RES / f"extension_coincidence_upstream_{TODAY}.csv", index=False)
    sub1 = coin_df[coin_df.min_syn == 1]
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(sub1))
    ax.bar(x - 0.2, sub1.frac_odor, 0.4, label="氣味(PN+KC)")
    ax.bar(x + 0.2, sub1.frac_reinf, 0.4, label="獎懲(PAM+PPL1)")
    ax.set_xticks(x, sub1.focus_label, rotation=45, ha="right")
    ax.set_title("上游巧合占比（≥1）")
    ax.legend(fontsize=8)
    savefig(FIG / f"coincidence_fractions_{TODAY}.png")
    log_append(
        log_path,
        [
            "## 3. 巧合偵測（Xia 2005 動機）",
            "```",
            sub1[
                [
                    "focus_label",
                    "syn_odor_PN_plus_KC",
                    "syn_reinf_PAM_plus_PPL1",
                    "frac_odor",
                    "frac_reinf",
                    "has_both_odor_and_reinf",
                ]
            ].to_string(index=False),
            "```",
            "白話：結構巧合，不是 NMDA 實驗。",
            "",
        ],
    )
    progress.write_text("done3\n", encoding="utf-8")

    # ----- 4 MB -----
    print("4 MB...", flush=True)
    loop_rows, mbon_rows = [], []
    detail = RES / f"extension_mbon_twohop_detail_{TODAY}.csv"
    if detail.exists():
        detail.unlink()
    for thr in (1, 5):
        for rid in focus_ids:
            to_kc = syn_between(conn, [rid], kc, thr)
            from_kc = syn_between(conn, kc, [rid], thr)
            u = agg_pairs(conn.loc[(conn.post_root_id == rid) & (conn.syn_count >= thr)])
            d = agg_pairs(conn.loc[(conn.pre_root_id == rid) & (conn.syn_count >= thr)])
            kc_in = set(u.loc[u.pre_root_id.isin(kc), "pre_root_id"].astype(int))
            kc_out = set(d.loc[d.post_root_id.isin(kc), "post_root_id"].astype(int))
            loop_rows.append(
                {
                    "min_syn": thr,
                    "focus_label": label_of[rid],
                    "root_id": rid,
                    "syn_KC_to_focus": from_kc,
                    "syn_focus_to_KC": to_kc,
                    "n_KC_in_partners": len(kc_in),
                    "n_KC_out_partners": len(kc_out),
                    "n_KC_reciprocal_same_id": len(kc_in & kc_out),
                    "has_KC_in_and_out": int(bool(kc_in) and bool(kc_out)),
                }
            )
            for mname, mids in [
                ("MBON14_alpha3", mbon14),
                ("MBON03_beta2mp", mbon03),
                ("MBON09_gamma3beta1", mbon09),
            ]:
                n2, b2, two = twohop_stats(conn, mids, [rid], thr, 20)
                mbon_rows.append(
                    {
                        "min_syn": thr,
                        "focus_label": label_of[rid],
                        "root_id": rid,
                        "mbon_set": mname,
                        "direct_syn_any_neuropil": syn_between(conn, mids, [rid], thr),
                        "direct_syn_CRE_SMP_SIP": syn_between(
                            conn, mids, [rid], thr, MBON_DAL_NEUROPILS
                        ),
                        "n_twohop_paths": n2,
                        "best_twohop_syn_path": b2,
                    }
                )
                if len(two):
                    tmp = two.copy()
                    tmp["min_syn"] = thr
                    tmp["focus_label"] = label_of[rid]
                    tmp["mbon_set"] = mname
                    tmp.to_csv(detail, mode="a", header=not detail.exists(), index=False)
    loop_df = pd.DataFrame(loop_rows)
    mbon_df = pd.DataFrame(mbon_rows)
    loop_df.to_csv(RES / f"extension_mb_loop_{TODAY}.csv", index=False)
    mbon_df.to_csv(RES / f"extension_mbon_paths_{TODAY}.csv", index=False)
    sub = loop_df[loop_df.min_syn == 1]
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(sub))
    ax.bar(x - 0.2, sub.syn_KC_to_focus, 0.4, label="KC→候選")
    ax.bar(x + 0.2, sub.syn_focus_to_KC, 0.4, label="候選→KC")
    ax.set_xticks(x, sub.focus_label, rotation=45, ha="right")
    ax.set_title("KC↔候選（≥1）")
    ax.legend()
    savefig(FIG / f"mb_loop_kc_{TODAY}.png")
    log_append(
        log_path,
        [
            "## 4. MB 迴路與 MBON 路徑",
            f"- Wu2017 第三種 MBON＝γ3,γ3β′1（GRASP 陰性）。Codex MBON09={len(mbon09)}；MBON08={mbon08_n}。",
            "```",
            sub[
                [
                    "focus_label",
                    "syn_KC_to_focus",
                    "syn_focus_to_KC",
                    "n_KC_reciprocal_same_id",
                    "has_KC_in_and_out",
                ]
            ].to_string(index=False),
            "```",
            "```",
            mbon_df[mbon_df.min_syn == 1][
                [
                    "focus_label",
                    "mbon_set",
                    "direct_syn_CRE_SMP_SIP",
                    "direct_syn_any_neuropil",
                    "n_twohop_paths",
                    "best_twohop_syn_path",
                ]
            ].to_string(index=False),
            "```",
            "",
        ],
    )
    progress.write_text("done4\n", encoding="utf-8")

    # ----- 5 sleep -----
    print("5 sleep...", flush=True)
    sleep_rows = []
    for thr in (1, 5):
        for rid in focus_ids:
            for sname, sids in [("FB1_FB2_as_vFB", list(vfb_ids)), ("ExR_related", list(exr_ids))]:
                d_vm = syn_between(conn, sids, mbon03, thr)
                n2vm, _, _ = twohop_stats(conn, sids, mbon03, thr, 0)
                d_mc = syn_between(conn, mbon03, [rid], thr)
                n2mc, _, _ = twohop_stats(conn, mbon03, [rid], thr, 0)
                chain = int(d_vm > 0 and d_mc > 0)
                if chain:
                    mf = set(
                        agg_pairs(conn.loc[conn.pre_root_id.isin(sids) & (conn.syn_count >= thr)])
                        .loc[lambda d: d.post_root_id.isin(mbon03), "post_root_id"]
                        .astype(int)
                    )
                    mt = set(
                        agg_pairs(conn.loc[conn.pre_root_id.isin(mbon03) & (conn.syn_count >= thr)])
                        .loc[lambda d: d.post_root_id == rid, "pre_root_id"]
                        .astype(int)
                    )
                    nbr = len(mf & mt)
                else:
                    nbr = 0
                sleep_rows.append(
                    {
                        "min_syn": thr,
                        "focus_label": label_of[rid],
                        "root_id": rid,
                        "source_set": sname,
                        "n_source_cells": len(sids),
                        "syn_source_to_MBON03": d_vm,
                        "n_twohop_source_to_MBON03": n2vm,
                        "syn_MBON03_to_focus": d_mc,
                        "n_twohop_MBON03_to_focus": n2mc,
                        "chain_direct_source_MBON03_focus": chain,
                        "n_MBON03_bridging_both": nbr,
                        "syn_MBON03_to_source_reverse_lit": syn_between(conn, mbon03, sids, thr),
                    }
                )
    sleep_df = pd.DataFrame(sleep_rows)
    sleep_df.to_csv(RES / f"extension_sleep_vfb_{TODAY}.csv", index=False)
    sub = sleep_df[(sleep_df.min_syn == 1) & (sleep_df.source_set == "FB1_FB2_as_vFB")]
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(sub))
    ax.bar(x - 0.2, sub.syn_source_to_MBON03, 0.4, label="FB1/2→MBON03")
    ax.bar(x + 0.2, sub.syn_MBON03_to_focus, 0.4, label="MBON03→候選")
    ax.set_xticks(x, sub.focus_label, rotation=45, ha="right")
    ax.set_title("睡眠路徑（FB1/2 操作定義）")
    ax.legend(fontsize=8)
    savefig(FIG / f"sleep_vfb_chain_{TODAY}.png")
    log_append(
        log_path,
        [
            "## 5. 睡眠路徑",
            "- Codex 無字面 vFB；用 FB1*/FB2*，ExR* 對照。",
            "```",
            sleep_df[sleep_df.min_syn == 1][
                [
                    "focus_label",
                    "source_set",
                    "syn_source_to_MBON03",
                    "syn_MBON03_to_focus",
                    "chain_direct_source_MBON03_focus",
                    "n_MBON03_bridging_both",
                    "syn_MBON03_to_source_reverse_lit",
                ]
            ].to_string(index=False),
            "```",
            "",
            "## 總結：做了什麼／有結果／空的",
            "",
            "### 做了什麼",
            "1. 左腦配對 2. 時鐘直接/2跳 3. 上游巧合 4. KC↔候選與 MBON 路徑 5. FB→MBON03→候選",
            "",
            "### 有結果",
            "- 各 CSV／PNG 已寫入 `results/`；細節以表為準。",
            "",
            "### 空的或弱的",
            "- 時鐘直接邊若上節為空 → 空。",
            "- 字面 vFB 名稱：空（改 FB1/2）。",
            f"- MBON08：{mbon08_n} 顆。",
            "- ≥5 常使弱邊消失；NBLAST／濕實驗仍未做。",
            "",
            "### 文獻",
            "- Yin 2025 https://doi.org/10.3389/fnins.2025.1623251",
            "- Wu 2017 https://doi.org/10.1038/s41598-017-07600-2（第三 MBON＝γ3,γ3β′1）",
            "- Xia 2005：僅動機",
            "",
        ],
    )
    (RES / "log_extension_latest.md").write_text(log_path.read_text(encoding="utf-8"), encoding="utf-8")
    progress.write_text("all_done\n", encoding="utf-8")
    print("DONE", log_path, flush=True)


if __name__ == "__main__":
    main()
