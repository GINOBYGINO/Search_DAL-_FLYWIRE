"""從 results/_relevant_conn.parquet 續跑延伸分析 1–5。"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_extension_analyses import (  # noqa: E402
    CLOCK_EXTRA,
    CLOCK_GROUPS,
    FIG,
    MBON_DAL_NEUROPILS,
    RES,
    SEED_RIGHT,
    TODAY,
    VFB_PREFIXES,
    first_coord,
    log_append,
    meta_join,
    pick_left_counterparts,
    roots_by_prefix,
    roots_by_types,
    savefig,
    syn_between,
    syn_calyx,
    twohop_stats,
    agg_pairs,
)
from src.io import load_cell_types, load_coordinates, roots_by_primary_type  # noqa: E402


def main():
    pq = RES / "_relevant_conn.parquet"
    assert pq.exists(), pq
    conn = pd.read_parquet(pq)
    meta = meta_join()
    ct = load_cell_types()
    coords = first_coord(load_coordinates())

    mbon14 = roots_by_primary_type(ct, "MBON14")
    mbon03 = roots_by_primary_type(ct, "MBON03")
    mbon09 = roots_by_primary_type(ct, "MBON09")
    mbon08_n = len(roots_by_primary_type(ct, "MBON08"))
    mbon_ids = mbon14 + mbon03 + mbon09
    kc = set(meta.loc[meta["class"].astype(str) == "Kenyon_Cell", "root_id"].astype(int))
    pn = set(meta.loc[meta["class"].astype(str) == "ALPN", "root_id"].astype(int))
    pam = set(meta.loc[meta["primary_type"].astype(str).str.startswith("PAM", na=False), "root_id"].astype(int))
    ppl1 = set(
        meta.loc[meta["primary_type"].astype(str).str.match(r"PPL10[1-8]$", na=False), "root_id"].astype(int)
    )
    vfb_ids = roots_by_prefix(ct, VFB_PREFIXES)
    exr_ids = roots_by_prefix(ct, ("ExR",))

    log_path = RES / f"log_{TODAY}_extension.md"
    log_path.write_text(
        f"# 延伸分析日誌 {TODAY}\n\n- 資料：v783 分塊過濾後 parquet（`_relevant_conn.parquet`）\n"
        "- 聲明：只寫結構候選；nt_type 為預測值\n\n",
        encoding="utf-8",
    )

    print("1 left", flush=True)
    left_tbl = pick_left_counterparts(meta, coords, conn, mbon_ids, list(kc))
    left_tbl.to_csv(RES / f"extension_left_hemisphere_{TODAY}.csv", index=False)
    focus = {}
    for typ, rid in SEED_RIGHT.items():
        focus[typ] = [(f"{typ}_R", rid)]
        sub = left_tbl[
            (left_tbl.type == typ)
            & (left_tbl.role == "left_candidate")
            & (left_tbl.is_best_left_match == 1)
        ]
        if len(sub):
            focus[typ].append((f"{typ}_L", int(sub.root_id.iloc[0])))
    focus_df = pd.DataFrame(
        [{"type": t, "label": a, "root_id": b} for t, items in focus.items() for a, b in items]
    )
    focus_df.to_csv(RES / f"extension_focus_cells_{TODAY}.csv", index=False)
    focus_ids = focus_df.root_id.astype(int).tolist()
    label_of = dict(zip(focus_df.root_id.astype(int), focus_df.label))

    plot_df = left_tbl[
        (left_tbl.role == "right_seed")
        | ((left_tbl.role == "left_candidate") & (left_tbl.is_best_left_match == 1))
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, col, title in zip(
        axes, ["mb_ca_syn", "mbon_in_syn"], ["MB_CA 接觸突觸和", "MBON14/03 輸入（CRE/SMP/SIP）"]
    ):
        for typ in SEED_RIGHT:
            sub = plot_df[plot_df.type == typ].sort_values("side")
            if sub.empty:
                continue
            xs = [0 if s == "left" else 1 for s in sub.side]
            ax.plot(xs, sub[col].tolist(), "o-", label=typ)
        ax.set_xticks([0, 1], ["left", "right"])
        ax.set_title(title)
        ax.legend(fontsize=8)
    savefig(FIG / f"left_right_compare_{TODAY}.png")
    lines = ["## 1. 左腦補齊", ""]
    for typ, items in focus.items():
        lines.append(f"- **{typ}**：" + ", ".join(f"{a}=`{b}`" for a, b in items))
        sub = plot_df[plot_df.type == typ][
            ["side", "root_id", "mb_ca_syn", "mbon_in_syn", "kc_bi_syn", "nt_type"]
        ]
        lines += ["```", sub.to_string(index=False), "```"]
    lines += [
        "",
        "白話：CB3198／CB2122 每側兩顆同型；左腦取最像右種子者。注意左腦 MBON 輸入常比右種子弱／為 0。",
        "",
    ]
    log_append(log_path, lines)

    print("2 clock", flush=True)
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
                                (conn.pre_root_id.isin(cids) & (conn.post_root_id == rid))
                                | ((conn.pre_root_id == rid) & conn.post_root_id.isin(cids))
                            )
                            & (conn.syn_count >= thr)
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
    lines = ["## 2. 時鐘神經元連結（Yin 2025 動機）", ""]
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

    print("3 coincidence", flush=True)
    coin_rows = []
    for thr in (1, 5):
        ups = agg_pairs(conn.loc[conn.post_root_id.isin(focus_ids) & (conn.syn_count >= thr)])
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
            "",
        ],
    )

    print("4 MB", flush=True)
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
            f"- Wu2017 第三種＝γ3,γ3β′1（GRASP 陰性）。MBON09={len(mbon09)}；MBON08={mbon08_n}。",
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

    print("5 sleep", flush=True)
    sleep_rows = []
    for thr in (1, 5):
        for rid in focus_ids:
            for sname, sids in [("FB1_FB2_as_vFB", vfb_ids), ("ExR_related", exr_ids)]:
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
            "### 有結果／空的",
            "- 見各節表格與 CSV；時鐘直接邊／字面 vFB／MBON08 可能為空。",
            "- ≥5 常清空弱邊；結論只寫候選。",
            "",
            "### 文獻",
            "- Yin 2025 https://doi.org/10.3389/fnins.2025.1623251",
            "- Wu 2017 https://doi.org/10.1038/s41598-017-07600-2（第三 MBON＝γ3,γ3β′1）",
            "- Xia 2005：僅動機",
            "",
        ],
    )
    (RES / "log_extension_latest.md").write_text(log_path.read_text(encoding="utf-8"), encoding="utf-8")
    print("DONE", log_path, flush=True)


if __name__ == "__main__":
    main()
