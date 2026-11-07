"""延伸分析最終版：預先抽出相關邊，避免大集合反覆 isin 造成崩潰。"""
from __future__ import annotations

from datetime import date
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.dal_filters import MBON_DAL_NEUROPILS, parse_position
from src.io import (
    load_cell_types,
    load_classification,
    load_coordinates,
    load_neurons,
    roots_by_primary_type,
)

TODAY = date.today().strftime("%Y%m%d")
RES, FIG = ROOT / "results", ROOT / "results" / "extension_figures"
FIG.mkdir(exist_ok=True)
SEED = {
    "CB3198": 720575940633320785,
    "CB2122": 720575940628347654,
    "CB0339": 720575940632389028,
}
CLOCK = {
    "s-LNv": ["s-LNv"],
    "l-LNv": ["l-LNv"],
    "LNd": ["LNd_CRY+", "LNd_CRY+_ITP+", "LNd_CRY-"],
    "DN1": ["DN1a", "DN1pA", "DN1pB", "DN1pC", "DN1pD", "DN1pE", "DN1-l"],
    "DN2": ["DN2"],
    "DN3": ["s-CPDN3A", "s-CPDN3B", "s-CPDN3C", "s-CPDN3D", "s-CPDN3E", "l-CPDN3", "APDN3"],
    "5th-LNv": ["5th-LNv"],
}


def savefig(p):
    plt.tight_layout()
    plt.savefig(p, dpi=150, bbox_inches="tight")
    plt.close()


def agg(sub):
    if sub.empty:
        return pd.DataFrame(columns=["pre_root_id", "post_root_id", "syn_count"])
    return sub.groupby(["pre_root_id", "post_root_id"], as_index=False)["syn_count"].sum()


def syn_set(pairs, pre, post):
    pre, post = set(map(int, pre)), set(map(int, post))
    m = pairs.pre_root_id.isin(pre) & pairs.post_root_id.isin(post)
    return int(pairs.loc[m, "syn_count"].sum()) if m.any() else 0


def twohop(pairs, sources, targets, top=0):
    sources, targets = list(map(int, sources)), list(map(int, targets))
    down = pairs[pairs.pre_root_id.isin(sources)]
    up = pairs[pairs.post_root_id.isin(targets)]
    if down.empty or up.empty:
        return 0, 0, pd.DataFrame()
    mids = set(down.post_root_id.astype(int)) & set(up.pre_root_id.astype(int))
    mids -= set(sources) | set(targets)
    if not mids:
        return 0, 0, pd.DataFrame()
    d1 = down[down.post_root_id.isin(mids)].rename(
        columns={"pre_root_id": "source", "post_root_id": "mid", "syn_count": "syn_1"}
    )
    d2 = up[up.pre_root_id.isin(mids)].rename(
        columns={"pre_root_id": "mid", "post_root_id": "target", "syn_count": "syn_2"}
    )
    two = d1.merge(d2, on="mid")
    if two.empty:
        return 0, 0, pd.DataFrame()
    two["syn_path"] = two.syn_1 + two.syn_2
    two = two.sort_values("syn_path", ascending=False)
    return len(two), int(two.syn_path.iloc[0]), two.head(top) if top else two.head(0)


print("load", flush=True)
conn = pd.read_parquet(RES / "_relevant_conn.parquet")
meta = (
    load_cell_types()[["root_id", "primary_type"]]
    .merge(load_classification()[["root_id", "class", "side"]], on="root_id", how="left")
    .merge(load_neurons()[["root_id", "nt_type"]], on="root_id", how="left")
)
ct = load_cell_types()
coords_raw = load_coordinates()
xyz = parse_position(coords_raw["position"])
coords = pd.concat([coords_raw[["root_id"]].reset_index(drop=True), xyz], axis=1).groupby(
    "root_id", as_index=False
).first()

mbon14 = roots_by_primary_type(ct, "MBON14")
mbon03 = roots_by_primary_type(ct, "MBON03")
mbon09 = roots_by_primary_type(ct, "MBON09")
mbon = mbon14 + mbon03 + mbon09
kc = set(meta.loc[meta["class"].astype(str) == "Kenyon_Cell", "root_id"].astype(int))
pn = set(meta.loc[meta["class"].astype(str) == "ALPN", "root_id"].astype(int))
pam = set(meta.loc[meta.primary_type.astype(str).str.startswith("PAM", na=False), "root_id"].astype(int))
ppl1 = set(meta.loc[meta.primary_type.astype(str).str.match(r"PPL10[1-8]$", na=False), "root_id"].astype(int))
cb_all = set(meta.loc[meta.primary_type.isin(SEED.keys()), "root_id"].astype(int))
clock_all = set()
for ts in CLOCK.values():
    clock_all |= set(ct.loc[ct.primary_type.isin(ts), "root_id"].astype(int))
vfb = ct.loc[ct.primary_type.astype(str).str.startswith(("FB1", "FB2"), na=False), "root_id"].astype(int).tolist()
exr = ct.loc[ct.primary_type.astype(str).str.startswith("ExR", na=False), "root_id"].astype(int).tolist()

# 預先加總邊（門檻稍後再切）
print("precompute pairs", flush=True)
pairs_all = agg(conn)
ca = conn[conn.neuropil.isin(("MB_CA_L", "MB_CA_R"))]
mbon_np = conn[conn.pre_root_id.isin(mbon) & conn.neuropil.isin(MBON_DAL_NEUROPILS)]
kc_e = conn[
    (conn.pre_root_id.isin(cb_all) & conn.post_root_id.isin(kc))
    | (conn.pre_root_id.isin(kc) & conn.post_root_id.isin(cb_all))
]

def thr_pairs(df, thr):
    return agg(df[df.syn_count >= thr]) if len(df) else agg(df)

# ---------- 1 ----------
print("1 left", flush=True)
rows = []
for typ, rid in SEED.items():
    same = meta[meta.primary_type == typ]
    lefts = same[same.side.astype(str).str.lower() == "left"]
    rights = same[same.side.astype(str).str.lower() == "right"]
    r_ca = int(ca.loc[(ca.pre_root_id == rid) | (ca.post_root_id == rid), "syn_count"].sum())
    r_mbon = int(mbon_np.loc[mbon_np.post_root_id == rid, "syn_count"].sum())
    r_kc = int(kc_e.loc[(kc_e.pre_root_id == rid) | (kc_e.post_root_id == rid), "syn_count"].sum())
    rxyz = coords[coords.root_id == rid]
    rx, ry, rz = [float(rxyz[c].iloc[0]) if len(rxyz) else np.nan for c in "xyz"]
    nt = meta.loc[meta.root_id == rid, "nt_type"]
    left_scored = []
    for lid in lefts.root_id.astype(int).tolist():
        l_ca = int(ca.loc[(ca.pre_root_id == lid) | (ca.post_root_id == lid), "syn_count"].sum())
        l_mbon = int(mbon_np.loc[mbon_np.post_root_id == lid, "syn_count"].sum())
        l_kc = int(kc_e.loc[(kc_e.pre_root_id == lid) | (kc_e.post_root_id == lid), "syn_count"].sum())
        lxyz = coords[coords.root_id == lid]
        lx, ly, lz = [float(lxyz[c].iloc[0]) if len(lxyz) else np.nan for c in "xyz"]
        mid = 512000.0
        pen = 0.0
        if np.isfinite(rx) and np.isfinite(lx):
            pen += abs(abs(rx - mid) - abs(lx - mid)) / 1000
        if np.isfinite(ry) and np.isfinite(ly):
            pen += abs(ry - ly) / 1000
        if np.isfinite(rz) and np.isfinite(lz):
            pen += abs(rz - lz) / 1000
        # 萼部權重提高：DAL 關鍵特徵
        score = 3 * abs(l_ca - r_ca) + 5 * abs(l_mbon - r_mbon) + abs(l_kc - r_kc) + pen
        ntl = meta.loc[meta.root_id == lid, "nt_type"]
        left_scored.append(
            dict(
                type=typ,
                role="left_candidate",
                root_id=lid,
                side="left",
                mb_ca_syn=l_ca,
                mbon_in_syn=l_mbon,
                kc_bi_syn=l_kc,
                x=lx,
                y=ly,
                z=lz,
                match_score_lower_better=score,
                paired_right_root_id=rid,
                nt_type=ntl.iloc[0] if len(ntl) else None,
            )
        )
    if left_scored:
        best = min(left_scored, key=lambda d: d["match_score_lower_better"])
        for d in left_scored:
            d["is_best_left_match"] = int(d["root_id"] == best["root_id"])
        rows.extend(left_scored)
    rows.append(
        dict(
            type=typ,
            role="right_seed",
            root_id=rid,
            side="right",
            mb_ca_syn=r_ca,
            mbon_in_syn=r_mbon,
            kc_bi_syn=r_kc,
            x=rx,
            y=ry,
            z=rz,
            match_score_lower_better=0.0,
            paired_right_root_id=rid,
            is_best_left_match=0,
            nt_type=nt.iloc[0] if len(nt) else None,
            n_type_left=len(lefts),
            n_type_right=len(rights),
        )
    )
    for oid in rights.root_id.astype(int).tolist():
        if oid == rid:
            continue
        oxyz = coords[coords.root_id == oid]
        nto = meta.loc[meta.root_id == oid, "nt_type"]
        rows.append(
            dict(
                type=typ,
                role="right_other_same_type",
                root_id=oid,
                side="right",
                mb_ca_syn=int(ca.loc[(ca.pre_root_id == oid) | (ca.post_root_id == oid), "syn_count"].sum()),
                mbon_in_syn=int(mbon_np.loc[mbon_np.post_root_id == oid, "syn_count"].sum()),
                kc_bi_syn=int(kc_e.loc[(kc_e.pre_root_id == oid) | (kc_e.post_root_id == oid), "syn_count"].sum()),
                x=float(oxyz.x.iloc[0]) if len(oxyz) else np.nan,
                y=float(oxyz.y.iloc[0]) if len(oxyz) else np.nan,
                z=float(oxyz.z.iloc[0]) if len(oxyz) else np.nan,
                match_score_lower_better=np.nan,
                paired_right_root_id=rid,
                is_best_left_match=0,
                nt_type=nto.iloc[0] if len(nto) else None,
            )
        )

left = pd.DataFrame(rows)
left.to_csv(RES / f"extension_left_hemisphere_{TODAY}.csv", index=False)
focus = []
for typ, rid in SEED.items():
    focus.append(dict(type=typ, label=f"{typ}_R", root_id=rid))
    sub = left[(left.type == typ) & (left.role == "left_candidate") & (left.is_best_left_match == 1)]
    if len(sub):
        focus.append(dict(type=typ, label=f"{typ}_L", root_id=int(sub.root_id.iloc[0])))
focus_df = pd.DataFrame(focus)
focus_df.to_csv(RES / f"extension_focus_cells_{TODAY}.csv", index=False)
focus_ids = focus_df.root_id.astype(int).tolist()
label = dict(zip(focus_df.root_id.astype(int), focus_df.label))
plot = left[(left.role == "right_seed") | ((left.role == "left_candidate") & (left.is_best_left_match == 1))]
fig, axes = plt.subplots(1, 2, figsize=(9, 4))
for ax, col, title in zip(axes, ["mb_ca_syn", "mbon_in_syn"], ["MB_CA 接觸", "MBON 輸入"]):
    for typ in SEED:
        sub = plot[plot.type == typ].sort_values("side")
        if sub.empty:
            continue
        ax.plot([0 if s == "left" else 1 for s in sub.side], sub[col].tolist(), "o-", label=typ)
    ax.set_xticks([0, 1], ["left", "right"])
    ax.set_title(title)
    ax.legend(fontsize=8)
savefig(FIG / f"left_right_compare_{TODAY}.png")
print("1 done", flush=True)

# ---------- 2–5 for thr ----------
def pairs_at(thr):
    return thr_pairs(conn, thr)

clock_rows, edge_rows, coin_rows, loop_rows, mbon_rows, sleep_rows = [], [], [], [], [], []
detail = RES / f"extension_mbon_twohop_detail_{TODAY}.csv"
if detail.exists():
    detail.unlink()

for thr in (1, 5):
    print("thr", thr, flush=True)
    P = pairs_at(thr)
    # clock
    for grp, types in CLOCK.items():
        cids = ct.loc[ct.primary_type.isin(types), "root_id"].astype(int).tolist()
        for rid in focus_ids:
            din = syn_set(P, cids, [rid])
            dout = syn_set(P, [rid], cids)
            n2a, b2a, _ = twohop(P, cids, [rid])
            n2b, b2b, _ = twohop(P, [rid], cids)
            clock_rows.append(
                dict(
                    min_syn=thr,
                    focus_label=label[rid],
                    root_id=rid,
                    clock_group=grp,
                    n_clock_cells=len(cids),
                    direct_clock_to_focus=din,
                    direct_focus_to_clock=dout,
                    n_twohop_clock_to_focus=n2a,
                    n_twohop_focus_to_clock=n2b,
                    best_twohop_syn_path=max(b2a, b2b),
                    has_direct=int(din + dout > 0),
                    has_within_2hop=int((din + dout > 0) or n2a or n2b),
                )
            )
            if din + dout:
                sub = P[
                    (P.pre_root_id.isin(cids) & (P.post_root_id == rid))
                    | ((P.pre_root_id == rid) & P.post_root_id.isin(cids))
                ]
                for r in sub.itertuples(index=False):
                    edge_rows.append(
                        dict(
                            min_syn=thr,
                            focus_label=label[rid],
                            clock_group=grp,
                            direction="clock→focus" if int(r.pre_root_id) in set(cids) else "focus→clock",
                            pre=int(r.pre_root_id),
                            post=int(r.post_root_id),
                            syn_count=int(r.syn_count),
                        )
                    )
    # coincidence
    ups = P[P.post_root_id.isin(focus_ids)]
    for rid in focus_ids:
        u = ups[ups.post_root_id == rid]
        total = int(u.syn_count.sum()) if len(u) else 0
        sp = int(u.loc[u.pre_root_id.isin(pn), "syn_count"].sum()) if len(u) else 0
        sk = int(u.loc[u.pre_root_id.isin(kc), "syn_count"].sum()) if len(u) else 0
        sa = int(u.loc[u.pre_root_id.isin(pam), "syn_count"].sum()) if len(u) else 0
        sl = int(u.loc[u.pre_root_id.isin(ppl1), "syn_count"].sum()) if len(u) else 0
        odor, reinf = sp + sk, sa + sl
        coin_rows.append(
            dict(
                min_syn=thr,
                focus_label=label[rid],
                root_id=rid,
                upstream_syn_total=total,
                syn_from_PN_ALPN=sp,
                syn_from_KC=sk,
                syn_odor_PN_plus_KC=odor,
                syn_from_PAM=sa,
                syn_from_PPL1=sl,
                syn_reinf_PAM_plus_PPL1=reinf,
                frac_odor=odor / total if total else 0.0,
                frac_reinf=reinf / total if total else 0.0,
                has_both_odor_and_reinf=int(odor > 0 and reinf > 0),
            )
        )
    # MB loop + MBON
    for rid in focus_ids:
        u = P[P.post_root_id == rid]
        d = P[P.pre_root_id == rid]
        kin = set(u.loc[u.pre_root_id.isin(kc), "pre_root_id"].astype(int))
        kout = set(d.loc[d.post_root_id.isin(kc), "post_root_id"].astype(int))
        loop_rows.append(
            dict(
                min_syn=thr,
                focus_label=label[rid],
                root_id=rid,
                syn_KC_to_focus=int(u.loc[u.pre_root_id.isin(kc), "syn_count"].sum()) if len(u) else 0,
                syn_focus_to_KC=int(d.loc[d.post_root_id.isin(kc), "syn_count"].sum()) if len(d) else 0,
                n_KC_in_partners=len(kin),
                n_KC_out_partners=len(kout),
                n_KC_reciprocal_same_id=len(kin & kout),
                has_KC_in_and_out=int(bool(kin) and bool(kout)),
            )
        )
        for name, mids in [
            ("MBON14_alpha3", mbon14),
            ("MBON03_beta2mp", mbon03),
            ("MBON09_gamma3beta1", mbon09),
        ]:
            n2, b2, two = twohop(P, mids, [rid], 20)
            # neuropil-limited direct from raw conn
            d_np = int(
                conn.loc[
                    conn.pre_root_id.isin(mids)
                    & (conn.post_root_id == rid)
                    & conn.neuropil.isin(MBON_DAL_NEUROPILS)
                    & (conn.syn_count >= thr),
                    "syn_count",
                ].sum()
            )
            mbon_rows.append(
                dict(
                    min_syn=thr,
                    focus_label=label[rid],
                    root_id=rid,
                    mbon_set=name,
                    direct_syn_any_neuropil=syn_set(P, mids, [rid]),
                    direct_syn_CRE_SMP_SIP=d_np,
                    n_twohop_paths=n2,
                    best_twohop_syn_path=b2,
                )
            )
            if len(two):
                tmp = two.copy()
                tmp["min_syn"] = thr
                tmp["focus_label"] = label[rid]
                tmp["mbon_set"] = name
                tmp.to_csv(detail, mode="a", header=not detail.exists(), index=False)
    # sleep
    for rid in focus_ids:
        for sname, sids in [("FB1_FB2_as_vFB", vfb), ("ExR_related", exr)]:
            dvm = syn_set(P, sids, mbon03)
            n2vm, _, _ = twohop(P, sids, mbon03)
            dmc = syn_set(P, mbon03, [rid])
            n2mc, _, _ = twohop(P, mbon03, [rid])
            chain = int(dvm > 0 and dmc > 0)
            if chain:
                mf = set(P.loc[P.pre_root_id.isin(sids) & P.post_root_id.isin(mbon03), "post_root_id"].astype(int))
                mt = set(P.loc[P.pre_root_id.isin(mbon03) & (P.post_root_id == rid), "pre_root_id"].astype(int))
                nbr = len(mf & mt)
            else:
                nbr = 0
            sleep_rows.append(
                dict(
                    min_syn=thr,
                    focus_label=label[rid],
                    root_id=rid,
                    source_set=sname,
                    n_source_cells=len(sids),
                    syn_source_to_MBON03=dvm,
                    n_twohop_source_to_MBON03=n2vm,
                    syn_MBON03_to_focus=dmc,
                    n_twohop_MBON03_to_focus=n2mc,
                    chain_direct_source_MBON03_focus=chain,
                    n_MBON03_bridging_both=nbr,
                    syn_MBON03_to_source_reverse_lit=syn_set(P, mbon03, sids),
                )
            )

cdf = pd.DataFrame(clock_rows)
cdf.to_csv(RES / f"extension_clock_summary_{TODAY}.csv", index=False)
pd.DataFrame(edge_rows).to_csv(RES / f"extension_clock_direct_edges_{TODAY}.csv", index=False)
coin = pd.DataFrame(coin_rows)
coin.to_csv(RES / f"extension_coincidence_upstream_{TODAY}.csv", index=False)
loop = pd.DataFrame(loop_rows)
loop.to_csv(RES / f"extension_mb_loop_{TODAY}.csv", index=False)
mbon_df = pd.DataFrame(mbon_rows)
mbon_df.to_csv(RES / f"extension_mbon_paths_{TODAY}.csv", index=False)
sleep = pd.DataFrame(sleep_rows)
sleep.to_csv(RES / f"extension_sleep_vfb_{TODAY}.csv", index=False)

# figures
heat = (
    cdf[cdf.min_syn == 1]
    .pivot_table(index="focus_label", columns="clock_group", values="has_within_2hop", aggfunc="max")
    .reindex(columns=list(CLOCK))
)
fig, ax = plt.subplots(figsize=(8, 4))
im = ax.imshow(heat.fillna(0).values, aspect="auto", cmap="Greens", vmin=0, vmax=1)
ax.set_xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
ax.set_yticks(range(heat.shape[0]), heat.index)
ax.set_title("時鐘↔候選 2跳內（≥1）")
fig.colorbar(im, ax=ax, fraction=0.046)
savefig(FIG / f"clock_twohop_heatmap_{TODAY}.png")

sub1 = coin[coin.min_syn == 1]
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(sub1))
ax.bar(x - 0.2, sub1.frac_odor, 0.4, label="氣味(PN+KC)")
ax.bar(x + 0.2, sub1.frac_reinf, 0.4, label="獎懲(PAM+PPL1)")
ax.set_xticks(x, sub1.focus_label, rotation=45, ha="right")
ax.legend(fontsize=8)
ax.set_title("上游巧合占比（≥1）")
savefig(FIG / f"coincidence_fractions_{TODAY}.png")

sub = loop[loop.min_syn == 1]
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(sub))
ax.bar(x - 0.2, sub.syn_KC_to_focus, 0.4, label="KC→候選")
ax.bar(x + 0.2, sub.syn_focus_to_KC, 0.4, label="候選→KC")
ax.set_xticks(x, sub.focus_label, rotation=45, ha="right")
ax.legend()
ax.set_title("KC↔候選（≥1）")
savefig(FIG / f"mb_loop_kc_{TODAY}.png")

sub = sleep[(sleep.min_syn == 1) & (sleep.source_set == "FB1_FB2_as_vFB")]
fig, ax = plt.subplots(figsize=(8, 4))
x = np.arange(len(sub))
ax.bar(x - 0.2, sub.syn_source_to_MBON03, 0.4, label="FB1/2→MBON03")
ax.bar(x + 0.2, sub.syn_MBON03_to_focus, 0.4, label="MBON03→候選")
ax.set_xticks(x, sub.focus_label, rotation=45, ha="right")
ax.legend(fontsize=8)
ax.set_title("睡眠路徑（FB1/2 操作定義）")
savefig(FIG / f"sleep_vfb_chain_{TODAY}.png")

# log
lines = [
    f"# 延伸分析日誌 {TODAY}",
    "",
    "- 資料：v783；分塊過濾後 `_relevant_conn.parquet`",
    "- 只寫結構候選；`nt_type` 為預測值",
    "",
    "## 1. 左腦補齊",
    "- Codex：CB3198/CB2122 各 4 顆（左右各 2）；CB0339 左右各 1。",
    "- 配對：同型 left 中，以 3×|ΔMB_CA|＋5×|ΔMBON|＋|ΔKC|＋位置差 選最佳。",
]
for typ in SEED:
    items = focus_df[focus_df.type == typ]
    lines.append("- **" + typ + "**：" + ", ".join(f"{r.label}=`{int(r.root_id)}`" for r in items.itertuples()))
    sub = plot[plot.type == typ][["side", "root_id", "mb_ca_syn", "mbon_in_syn", "kc_bi_syn", "nt_type"]]
    lines += ["```", sub.to_string(index=False), "```"]
lines += [
    "白話：左右對稱性要看表。若左腦 MB_CA 或 MBON 明顯弱於右種子，不支持「每側一顆且同質」的 DAL 故事。",
    "",
    "## 2. 時鐘（Yin 2025 動機；Codex 名稱已核）",
    f"- ≥1：直接組合={int(cdf[cdf.min_syn==1].has_direct.sum())}；2跳內={int(cdf[cdf.min_syn==1].has_within_2hop.sum())}",
    f"- ≥5：直接組合={int(cdf[cdf.min_syn==5].has_direct.sum())}；2跳內={int(cdf[cdf.min_syn==5].has_within_2hop.sum())}",
]
dh = cdf[cdf.has_direct == 1]
lines += (
    [
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
    if len(dh)
    else ["- **直接邊：空的**"]
)
lines += [
    "白話：細胞內時鐘基因≠一定有很多中樞時鐘突觸。",
    "",
    "## 3. 巧合偵測（Xia 2005 動機）",
    "```",
    coin[coin.min_syn == 1][
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
    "## 4. MB 迴路／MBON（Wu 2017 第三種＝γ3,γ3β′1，GRASP 陰性；Codex MBON09 有、MBON08=0）",
    "```",
    loop[loop.min_syn == 1][
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
    "## 5. 睡眠路徑（無字面 vFB→用 FB1/FB2；ExR 對照）",
    "```",
    sleep[sleep.min_syn == 1][
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
    "1. 左腦配對與左右比較圖",
    "2. 時鐘六類＋5th-LNv 直接／2跳（≥1、≥5）",
    "3. 上游氣味 vs 獎懲占比",
    "4. KC↔候選；MBON14/03/09→候選",
    "5. FB1/2（與 ExR）→MBON03→候選",
    "",
    "### 有結果",
    "- 左腦、巧合、MB/MBON、睡眠鏈的 CSV／PNG 皆已產出（見 `results/extension_*`）。",
    "- 時鐘至少有 2 跳結構訊號的組合數見上節。",
    "",
    "### 空的或弱的",
    "- 時鐘**直接邊**若上節標空 → 空。",
    "- 字面 `vFB` primary_type：空。",
    "- MBON08：0。",
    "- 多數左腦種子的 MBON 輸入偏弱／0；≥5 常砍掉弱邊。",
    "- NBLAST／濕實驗／NT 驗證：未做。",
    "",
    "### 文獻",
    "- Yin 2025 https://doi.org/10.3389/fnins.2025.1623251",
    "- Wu 2017 https://doi.org/10.1038/s41598-017-07600-2",
    "- Xia 2005：僅動機",
    "",
]
text = "\n".join(lines)
(RES / f"log_{TODAY}_extension.md").write_text(text, encoding="utf-8")
(RES / "log_extension_latest.md").write_text(text, encoding="utf-8")
print("ALL DONE", flush=True)
