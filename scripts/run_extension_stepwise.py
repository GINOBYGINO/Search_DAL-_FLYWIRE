"""延伸分析分步版：CSV＋PIL 圖＋log（避開本機 matplotlib savefig 崩潰）。"""
from __future__ import annotations

from datetime import date
from pathlib import Path
import gc
import sys

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

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
COLORS = [(31, 119, 180), (255, 127, 14), (44, 160, 44), (214, 39, 40), (148, 103, 189), (140, 86, 75)]


def font(size=14):
    for name in ("arial.ttf", "calibri.ttf", "C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/calibri.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            pass
    return ImageFont.load_default()


def _text(draw, xy, text, f, fill=(20, 20, 20)):
    draw.text(xy, text, font=f, fill=fill)


def save_paired_lines(path, series_by_label, title, ylab):
    """series_by_label: {label: (left_val, right_val)}"""
    W, H, pad = 900, 420, 60
    im = Image.new("RGB", (W, H), (255, 255, 255))
    dr = ImageDraw.Draw(im)
    f, fs = font(16), font(12)
    _text(dr, (pad, 15), title, f)
    vals = [v for pair in series_by_label.values() for v in pair]
    vmax = max(vals) * 1.15 if vals and max(vals) > 0 else 1.0
    x0, x1, y0, y1 = pad + 40, W - pad - 160, pad + 30, H - pad
    dr.line((x0, y1, x1, y1), fill=(0, 0, 0), width=2)
    dr.line((x0, y0, x0, y1), fill=(0, 0, 0), width=2)
    _text(dr, (x0 - 10, y1 + 8), "left", fs)
    _text(dr, (x1 - 20, y1 + 8), "right", fs)
    _text(dr, (10, y0), ylab, fs)
    for i, (lab, (lv, rv)) in enumerate(series_by_label.items()):
        col = COLORS[i % len(COLORS)]
        py = lambda v: y1 - int((v / vmax) * (y1 - y0))
        p1, p2 = (x0, py(lv)), (x1, py(rv))
        dr.line([p1, p2], fill=col, width=3)
        for p, v in ((p1, lv), (p2, rv)):
            dr.ellipse((p[0] - 5, p[1] - 5, p[0] + 5, p[1] + 5), fill=col)
            _text(dr, (p[0] + 6, p[1] - 12), str(int(v)), fs, col)
        _text(dr, (x1 + 20, y0 + 20 + i * 22), lab, fs, col)
    im.save(path)


def save_grouped_bars(path, labels, series, title):
    """series: [(name, values_list, color)]"""
    W, H, pad = 1000, 420, 60
    im = Image.new("RGB", (W, H), (255, 255, 255))
    dr = ImageDraw.Draw(im)
    f, fs = font(16), font(11)
    _text(dr, (pad, 12), title, f)
    n, g = len(labels), len(series)
    if n == 0:
        im.save(path)
        return
    allv = [v for _, vals, _ in series for v in vals]
    vmax = max(allv) * 1.15 if allv and max(allv) > 0 else 1.0
    x0, x1, y0, y1 = pad + 20, W - pad, pad + 40, H - pad - 40
    group_w = (x1 - x0) / n
    bar_w = max(6, int(group_w / (g + 1.5)))
    for i, lab in enumerate(labels):
        gx = x0 + i * group_w
        for j, (name, vals, col) in enumerate(series):
            v = vals[i] if i < len(vals) else 0
            bh = int((v / vmax) * (y1 - y0))
            bx = int(gx + (j + 0.5) * bar_w)
            dr.rectangle((bx, y1 - bh, bx + bar_w - 2, y1), fill=col)
        _text(dr, (int(gx + 4), y1 + 6), str(lab)[:12], fs)
    for j, (name, _, col) in enumerate(series):
        _text(dr, (x0 + j * 160, 36), name, fs, col)
    dr.line((x0, y1, x1, y1), fill=(0, 0, 0), width=1)
    im.save(path)


def save_heatmap(path, rows, cols, mat, title):
    W, H = 140 + 70 * len(cols), 80 + 36 * len(rows)
    im = Image.new("RGB", (max(W, 500), max(H, 240)), (255, 255, 255))
    dr = ImageDraw.Draw(im)
    f, fs = font(15), font(11)
    _text(dr, (20, 10), title, f)
    x0, y0, cw, rh = 110, 50, 60, 30
    for j, c in enumerate(cols):
        _text(dr, (x0 + j * cw + 5, y0 - 22), str(c), fs)
    for i, r in enumerate(rows):
        _text(dr, (10, y0 + i * rh + 8), str(r), fs)
        for j in range(len(cols)):
            v = int(mat[i][j]) if j < len(mat[i]) else 0
            col = (210, 245, 210) if v else (240, 240, 240)
            if v:
                col = (40, 140, 70)
            x, y = x0 + j * cw, y0 + i * rh
            dr.rectangle((x, y, x + cw - 4, y + rh - 4), fill=col, outline=(180, 180, 180))
            _text(dr, (x + 20, y + 6), str(v), fs, (255, 255, 255) if v else (80, 80, 80))
    im.save(path)


def agg(sub):
    if sub is None or len(sub) == 0:
        return pd.DataFrame(columns=["pre_root_id", "post_root_id", "syn_count"])
    return sub.groupby(["pre_root_id", "post_root_id"], as_index=False)["syn_count"].sum()


def syn_set(pairs, pre, post):
    pre_s, post_s = set(map(int, pre)), set(map(int, post))
    if not pre_s or not post_s or pairs.empty:
        return 0
    m = pairs.pre_root_id.isin(pre_s) & pairs.post_root_id.isin(post_s)
    return int(pairs.loc[m, "syn_count"].sum()) if m.any() else 0


def twohop(pairs, sources, targets, top=0):
    sources = set(map(int, sources))
    targets = set(map(int, targets))
    if pairs.empty or not sources or not targets:
        return 0, 0, pd.DataFrame()
    down = pairs[pairs.pre_root_id.isin(sources)]
    up = pairs[pairs.post_root_id.isin(targets)]
    if down.empty or up.empty:
        return 0, 0, pd.DataFrame()
    mids = set(down.post_root_id.astype(int)) & set(up.pre_root_id.astype(int))
    mids -= sources | targets
    if not mids:
        return 0, 0, pd.DataFrame()
    d1 = down[down.post_root_id.isin(mids)]
    d2 = up[up.pre_root_id.isin(mids)]
    if top <= 0:
        g1 = d1.groupby("post_root_id")["syn_count"].agg(["count", "max"])
        g2 = d2.groupby("pre_root_id")["syn_count"].agg(["count", "max"])
        common = g1.index.intersection(g2.index)
        if len(common) == 0:
            return 0, 0, pd.DataFrame()
        n_paths = int((g1.loc[common, "count"] * g2.loc[common, "count"]).sum())
        best = int((g1.loc[common, "max"] + g2.loc[common, "max"]).max())
        return n_paths, best, pd.DataFrame()
    d1 = d1.rename(columns={"pre_root_id": "source", "post_root_id": "mid", "syn_count": "syn_1"})
    d2 = d2.rename(columns={"pre_root_id": "mid", "post_root_id": "target", "syn_count": "syn_2"})
    two = d1.merge(d2, on="mid")
    if two.empty:
        return 0, 0, pd.DataFrame()
    two["syn_path"] = two.syn_1 + two.syn_2
    two = two.sort_values("syn_path", ascending=False)
    return len(two), int(two.syn_path.iloc[0]), two.head(top)


print("load meta", flush=True)
ct = load_cell_types()
meta = (
    ct[["root_id", "primary_type"]]
    .merge(load_classification()[["root_id", "class", "side"]], on="root_id", how="left")
    .merge(load_neurons()[["root_id", "nt_type"]], on="root_id", how="left")
)
coords_raw = load_coordinates()
xyz = parse_position(coords_raw["position"])
coords = pd.concat([coords_raw[["root_id"]].reset_index(drop=True), xyz], axis=1).groupby(
    "root_id", as_index=False
).first()
del coords_raw, xyz
gc.collect()

mbon14 = roots_by_primary_type(ct, "MBON14")
mbon03 = roots_by_primary_type(ct, "MBON03")
mbon09 = roots_by_primary_type(ct, "MBON09")
mbon = set(mbon14 + mbon03 + mbon09)
kc = set(meta.loc[meta["class"].astype(str) == "Kenyon_Cell", "root_id"].astype(int))
pn = set(meta.loc[meta["class"].astype(str) == "ALPN", "root_id"].astype(int))
pam = set(meta.loc[meta.primary_type.astype(str).str.startswith("PAM", na=False), "root_id"].astype(int))
ppl1 = set(meta.loc[meta.primary_type.astype(str).str.match(r"PPL10[1-8]$", na=False), "root_id"].astype(int))
cb_all = set(meta.loc[meta.primary_type.isin(SEED.keys()), "root_id"].astype(int))
vfb = ct.loc[ct.primary_type.astype(str).str.startswith(("FB1", "FB2"), na=False), "root_id"].astype(int).tolist()
exr = ct.loc[ct.primary_type.astype(str).str.startswith("ExR", na=False), "root_id"].astype(int).tolist()
print(f"sizes kc={len(kc)} pn={len(pn)} pam={len(pam)} ppl1={len(ppl1)} vfb={len(vfb)} exr={len(exr)}", flush=True)

print("load conn", flush=True)
conn = pd.read_parquet(RES / "_relevant_conn.parquet")
ca = conn[conn.neuropil.isin(("MB_CA_L", "MB_CA_R"))].copy()
mbon_np = conn[conn.pre_root_id.isin(mbon) & conn.neuropil.isin(MBON_DAL_NEUROPILS)].copy()
kc_mask = conn.pre_root_id.isin(cb_all) | conn.post_root_id.isin(cb_all)
kc_cand = conn.loc[kc_mask]
kc_e = kc_cand[
    (kc_cand.pre_root_id.isin(cb_all) & kc_cand.post_root_id.isin(kc))
    | (kc_cand.pre_root_id.isin(kc) & kc_cand.post_root_id.isin(cb_all))
].copy()
del kc_cand
gc.collect()

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
series_ca, series_mbon = {}, {}
for typ in SEED:
    sub = plot[plot.type == typ]
    if len(sub) < 2:
        continue
    lv = int(sub.loc[sub.side == "left", "mb_ca_syn"].iloc[0])
    rv = int(sub.loc[sub.side == "right", "mb_ca_syn"].iloc[0])
    series_ca[typ] = (lv, rv)
    lv = int(sub.loc[sub.side == "left", "mbon_in_syn"].iloc[0])
    rv = int(sub.loc[sub.side == "right", "mbon_in_syn"].iloc[0])
    series_mbon[typ] = (lv, rv)
save_paired_lines(FIG / f"left_right_compare_ca_{TODAY}.png", series_ca, "MB_CA contact L vs R", "syn")
save_paired_lines(FIG / f"left_right_compare_mbon_{TODAY}.png", series_mbon, "MBON input L vs R", "syn")
# also combined side-by-side via reusing ca as main deliverable name expected
save_paired_lines(FIG / f"left_right_compare_{TODAY}.png", series_ca, "MB_CA contact L vs R (see also mbon fig)", "syn")
print("1 done", focus_ids, flush=True)
del ca, mbon_np, kc_e
gc.collect()

print("agg pairs", flush=True)
pairs_raw = agg(conn)
print("pairs", pairs_raw.shape, flush=True)

clock_rows, edge_rows, coin_rows, loop_rows, mbon_rows, sleep_rows = [], [], [], [], [], []
detail = RES / f"extension_mbon_twohop_detail_{TODAY}.csv"
if detail.exists():
    detail.unlink()

for thr in (1, 5):
    print("thr", thr, flush=True)
    P = pairs_raw[pairs_raw.syn_count >= thr].copy()

    print("  clock", flush=True)
    for grp, types in CLOCK.items():
        cids = ct.loc[ct.primary_type.isin(types), "root_id"].astype(int).tolist()
        cset = set(cids)
        for rid in focus_ids:
            din = syn_set(P, cids, [rid])
            dout = syn_set(P, [rid], cids)
            n2a, b2a, _ = twohop(P, cids, [rid], 0)
            n2b, b2b, _ = twohop(P, [rid], cids, 0)
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
                    (P.pre_root_id.isin(cset) & (P.post_root_id == rid))
                    | ((P.pre_root_id == rid) & P.post_root_id.isin(cset))
                ]
                for r in sub.itertuples(index=False):
                    edge_rows.append(
                        dict(
                            min_syn=thr,
                            focus_label=label[rid],
                            clock_group=grp,
                            direction="clock->focus" if int(r.pre_root_id) in cset else "focus->clock",
                            pre=int(r.pre_root_id),
                            post=int(r.post_root_id),
                            syn_count=int(r.syn_count),
                        )
                    )

    print("  coincidence", flush=True)
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

    print("  mb loop / mbon", flush=True)
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

    print("  sleep", flush=True)
    for rid in focus_ids:
        for sname, sids in [("FB1_FB2_as_vFB", vfb), ("ExR_related", exr)]:
            dvm = syn_set(P, sids, mbon03)
            n2vm, _, _ = twohop(P, sids, mbon03, 0)
            dmc = syn_set(P, mbon03, [rid])
            n2mc, _, _ = twohop(P, mbon03, [rid], 0)
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
    del P
    gc.collect()
    pd.DataFrame(clock_rows).to_csv(RES / f"extension_clock_summary_{TODAY}.csv", index=False)
    pd.DataFrame(edge_rows).to_csv(RES / f"extension_clock_direct_edges_{TODAY}.csv", index=False)
    pd.DataFrame(coin_rows).to_csv(RES / f"extension_coincidence_upstream_{TODAY}.csv", index=False)
    pd.DataFrame(loop_rows).to_csv(RES / f"extension_mb_loop_{TODAY}.csv", index=False)
    pd.DataFrame(mbon_rows).to_csv(RES / f"extension_mbon_paths_{TODAY}.csv", index=False)
    pd.DataFrame(sleep_rows).to_csv(RES / f"extension_sleep_vfb_{TODAY}.csv", index=False)
    print("  thr saved", thr, flush=True)

cdf = pd.DataFrame(clock_rows)
coin = pd.DataFrame(coin_rows)
loop = pd.DataFrame(loop_rows)
mbon_df = pd.DataFrame(mbon_rows)
sleep = pd.DataFrame(sleep_rows)

print("figures", flush=True)
heat = (
    cdf[cdf.min_syn == 1]
    .pivot_table(index="focus_label", columns="clock_group", values="has_within_2hop", aggfunc="max")
    .reindex(columns=list(CLOCK))
)
save_heatmap(
    FIG / f"clock_twohop_heatmap_{TODAY}.png",
    list(heat.index),
    list(heat.columns),
    heat.fillna(0).astype(int).values.tolist(),
    "clock <-> focus within 2 hops (>=1)",
)

sub1 = coin[coin.min_syn == 1].copy()
save_grouped_bars(
    FIG / f"coincidence_fractions_{TODAY}.png",
    sub1.focus_label.tolist(),
    [
        ("odor(PN+KC)", sub1.frac_odor.tolist(), COLORS[0]),
        ("reinf(PAM+PPL1)", sub1.frac_reinf.tolist(), COLORS[1]),
    ],
    "upstream coincidence fractions (>=1)",
)

sub = loop[loop.min_syn == 1].copy()
save_grouped_bars(
    FIG / f"mb_loop_kc_{TODAY}.png",
    sub.focus_label.tolist(),
    [
        ("KC->focus", sub.syn_KC_to_focus.tolist(), COLORS[0]),
        ("focus->KC", sub.syn_focus_to_KC.tolist(), COLORS[1]),
    ],
    "KC <-> focus (>=1)",
)

sub = sleep[(sleep.min_syn == 1) & (sleep.source_set == "FB1_FB2_as_vFB")].copy()
save_grouped_bars(
    FIG / f"sleep_vfb_chain_{TODAY}.png",
    sub.focus_label.tolist(),
    [
        ("FB1/2->MBON03", sub.syn_source_to_MBON03.tolist(), COLORS[0]),
        ("MBON03->focus", sub.syn_MBON03_to_focus.tolist(), COLORS[1]),
    ],
    "sleep path (FB1/2 operational)",
)

print("log", flush=True)
left_csv = pd.read_csv(RES / f"extension_left_hemisphere_{TODAY}.csv")
lines = [
    f"# 延伸分析日誌 {TODAY}",
    "",
    "- 資料：v783；分塊過濾後 `_relevant_conn.parquet`",
    "- 只寫結構候選；`nt_type` 為預測值",
    "- 圖：本機 matplotlib savefig 崩潰，改以 PIL 輸出 PNG",
    "",
    "## 1. 左腦補齊",
    "- Codex：CB3198/CB2122 各 4 顆（左右各 2）；CB0339 左右各 1。",
    "- 配對：同型 left 中，以 3×|ΔMB_CA|＋5×|ΔMBON|＋|ΔKC|＋位置差 選最佳。",
]
for typ in SEED:
    items = focus_df[focus_df.type == typ]
    lines.append("- **" + typ + "**：" + ", ".join(f"{r.label}=`{int(r.root_id)}`" for r in items.itertuples()))
    subp = left_csv[
        (left_csv.type == typ)
        & ((left_csv.role == "right_seed") | ((left_csv.role == "left_candidate") & (left_csv.is_best_left_match == 1)))
    ][["side", "root_id", "mb_ca_syn", "mbon_in_syn", "kc_bi_syn", "nt_type"]]
    lines += ["```", subp.to_string(index=False), "```"]
lines += [
    "白話：左右對稱性看表。若左腦 MB_CA 或 MBON 明顯弱於右種子，不支持「每側一顆且同質」的 DAL 故事。",
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
]
n_direct = int(cdf.has_direct.sum())
n_2hop1 = int(cdf[cdf.min_syn == 1].has_within_2hop.sum())
n_both = int(coin[(coin.min_syn == 1) & (coin.has_both_odor_and_reinf == 1)].shape[0])
n_loop = int(loop[(loop.min_syn == 1) & (loop.has_KC_in_and_out == 1)].shape[0])
n_chain = int(sleep[(sleep.min_syn == 1) & (sleep.chain_direct_source_MBON03_focus == 1)].shape[0])
cb2122_l = focus_df[(focus_df.type == "CB2122") & (focus_df.label.str.endswith("_L"))]
weak_ca = ""
if len(cb2122_l):
    lid = int(cb2122_l.root_id.iloc[0])
    lrow = left_csv[left_csv.root_id == lid]
    if len(lrow) and int(lrow.mb_ca_syn.iloc[0]) == 0:
        weak_ca = "- CB2122 左腦最佳配對 **MB_CA=0**（相對右種子弱），僅結構配對候選。"
lines += [
    f"- 時鐘 ≥1 有 2 跳訊號組合數={n_2hop1}；直接邊組合數={n_direct}",
    f"- 巧合同時氣味+獎懲（≥1）細胞數={n_both}",
    f"- KC 雙向（≥1）細胞數={n_loop}",
    f"- FB1/2→MBON03→候選直接鏈（≥1）組合數={n_chain}",
    "",
    "### 空的或弱的",
]
if n_direct == 0:
    lines.append("- 時鐘**直接邊：空的**")
if weak_ca:
    lines.append(weak_ca)
lines += [
    "- 字面 `vFB` primary_type：空（改用 FB1*/FB2* 操作定義）。",
    "- MBON08：0；Wu 2017 第三種以 γ3,γ3β′1 → Codex **MBON09**。",
    "- ≥5 門檻常砍掉弱邊。",
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
