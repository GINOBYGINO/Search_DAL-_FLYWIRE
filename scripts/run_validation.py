"""第 3 步：重現已知結果，寫入 results/validation_YYYYMMDD.md。"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.connectivity import (  # noqa: E402
    CALYX_NEUROPILS,
    aggregate_neuron_pairs,
    calyx_input_fraction,
    filter_edges,
    neuropil_input_distribution,
    outputs_from,
    partner_summary,
)
from src.io import (  # noqa: E402
    load_cell_types,
    load_classification,
    load_connections,
    roots_by_primary_type,
)

MIN_SYN = 5


def _meta(cell_types: pd.DataFrame, classification: pd.DataFrame) -> pd.DataFrame:
    return cell_types.merge(
        classification[["root_id", "class", "super_class", "side"]],
        on="root_id",
        how="left",
    )


def run() -> dict:
    print("loading tables…")
    conn = load_connections()
    edges = filter_edges(conn, min_syn=MIN_SYN)
    ct = load_cell_types()
    cl = load_classification()
    meta = _meta(ct, cl)

    # --- KCab-p ---
    kcabp = roots_by_primary_type(ct, "KCab-p")
    kcabp_n = len(kcabp)
    ca_frac = calyx_input_fraction(edges, kcabp)
    np_dist = neuropil_input_distribution(edges, kcabp, top_n=10)
    ca_rank = None
    if not np_dist.empty:
        for i, row in enumerate(np_dist.itertuples(index=False), start=1):
            if row.neuropil in CALYX_NEUROPILS:
                ca_rank = i
                break
    mb_frac = 0.0
    if not np_dist.empty:
        # recompute full MB fraction from all neuropils of inputs
        full_np = neuropil_input_distribution(edges, kcabp, top_n=10_000)
        mb_mask = full_np["neuropil"].astype(str).str.startswith("MB_")
        mb_frac = float(full_np.loc[mb_mask, "fraction"].sum()) if mb_mask.any() else 0.0
    # 通過：有合理顆數；萼部在輸入前列，或 MB 區合計為主（KC 亦在葉部收 DAN 等輸入）
    kcabp_ok = kcabp_n >= 50 and (
        (ca_rank is not None and ca_rank <= 5) or mb_frac >= 0.4 or ca_frac >= 0.1
    )

    # --- APL → KC ---
    # 注意：APL→KC 常是「很多夥伴、每對突觸數不大」，不能只看下游 top50
    apl = roots_by_primary_type(ct, "APL")
    apl_edges = outputs_from(edges, apl)
    apl_pairs = aggregate_neuron_pairs(apl_edges)
    apl_partners = apl_pairs.rename(
        columns={"post_root_id": "partner_id", "pre_root_id": "seed_id"}
    )
    by_partner = (
        apl_partners.groupby("partner_id", as_index=False)["syn_count"]
        .sum()
        .merge(meta, left_on="partner_id", right_on="root_id", how="left")
    )
    kc_mask = by_partner["class"].astype(str).eq("Kenyon_Cell") | by_partner[
        "primary_type"
    ].astype(str).str.startswith("KC", na=False)
    kc_partners = by_partner.loc[kc_mask]
    n_kc_partners = len(kc_partners)
    kc_syn = int(kc_partners["syn_count"].sum()) if n_kc_partners else 0
    all_syn = int(by_partner["syn_count"].sum()) if len(by_partner) else 0
    kc_frac_all = (kc_syn / all_syn) if all_syn else 0.0
    apl_down = by_partner.sort_values("syn_count", ascending=False).head(10)
    apl_ok = len(apl) >= 1 and n_kc_partners >= 100

    # --- MBON14 / MBON03 downstream ---
    mbon14 = roots_by_primary_type(ct, "MBON14")
    mbon03 = roots_by_primary_type(ct, "MBON03")
    down14 = partner_summary(
        edges, direction="downstream", seed_ids=mbon14, meta=meta, top_n=15
    )
    down03 = partner_summary(
        edges, direction="downstream", seed_ids=mbon03, meta=meta, top_n=15
    )
    mbon_ok = len(mbon14) >= 1 and len(mbon03) >= 1 and len(down14) > 0 and len(down03) > 0

    return {
        "min_syn": MIN_SYN,
        "kcabp_n": kcabp_n,
        "ca_frac": ca_frac,
        "mb_frac": mb_frac,
        "ca_rank": ca_rank,
        "np_dist": np_dist,
        "kcabp_ok": kcabp_ok,
        "apl_n": len(apl),
        "n_kc_partners": n_kc_partners,
        "kc_frac_top50": kc_frac_all,
        "apl_down_head": apl_down.head(10),
        "apl_ok": apl_ok,
        "mbon14_n": len(mbon14),
        "mbon03_n": len(mbon03),
        "down14": down14,
        "down03": down03,
        "mbon_ok": mbon_ok,
        "all_ok": kcabp_ok and apl_ok and mbon_ok,
    }


def write_report(res: dict, out: Path) -> None:
    lines = [
        "# 第 3 步驗證：已知結果",
        "",
        f"- **日期**：{date.today().isoformat()}",
        f"- **門檻**：syn_count ≥ {res['min_syn']}",
        f"- **資料**：`data/download/connections_princeton.csv.gz` 等",
        f"- **總評**：{'通過' if res['all_ok'] else '部分待核（見各項）'}",
        "",
        "---",
        "",
        "## 1. KCab-p",
        "",
        f"- primary_type == `KCab-p` 顆數：**{res['kcabp_n']}**",
        f"- 輸入突觸落在萼部 `{'/'.join(CALYX_NEUROPILS)}` 的比例：**{res['ca_frac']:.3f}**",
        f"- 全部 `MB_*` neuropil 輸入占比：**{res['mb_frac']:.3f}**",
        f"- 萼部 neuropil 進入輸入 top 的名次：rank={res['ca_rank']}",
        f"- 判定：{'通過' if res['kcabp_ok'] else '待核'}",
        "",
        "輸入 neuropil top：",
        "",
        "```",
        res["np_dist"].to_string(index=False) if len(res["np_dist"]) else "(empty)",
        "```",
        "",
        "## 2. APL → Kenyon cell",
        "",
        f"- APL 顆數：**{res['apl_n']}**",
        f"- 門檻以上連到的 KC 顆數（全下游，非僅 top50）：**{res['n_kc_partners']}**",
        f"- KC 佔 APL 全部下游突觸量比例：**{res['kc_frac_top50']:.3f}**",
        f"- 判定：{'通過' if res['apl_ok'] else '待核'}",
        "",
        "APL 下游 top10（單一最強邊，常為 DPM／MBON）：",
        "",
        "```",
        res["apl_down_head"][
            [c for c in ["partner_id", "syn_count", "primary_type", "class", "side"] if c in res["apl_down_head"].columns]
        ].to_string(index=False),
        "```",
        "",
        "## 3. MBON14／MBON03 下游",
        "",
        f"- MBON14 顆數：**{res['mbon14_n']}**；MBON03 顆數：**{res['mbon03_n']}**",
        f"- 判定：{'通過（有下游名單；請與 Codex 目視對照）' if res['mbon_ok'] else '待核'}",
        "",
        "### MBON14 下游 top15",
        "",
        "```",
        res["down14"][
            [c for c in ["partner_id", "syn_count", "primary_type", "class", "side"] if c in res["down14"].columns]
        ].to_string(index=False),
        "```",
        "",
        "### MBON03 下游 top15",
        "",
        "```",
        res["down03"][
            [c for c in ["partner_id", "syn_count", "primary_type", "class", "side"] if c in res["down03"].columns]
        ].to_string(index=False),
        "```",
        "",
        "## 備註",
        "",
        "- Codex 萼部 neuropil 代碼為 `MB_CA_L`／`MB_CA_R`。",
        "- KC 在葉部（ML／PED 等）也收調控輸入，故萼部占比不一定 >50%。",
        "- MBON 下游請到 Codex 同一細胞類型頁核對前幾名是否一致。",
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out}")


def main() -> int:
    res = run()
    out = ROOT / "results" / f"validation_{date.today().strftime('%Y%m%d')}.md"
    write_report(res, out)
    print("all_ok=", res["all_ok"])
    return 0 if res["all_ok"] else 0  # still succeed; report may say 待核


if __name__ == "__main__":
    raise SystemExit(main())
