"""查 CB 候選 → PPL1-α2α′2（PPL105）／PPL1-α3（PPL106）。

DAL2（Feng et al. 2021）輸出到這兩類；有邊較像 DAL2，無邊較留在 DAL 候選。
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.connectivity import filter_edges  # noqa: E402
from src.io import load_cell_types, load_connections, roots_by_primary_type  # noqa: E402

# Codex labels：PPL105 = a`2a2（α′2α2／α2α′2）；PPL106 = a3（α3）
TARGET_MAP = {
    "PPL1-α2α′2": "PPL105",
    "PPL1-α3": "PPL106",
}

# r3 KCab* 結構三分排名前 9
CANDIDATES = {
    "1_LHPV5e1_R": 720575940622999688,
    "2_LHPV5e1_L": 720575940624155672,
    "3_AL-MBDL1_R": 720575940624206472,
    "4_AL-MBDL1_L": 720575940648669060,
    "5_OA-VPM3": 720575940640356149,
    "6_SMP194": 720575940642272456,
    "7_CB2122": 720575940628347654,
    "8_CB0339": 720575940632389028,
    "9_CB3198": 720575940633320785,
}


def main() -> int:
    ct = load_cell_types()
    targets_rows = []
    for group, ptype in TARGET_MAP.items():
        ids = roots_by_primary_type(ct, ptype)
        for rid in ids:
            targets_rows.append(
                {
                    "dal2_target_group": group,
                    "primary_type": ptype,
                    "root_id": rid,
                }
            )
    targets = pd.DataFrame(targets_rows)
    print("targets:")
    print(targets.to_string(index=False))

    target_ids = set(targets["root_id"].astype(int))
    group_of = {
        int(r.root_id): r.dal2_target_group for _, r in targets.iterrows()
    }

    print("loading no_threshold connections…")
    conn = load_connections(no_threshold=True)
    edges = filter_edges(conn, min_syn=1)
    cand_ids = list(CANDIDATES.values())

    # 主檢驗：候選 → 目標
    sub = edges.loc[
        edges["pre_root_id"].isin(cand_ids) & edges["post_root_id"].isin(target_ids)
    ].copy()
    print(f"edges cand→PPL105/106: {len(sub)}")

    summary_rows = []
    detail_rows = []
    for name, rid in CANDIDATES.items():
        for group, ptype in TARGET_MAP.items():
            tids = set(roots_by_primary_type(ct, ptype))
            e = edges.loc[
                (edges["pre_root_id"] == rid) & (edges["post_root_id"].isin(tids))
            ]
            total = int(e["syn_count"].sum()) if len(e) else 0
            n_partners = int(e["post_root_id"].nunique()) if len(e) else 0
            summary_rows.append(
                {
                    "candidate": name,
                    "root_id": rid,
                    "target_group": group,
                    "target_primary_type": ptype,
                    "n_partners": n_partners,
                    "syn_count_sum": total,
                    "has_output": total > 0,
                }
            )
            if len(e):
                by = (
                    e.groupby(["post_root_id", "neuropil"], as_index=False)["syn_count"]
                    .sum()
                    .sort_values("syn_count", ascending=False)
                )
                by["candidate"] = name
                by["pre_root_id"] = rid
                by["target_group"] = group
                by["post_primary_type"] = ptype
                detail_rows.append(by)

    summary = pd.DataFrame(summary_rows)
    detail = pd.concat(detail_rows, ignore_index=True) if detail_rows else pd.DataFrame()

    # 反向也記一筆（GRASP 不指定方向）
    rev_rows = []
    for name, rid in CANDIDATES.items():
        for group, ptype in TARGET_MAP.items():
            tids = set(roots_by_primary_type(ct, ptype))
            e = edges.loc[
                (edges["post_root_id"] == rid) & (edges["pre_root_id"].isin(tids))
            ]
            rev_rows.append(
                {
                    "candidate": name,
                    "root_id": rid,
                    "target_group": group,
                    "direction": "PPL1→cand",
                    "syn_count_sum": int(e["syn_count"].sum()) if len(e) else 0,
                    "n_partners": int(e["pre_root_id"].nunique()) if len(e) else 0,
                }
            )
    reverse = pd.DataFrame(rev_rows)

    # 診斷：候選 → 任何 PPL*
    ppl_ids = set(
        ct.loc[ct["primary_type"].astype(str).str.startswith("PPL", na=False), "root_id"]
        .astype(int)
    )
    diag = edges.loc[
        edges["pre_root_id"].isin(cand_ids) & edges["post_root_id"].isin(ppl_ids)
    ].copy()
    if len(diag):
        diag = diag.merge(
            ct[["root_id", "primary_type"]],
            left_on="post_root_id",
            right_on="root_id",
            how="left",
            suffixes=("", "_y"),
        )
        diag_sum = (
            diag.groupby(["pre_root_id", "primary_type"], as_index=False)["syn_count"]
            .sum()
            .sort_values("syn_count", ascending=False)
        )
        diag_sum["candidate"] = diag_sum["pre_root_id"].map(
            {v: k for k, v in CANDIDATES.items()}
        )
    else:
        diag_sum = pd.DataFrame()

    by_cand = (
        summary.groupby(["candidate", "root_id"], as_index=False)
        .agg(syn_to_either=("syn_count_sum", "sum"), n_groups_hit=("has_output", "sum"))
    )
    def _interp(x: int) -> str:
        if x >= 10:
            return "較像 DAL2（雙標的輸出較重）"
        if x > 0:
            return "弱接觸→PPL1（1–數突觸；DAL2 證據弱）"
        return "無→PPL105/106（較留在 DAL 候選）"

    by_cand["interpretation"] = by_cand["syn_to_either"].map(_interp)

    today = date.today().strftime("%Y%m%d")
    summary.to_csv(
        ROOT / "results" / f"dal2_ppl1_check_{today}.csv",
        index=False,
        encoding="utf-8-sig",
    )
    targets.to_csv(
        ROOT / "results" / f"dal2_ppl1_targets_{today}.csv",
        index=False,
        encoding="utf-8-sig",
    )
    if len(detail):
        detail.to_csv(
            ROOT / "results" / f"dal2_ppl1_check_{today}_detail.csv",
            index=False,
            encoding="utf-8-sig",
        )
    reverse.to_csv(
        ROOT / "results" / f"dal2_ppl1_check_{today}_reverse.csv",
        index=False,
        encoding="utf-8-sig",
    )

    log = "\n".join(
        [
            f"# DAL2 檢驗：CB → PPL1-α2α′2／PPL1-α3（{today}）",
            "",
            "- **依據**：Feng et al. 2021 — DAL2 軸突到 PPL1-α2α′2、PPL1-α3。",
            "- **FlyWire 對照**（Codex labels）：`PPL105` = a\`2a2（α2α′2）；`PPL106` = a3（α3）。",
            "- **資料**：`connections_princeton_no_threshold.csv.gz`，突觸門檻 ≥1。",
            "- **候選**：r3 前 9 顆（LHPV5e1×2、AL-MBDL1×2、OA-VPM3、SMP194、CB2122、CB0339、CB3198）。",
            "- **注意**：不宣稱 Feng 2021 已在 hemibrain 鑑定 DAL2；ACh 預測不排除 DAL。",
            "",
            "## 目標細胞",
            "",
            "```",
            targets.to_string(index=False),
            "```",
            "",
            "## 候選 → 目標（主檢驗）",
            "",
            "```",
            summary.to_string(index=False),
            "```",
            "",
            "## 判定",
            "",
            "```",
            by_cand.to_string(index=False),
            "```",
            "",
            "## 反向（PPL1 → 候選，僅供參考）",
            "",
            "```",
            reverse.to_string(index=False),
            "```",
            "",
            "## 診斷：候選 → 任何 PPL*",
            "",
            "```",
            diag_sum.to_string(index=False) if len(diag_sum) else "(無)",
            "```",
            "",
            "## 細節（有邊）",
            "",
            "```",
            detail.to_string(index=False) if len(detail) else "(無 cand→PPL105/106 突觸)",
            "```",
            "",
        ]
    )
    out_log = ROOT / "results" / f"log_{today}_dal2_ppl1.md"
    out_log.write_text(log, encoding="utf-8")
    print(summary.to_string(index=False))
    print(by_cand.to_string(index=False))
    print("wrote", out_log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
