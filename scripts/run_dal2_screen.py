"""DAL2 全庫連結性篩選（僅 PPL105＋PPL106 雙輸出硬條件）。

依據 Feng et al. 2021：DAL2 輸出到 PPL1-α2α′2（FlyWire PPL105）與 PPL1-α3（PPL106）。
這是連結性候選，不能宣稱找到 DAL2；NPF 無法用 Codex NT 驗證。
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.connectivity import filter_edges  # noqa: E402
from src.dal_filters import parse_position  # noqa: E402
from src.io import (  # noqa: E402
    load_cell_types,
    load_classification,
    load_connections,
    load_coordinates,
    load_neurons,
    roots_by_primary_type,
)

# Codex：PPL105 = a`2a2；PPL106 = a3
PPL105 = "PPL105"
PPL106 = "PPL106"

# 先前 DAL 篩選排除＋所有 PPL*
EXCLUDE_CLASSES = ("Kenyon_Cell", "MBON", "DAN", "ALPN", "LHCENT")
EXCLUDE_SUPER = ("visual_projection",)
EXCLUDE_PRIMARY = ("APL", "DPM", "LHMB1")

# 背前外側粗框（與先前 DAL 相同）
MID, LAT = 512_000.0, 60_000.0
Y_LO, Y_HI, Z_LO = 110_000.0, 240_000.0, 100_000.0

# 先前關注的細胞
WATCH = {
    720575940642272456: "SMP194",
    720575940628347654: "CB2122",
    720575940633320785: "CB3198",
}


def _weights_to_set(edges: pd.DataFrame, post_ids: list[int]) -> pd.Series:
    """pre → 對 post 集合的 syn 總和。"""
    sub = edges.loc[edges["post_root_id"].isin(set(int(x) for x in post_ids))]
    if sub.empty:
        return pd.Series(dtype=float)
    return sub.groupby("pre_root_id")["syn_count"].sum()


def screen(edges: pd.DataFrame, min_syn_to_each: int) -> pd.DataFrame:
    ct = load_cell_types()
    cl = load_classification()
    neu = load_neurons()

    ids_105 = roots_by_primary_type(ct, PPL105)
    ids_106 = roots_by_primary_type(ct, PPL106)
    if not ids_105 or not ids_106:
        raise RuntimeError(f"缺少目標：PPL105={len(ids_105)} PPL106={len(ids_106)}")

    w105 = _weights_to_set(edges, ids_105)
    w106 = _weights_to_set(edges, ids_106)

    # 硬條件：兩邊都 ≥ 門檻
    common = sorted(set(w105[w105 >= min_syn_to_each].index.astype(int)) & set(w106[w106 >= min_syn_to_each].index.astype(int)))
    out = pd.DataFrame({"root_id": common})
    out["syn_to_PPL105"] = out["root_id"].map(w105).fillna(0).astype(int)
    out["syn_to_PPL106"] = out["root_id"].map(w106).fillna(0).astype(int)
    out["syn_sum_both"] = out["syn_to_PPL105"] + out["syn_to_PPL106"]

    out = out.merge(ct[["root_id", "primary_type"]], on="root_id", how="left")
    out = out.merge(
        cl[["root_id", "super_class", "class", "side"]],
        on="root_id",
        how="left",
    )
    out = out.merge(neu[["root_id", "nt_type"]], on="root_id", how="left")
    out = out.rename(columns={"primary_type": "type", "nt_type": "predicted_NT"})

    # 排除
    out = out.loc[~out["class"].astype(str).isin(EXCLUDE_CLASSES)]
    out = out.loc[~out["super_class"].astype(str).isin(EXCLUDE_SUPER)]
    out = out.loc[~out["type"].astype(str).isin(EXCLUDE_PRIMARY)]
    # 所有 PPL*
    out = out.loc[~out["type"].astype(str).str.startswith("PPL", na=False)]

    # soma
    coords = load_coordinates()
    xyz = parse_position(coords["position"])
    coords = (
        pd.concat([coords[["root_id"]], xyz], axis=1)
        .dropna(subset=["x", "y", "z"])
        .drop_duplicates("root_id", keep="first")
    )
    out = out.merge(coords, on="root_id", how="left")
    has = out[["x", "y", "z"]].notna().all(axis=1)
    out["in_DAL_box"] = False
    out.loc[has, "in_DAL_box"] = (
        ((out.loc[has, "x"] - MID).abs() >= LAT)
        & (out.loc[has, "y"] >= Y_LO)
        & (out.loc[has, "y"] <= Y_HI)
        & (out.loc[has, "z"] >= Z_LO)
    )

    out["watch"] = out["root_id"].map(WATCH).fillna("")
    out["rank"] = out["syn_sum_both"].rank(method="min", ascending=False).astype(int)
    out = out.sort_values(
        ["syn_sum_both", "syn_to_PPL105", "syn_to_PPL106"],
        ascending=False,
    ).reset_index(drop=True)
    # 重新依排序給 1-based rank
    out["rank"] = range(1, len(out) + 1)

    cols = [
        "rank",
        "type",
        "root_id",
        "side",
        "syn_to_PPL105",
        "syn_to_PPL106",
        "syn_sum_both",
        "predicted_NT",
        "x",
        "y",
        "z",
        "in_DAL_box",
        "watch",
        "class",
        "super_class",
    ]
    return out[[c for c in cols if c in out.columns]]


def _side_counts(df: pd.DataFrame) -> str:
    s = df["side"].astype(str).str.lower().value_counts()
    left = int(s.get("left", 0))
    right = int(s.get("right", 0))
    other = int(len(df) - left - right)
    return f"left={left}, right={right}, other/NA={other}"


def main() -> int:
    print("loading connections_princeton_no_threshold…")
    conn = load_connections(no_threshold=True)
    today = date.today().strftime("%Y%m%d")
    results: dict[int, pd.DataFrame] = {}

    for thr in (5, 1):
        print(f"screening min_syn_to_each≥{thr}…")
        edges = filter_edges(conn, min_syn=1)  # 用原始邊；門檻加在「對目標集合加總」
        # 實際上對集合加總後再比 thr；edge 先保留 ≥1 即可
        df = screen(edges, min_syn_to_each=thr)
        results[thr] = df
        out_csv = ROOT / "results" / f"dal2_screen_min{thr}.csv"
        df.to_csv(out_csv, index=False, encoding="utf-8-sig")
        print(f"wrote {out_csv} ({len(df)} rows)")

    # log
    lines = [
        f"# DAL2 全庫篩選（僅 PPL105＋PPL106 雙輸出）{today}",
        "",
        "- **資料**：FlyWire v783 `connections_princeton_no_threshold.csv.gz`",
        "- **硬條件**：同時輸出到 PPL105（α2α′2）與 PPL106（α3）；各達門檻",
        "- **排除**：Kenyon_Cell／MBON／DAN／ALPN／LHCENT／visual_projection／APL／DPM／LHMB1／所有 PPL*",
        "- **soma 粗框**：|x−512k|≥60k、y∈[110k,240k]、z≥100k（僅標記，不刪除）",
        "- **聲明**：連結性候選 ≠ 找到 DAL2；NPF 無法用 Codex NT 驗證。",
        "",
    ]
    for thr, df in results.items():
        watch_rows = df.loc[df["watch"] != ""]
        smp = df.loc[df["watch"] == "SMP194"]
        smp_rank = int(smp["rank"].iloc[0]) if len(smp) else None
        smp_sum = int(smp["syn_sum_both"].iloc[0]) if len(smp) else None
        lines.extend(
            [
                f"## 門檻 ≥{thr}（對 PPL105、PPL106 各自加總）",
                "",
                f"- **通過數**：{len(df)}",
                f"- **側別**：{_side_counts(df)}",
                f"- **背前外側粗框內**：{int(df['in_DAL_box'].sum())}",
                f"- **SMP194 排名**：{smp_rank}（syn_sum={smp_sum}）" if smp_rank else "- **SMP194**：未通過本門檻",
                "",
                "### 關注細胞",
                "",
                "```",
                watch_rows.to_string(index=False) if len(watch_rows) else "(無)",
                "```",
                "",
                "### 前 20 名",
                "",
                "```",
                df.head(20).to_string(index=False),
                "```",
                "",
            ]
        )

    out_log = ROOT / "results" / f"log_{today}_dal2_screen.md"
    # 固定檔名也寫一份方便找
    out_log2 = ROOT / "results" / "log_dal2_screen.md"
    text = "\n".join(lines)
    out_log.write_text(text, encoding="utf-8")
    out_log2.write_text(text, encoding="utf-8")
    print(f"wrote {out_log}")

    # 終端摘要給使用者
    for thr, df in results.items():
        print(f"\n===== min{thr} n={len(df)} sides={_side_counts(df)} =====")
        smp = df.loc[df["watch"] == "SMP194"]
        if len(smp):
            print(f"SMP194 rank={int(smp['rank'].iloc[0])} sum={int(smp['syn_sum_both'].iloc[0])}")
        else:
            print("SMP194: not in list")
        print(df.head(20)[["rank", "type", "root_id", "side", "syn_to_PPL105", "syn_to_PPL106", "syn_sum_both", "predicted_NT", "in_DAL_box", "watch"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
