"""第 4 步：跑 DAL 結構候選篩選，寫 CSV 與研究日誌。"""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.dal_filters import (  # noqa: E402
    MBON_DAL_NEUROPILS,
    ScreenParams,
    disable_soma,
    enable_soma_dal_box,
    params_to_markdown,
    run_screen,
)
from src.io import load_connections  # noqa: E402


def _result_block(df, cols: list[str], top: int) -> dict:
    core = df.loc[df["score_structural"] == 3]
    return {
        "n": len(df),
        "core": core,
        "pair_km": int(((df["pass_kcabp"] == 1) & (df["pass_mbon"] == 1)).sum()),
        "pair_kc": int(((df["pass_kcabp"] == 1) & (df["pass_calyx"] == 1)).sum()),
        "pair_mc": int(((df["pass_mbon"] == 1) & (df["pass_calyx"] == 1)).sum()),
        "n_kcabp": int(df["pass_kcabp"].sum()),
        "n_to": int(df["pass_to_kcabp"].sum()) if "pass_to_kcabp" in df.columns else None,
        "n_from": int(df["pass_from_kcabp"].sum()) if "pass_from_kcabp" in df.columns else None,
        "n_mbon": int(df["pass_mbon"].sum()),
        "n_calyx": int(df["pass_calyx"].sum()),
        "head": df.head(top),
        "cols": cols,
        "n_kc_targets": int(df["n_kc_targets"].iloc[0]) if len(df) and "n_kc_targets" in df.columns else None,
    }


def _cols(df) -> list[str]:
    return [
        c
        for c in [
            "root_id",
            "primary_type",
            "class",
            "super_class",
            "side",
            "nt_type",
            "x",
            "y",
            "z",
            "pass_kcabp",
            "pass_to_kcabp",
            "pass_from_kcabp",
            "pass_mbon",
            "pass_calyx",
            "pass_nt",
            "pass_soma",
            "score_to_kcabp",
            "score_from_kcabp",
            "score_kcabp",
            "score_from_mbon",
            "score_calyx",
            "score_structural",
            "score_pass_sum",
            "score_total",
            "reason",
        ]
        if c in df.columns
    ]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--min-syn", type=int, default=None)
    p.add_argument("--min-syn-edge", type=int, default=None)
    p.add_argument("--min-syn-to-kc", type=int, default=None)
    p.add_argument("--min-syn-from-mbon", type=int, default=None)
    p.add_argument("--min-ca-syn", type=int, default=None)
    p.add_argument("--thresholds", nargs="+", type=int, default=None)
    p.add_argument("--kc-mode", choices=["kcabp", "kcab"], default="kcabp")
    p.add_argument(
        "--kc-modes",
        nargs="+",
        choices=["kcabp", "kcab"],
        default=None,
        help="一次跑多個 KC 模式",
    )
    p.add_argument("--no-threshold", action="store_true", default=True)
    p.add_argument("--use-thresholded", action="store_true", help="改用 connections_princeton.csv.gz")
    p.add_argument("--mbon-neuropil", action="store_true", default=True, help="MBON 限 CRE/SMP/SIP")
    p.add_argument("--mbon-any-neuropil", action="store_true", help="MBON 不限 neuropil")
    p.add_argument("--soma", action="store_true", help="開啟背前外側 soma 粗框")
    p.add_argument("--no-soma", action="store_true", help="關閉 soma（預設）")
    p.add_argument("--one-way-kc", action="store_true")
    p.add_argument("--nt", nargs="*", default=[])
    p.add_argument("--top", type=int, default=40)
    p.add_argument("--tag", type=str, default="r3")
    args = p.parse_args()

    no_threshold = not args.use_thresholded
    mbon_neuropils: tuple[str, ...] = () if args.mbon_any_neuropil else MBON_DAL_NEUROPILS
    use_soma = args.soma and not args.no_soma

    thresholds = args.thresholds
    if thresholds is None:
        base = args.min_syn if args.min_syn is not None else 2
        thresholds = [base]
    kc_modes = args.kc_modes or [args.kc_mode]

    today = date.today().strftime("%Y%m%d")
    print(
        f"loading connections ({'no_threshold' if no_threshold else 'princeton≥filter'})…"
    )
    conn = load_connections(no_threshold=no_threshold)
    print(f"connections rows={len(conn)} syn_min={int(conn.syn_count.min())}")

    blocks: list[tuple[str, int, ScreenParams, object]] = []

    for kc_mode in kc_modes:
        for thr in thresholds:
            params = ScreenParams(
                min_syn_edge=args.min_syn_edge if args.min_syn_edge is not None else thr,
                min_syn_to_kc=args.min_syn_to_kc if args.min_syn_to_kc is not None else thr,
                min_syn_from_mbon=(
                    args.min_syn_from_mbon if args.min_syn_from_mbon is not None else thr
                ),
                min_ca_syn=args.min_ca_syn if args.min_ca_syn is not None else thr,
                no_threshold=no_threshold,
                kc_mode=kc_mode,
                mbon_neuropils=mbon_neuropils,
                nt_whitelist=tuple(args.nt),
                bidirectional_kc=not args.one_way_kc,
            )
            if use_soma:
                enable_soma_dal_box(params)
            else:
                disable_soma(params)

            label = f"kc={kc_mode} min_syn={thr}"
            print(f"running DAL screen ({label})…")
            df = run_screen(params, connections=conn)
            out_csv = (
                ROOT
                / "results"
                / f"dal_candidates_{today}_{args.tag}_{kc_mode}_min{thr}.csv"
            )
            cols = _cols(df)
            df[cols].to_csv(out_csv, index=False)
            print(f"wrote {out_csv} ({len(df)} rows)")
            blocks.append((kc_mode, thr, params, _result_block(df, cols, args.top)))

    out_log = ROOT / "results" / f"log_{today}_{args.tag}.md"
    lines = [
        f"# 研究日誌：DAL 篩選 {today}（{args.tag}）",
        "",
        f"- **日期**：{date.today().isoformat()}",
        f"- **資料**：`connections_princeton{'_no_threshold' if no_threshold else ''}.csv.gz`",
        "- **腳本**：`scripts/run_dal_screen.py`",
        "- **本輪變更**：改用 no_threshold；關閉 soma 粗框；MBON 限 CRE/SMP/SIP；"
        "分別測 KCab-p 與全部 KCab*；門檻 1–2",
        "",
    ]
    for kc_mode, thr, params, b in blocks:
        lines.extend(
            [
                f"## kc_mode={kc_mode}；min_syn={thr}",
                "",
                f"- **KC 目標數**：{b['n_kc_targets']}",
                f"- **候選總列數**：{b['n']}",
                f"- **結構三分皆過**：{len(b['core'])}",
                f"- **成對**：KC∩MBON={b['pair_km']}；KC∩calyx={b['pair_kc']}；MBON∩calyx={b['pair_mc']}",
                f"- **單條件**：KC={b['n_kcabp']}（→KC={b['n_to']}；KC→={b['n_from']}）；"
                f"MBON={b['n_mbon']}；calyx={b['n_calyx']}",
                "",
                "### 參數",
                "",
                params_to_markdown(params),
                "",
                "### 核心候選（結構三分）",
                "",
                "```",
                (
                    b["core"].head(args.top)[b["cols"]].to_string(index=False)
                    if len(b["core"])
                    else "(無)"
                ),
                "```",
                "",
                "### 總分前幾名",
                "",
                "```",
                b["head"][b["cols"]].to_string(index=False),
                "```",
                "",
            ]
        )
    lines.extend(
        [
            "## 備註",
            "",
            "- 結論僅為**結構候選**，非 Chen 2012 功能定義的 DAL。",
            "- MBON14=α3、MBON03=β′2mp（Codex labels 已核對）；MBON 突觸限 CRE/SMP/SIP。",
            "- 未開 `nt_whitelist`；未做骨架／NBLAST。",
            "",
        ]
    )
    out_log.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {out_log}")
    for kc_mode, thr, _, b in blocks:
        print(
            f"{kc_mode} min{thr}: n={b['n']} core={len(b['core'])} "
            f"KM={b['pair_km']} KC={b['n_kcabp']} M={b['n_mbon']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
