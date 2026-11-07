# FlyWire 找 DAL

在 FlyWire v783 公開接線圖上，依 Chen et al. 2012 的結構條件尋找 DAL 候選神經元。只讀公開連線與註解表，不做濕實驗。

## 環境（conda）

```bash
cd FlyWire_找DAL
conda env create -f environment.yml
conda activate flywire-dal
```

若環境已存在、只要更新套件：

```bash
conda env update -f environment.yml --prune
```

主要套件：`pandas`、`networkx`、`navis`、`fafbseg`、`matplotlib`、`jupyter`。

## 目錄

| 路徑 | 用途 |
|---|---|
| `data/raw/` | 原始下載，永不手改 |
| `data/processed/` | 清洗／合併後的中間表 |
| `notebooks/` | 探索與驗證用 notebook |
| `src/` | 可重跑的篩選與工具函式 |
| `results/` | 候選清單、圖表、日誌摘要 |
| `scripts/` | 下載等一次性腳本 |

資料出處、版本與下載日期見 [`data/README.md`](data/README.md)。整體步驟見 [`專案規劃.md`](專案規劃.md)。

## 下載資料

主分析資料放在 `data/download/`（Codex 手動下載）。備用腳本：

```bash
conda activate flywire-dal
python scripts/download_data.py
```

## 第 3–4 步

```bash
conda activate flywire-dal
python scripts/run_validation.py          # 已知結果驗證 → results/validation_*.md
python scripts/run_dal_screen.py          # DAL 篩選 → results/dal_candidates_*.csv
# 放寬門檻例：
python scripts/run_dal_screen.py --min-syn-to-kcabp 3 --min-syn-from-mbon 3
```

## 授權提醒

使用 FlyWire 資料須遵守 [FlyWire Principles](https://flywire.ai/)。引用時請同時註明 materialization **v783** 與本專案 `data/README.md` 中的檔案出處。
