# FlyWire 找 DAL

在 FlyWire / Codex **FAFB materialization v783** 公開接線圖上，依 Chen et al. (2012) 的**結構條件**尋找 DAL **結構候選**。  
只讀公開連線與註解表，**不做**濕實驗、行為或電生理。

> 結論語境：本專案輸出的是「在這隻腦、這套門檻下最可能的結構候選」，**不是**「已鑑定為 Chen 2012 的 DAL」。

---

## 誰適合看這個 repo？

| 你想做的事 | 建議從哪裡開始 |
|---|---|
| 沿用同一套 DAL 結構篩選 | 下方「重現分析」→ `src/dal_filters.py`、`scripts/run_dal_screen.py` |
| 只想看結果／候選是誰 | [`results/`](results/)：先看 `core9_four_criteria.md`、`log_extension_latest.md` |
| 想改門檻／KC 定義做敏感度 | `run_dal_screen.py` 的 CLI（見下） |
| 想做類似「已知細胞驗證再篩選」 | `scripts/run_validation.py` → 再寫自己的 filter |
| 想對照 DAL2（負向對照） | `scripts/check_dal2_ppl1.py`、`scripts/run_dal2_screen.py` |
| 寫論文／報告引用方法 | [`專案規劃.md`](專案規劃.md)、[`data/README.md`](data/README.md)、本頁「文獻與命名」 |

---

## 研究問題（一句話）

哪些神經元同時（或近似）滿足：接收 MBON-α3／MBON-β′2mp 輸入、與萼部（`MB_CA_*`）Kenyon cells 有接觸、位置落在背前外側粗框——且可對突觸門檻做敏感度分析？

**目前最強關注候選（結構意義）：** `CB3198`、`CB2122`（其次 `CB0339`）。  
細節與左右腦、時鐘／MB 迴路延伸見 `results/log_extension_latest.md`。

---

## 倉庫裡有什麼／沒有什麼

| 有（可直接 clone） | 沒有（須自行下載，勿進 git） |
|---|---|
| 篩選與驗證程式 `src/`、`scripts/` | Codex 連線大表（數 GB～十數 GB） |
| 候選 CSV、研究日誌、延伸圖 `results/` | 神經骨架（skeletons，約 13 GB） |
| 資料出處與欄位說明 `data/README.md` | 可配準的 DAL 免疫／共聚焦影像（故未做 NBLAST） |
| conda 環境 `environment.yml` | |

大檔請放到本機 `data/download/`（目錄已預留；內容被 `.gitignore` 排除）。

---

## 環境

```bash
git clone https://github.com/GINOBYGINO/Search_DAL-_FLYWIRE.git
cd Search_DAL-_FLYWIRE

conda env create -f environment.yml
conda activate flywire-dal
```

更新既有環境：

```bash
conda env update -f environment.yml --prune
```

主要套件：`pandas`、`pyarrow`、`networkx`、`matplotlib`、`jupyter`；pip：`navis`、`fafbseg`（形態／NBLAST 預留，目前主流程不依賴骨架）。

---

## 資料（重現前必做）

1. 從 Codex 下載 FAFB v783 資料包：  
   https://codex.flywire.ai/api/download?dataset=fafb  
2. 解壓／放到 `data/download/`。主流程至少需要：

| 檔案 | 用途 |
|---|---|
| `connections_princeton_no_threshold.csv.gz` | 預設連線表（含弱邊；篩選預設） |
| `connections_princeton.csv.gz` | 有門檻連線（對照；`--use-thresholded`） |
| `neurons.csv.gz` | `nt_type` 等 |
| `classification.csv.gz` | 階層分類 |
| `consolidated_cell_types.csv.gz` | `primary_type` |
| `coordinates.csv.gz` | 細胞本體粗位置（soma 框） |

完整清單、欄位與「勿用」檔案見 [`data/README.md`](data/README.md)。  
備用註解母表腳本：`python scripts/download_data.py`。

使用資料須遵守 [FlyWire Principles](https://flywire.ai/)。引用時請註明 **v783** 與下載日期／檔名。

---

## 目錄

| 路徑 | 用途 |
|---|---|
| `data/raw/` | 原始下載，永不手改 |
| `data/processed/` | 清洗／合併後的中間表 |
| `data/download/` | Codex 主分析表（本機；未進版控） |
| `src/` | 可重用函式：`io`、`connectivity`、`dal_filters` |
| `scripts/` | 可重跑入口（驗證／篩選／DAL2／延伸） |
| `notebooks/` | 探索用 notebook |
| `results/` | 候選清單、日誌、圖表 |

整體設計步驟見 [`專案規劃.md`](專案規劃.md)。

---

## 重現分析（建議順序）

在專案根目錄、已 `conda activate flywire-dal`、且 `data/download/` 就緒後：

```bash
# 1) 已知結果驗證（確認讀檔與細胞類型沒搞錯）
python scripts/run_validation.py
# → results/validation_*.md

# 2) DAL 結構篩選（門檻可調）
python scripts/run_dal_screen.py
# → results/dal_candidates_*.csv

# 放寬／收緊例：
python scripts/run_dal_screen.py --min-syn-to-kc 1 --min-syn-from-mbon 1 --kc-mode kcab --tag r3
python scripts/run_dal_screen.py --use-thresholded --thresholds 5 --tag thresh5

# 3) DAL2 對照（勿與 DAL 身分混稱）
python scripts/check_dal2_ppl1.py
python scripts/run_dal2_screen.py

# 4) 延伸分析（左右腦、時鐘、巧合偵測、MBON 路徑等）
# 需先有相關連線子集；見 scripts/run_extension_*.py 與 results/log_extension_latest.md
```

### 篩選預設假設（沿用前請先讀）

定義在 `src/dal_filters.py` 的 `ScreenParams`，重點：

- 突觸邊預設讀 **無門檻** 連線表；各條件另有 `min_syn_*`
- KC：預設 `KCab-p`；可用 `--kc-mode kcab` 扩成所有 `KCab*`
- MBON：Codex 名稱 **`MBON14`**＝文獻 α3、**`MBON03`**＝β′2mp；預設限 CRE／SMP／SIP neuropil
- 排除：Kenyon_Cell、MBON、DAN、ALPN、LHCENT、visual_projection、APL／DPM／LHMB1 等
- soma 背前外側粗框：**預設關閉**（避免誤殺）；需要時加 `--soma`
- `nt_type` 為 Codex **預測**，預設不強制白名單

改門檻後請另存 `results/` 並在日誌寫參數（本專案以檔名 tag／日期區分各輪）。

---

## 結果怎麼讀（給只想參考的人）

建議閱讀順序：

1. **`results/validation_20261019.md`** — 資料讀取是否可信  
2. **`results/log_20261019_r3.md`** — 主篩選輪次與參數  
3. **`results/core9_four_criteria.md`** — 九顆候選：位置／形態／MB_CA／NT 四點對照  
4. **`results/log_extension_latest.md`** — 左腦、時鐘、巧合、MB 迴路、睡眠路徑摘要  
5. 對應 CSV：`dal_candidates_*.csv`、`extension_*.csv`、圖在 `results/extension_figures/`

**CSV 欄位提示：** `pass_*`＝是否通過該條件；`score_*`＝突觸量或加總分數；`primary_type`／`root_id` 可回 Codex 查形態。

---

## 文獻與命名（引用／改方法時）

| 主題 | 參考 |
|---|---|
| DAL 結構條件來源 | Chen et al., 2012（長期記憶形成可視化；DAL） |
| 全腦接線圖 v783 | Dorkenwald et al., 2024, *Nature* |
| 細胞類型註解 | Schlegel et al., 2024, *Nature* |
| MBON–DAL／GRASP 動機 | Wu et al., 2017；路徑寫法可對照 Giakoumas et al.（咽／腸感覺圖譜） |
| 延伸分析動機 | Yin 2025（時鐘）、Xia 2005（巧合偵測動機）、Wu 2017（第三種 MBON → Codex **MBON09**） |
| DAL2 對照 | Feng et al., 2021；本專案僅作方法／負向對照 |

Codex 名稱與文獻別稱常不一致（例如 MBON14／03／09）。改腳本時請以 `consolidated_cell_types.csv.gz` 的 `primary_type` 為準，並在方法寫明對照表。

---

## 已知限制（請一併寫進你的方法／討論）

- 單一隻腦（v783）；結構可達 ≠ 功能必要  
- `nt_type` 預測與文獻 DAL 的 DDC 陽性**可能衝突**（關注候選多為 ACH）  
- 無合法可配準 DAL 影像前，**不宣稱 NBLAST 鑑定**  
- 突觸門檻、soma 框、排除類型皆為作者定義假設 → 需敏感度  
- ≥5 門檻常切掉弱邊；左右腦對稱性並非每型都成立（見延伸日誌）

---

## 授權與引用

- FlyWire 資料：[FlyWire Principles](https://flywire.ai/)  
- 請同時註明 materialization **v783**、Codex 檔名與本 repo 的 `data/README.md`／對應 `results/` 檔  
- 若沿用本專案的篩選定義或候選清單，請連結本 repository，並區分「結構候選」與「功能鑑定」

---

## 快速對照：腳本入口

| 腳本 | 作用 |
|---|---|
| `scripts/run_validation.py` | KCab-p／APL／MBON 下游等已知檢查 |
| `scripts/run_dal_screen.py` | DAL 結構篩選主程式 |
| `scripts/check_dal2_ppl1.py` | DAL2 ↔ PPL1 連結檢查 |
| `scripts/run_dal2_screen.py` | DAL2 式雙輸出篩選對照 |
| `scripts/run_extension_*.py` | 左右腦、時鐘、巧合、MBON／睡眠路徑延伸 |
| `src/dal_filters.py` | 條件函式與 `ScreenParams`（改邏輯從這裡） |
| `src/connectivity.py` / `src/io.py` | 連線過濾與讀檔 |
