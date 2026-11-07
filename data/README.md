# 資料出處紀錄

專案 materialization：**v783**  
註解對照：**flywire_annotations v2.1.0**（Nature 2024）＋ **Codex FAFB v783 下載包**  
本 README 最近更新：**2026-10-19**

原始檔**不要手改**。報告引用請依下表。

---

## 主資料：`download/`（Codex 手動下載，第 2 步以此為準）

來源入口：https://codex.flywire.ai/api/download?dataset=fafb  
（使用者於 2026-10-19 下載；未含神經骨架／skeletons，約 13 GB，暫不需要。）

| 本機檔名 | 用途（對應第 2 步） |
|---|---|
| `connections_princeton.csv.gz` | 主連線表（pre／post／neuropil／syn_count／nt_type） |
| `connections_princeton_no_threshold.csv.gz` | 無門檻連線（對照用） |
| `neurons.csv.gz` | 神經元表；傳導物質欄為 **`nt_type`**（等同計畫中的 top_nt） |
| `classification.csv.gz` | 階層分類（flow／super_class／class／…） |
| `consolidated_cell_types.csv.gz` | 細胞類型（primary_type 等） |
| `fafb_v783_princeton_synapse_table.csv.gz` | 突觸級大表（約 2.5 GB；第 3 步若需要再詳用） |
| `synapse_coordinates.csv.gz` 等其餘檔 | 座標、標籤、neuropil 統計等輔助 |

欄位抽樣確認（2026-10-19）：

- `neurons.csv.gz`：`root_id`, `nt_type`, `nt_type_score`, …
- `classification.csv.gz`：`root_id`, `flow`, `super_class`, `class`, …
- `consolidated_cell_types.csv.gz`：`root_id`, `primary_type`, …
- `connections_princeton.csv.gz`：`pre_root_id`, `post_root_id`, `neuropil`, `syn_count`, `nt_type`

---

## 備用：`raw/` 與 `processed/`

| 本機路徑 | 狀態 | 來源 | 用途 |
|---|---|---|---|
| `raw/Supplemental_file1_neuron_annotations.tsv` | 完整 | [flywire_annotations v2.1.0](https://raw.githubusercontent.com/flyconnectome/flywire_annotations/v2.1.0/supplemental_files/Supplemental_file1_neuron_annotations.tsv) | 含 `top_nt` 的註解母表；可與 Codex 交叉核對 |
| `processed/neurons_783.csv` 等 | 由母表拆出 | 同上 | 備用；正式分析優先用 `download/` |
| `raw/proofread_connections_783.feather` | **不完整／勿用**（約 370 MB／852 MB） | Zenodo 10676866 | 已改以 Codex `connections_princeton` 取代 |

---

## 刻意未下載

- **神經骨架（skeletons）**：約 13 GB。第 3–4 步若要做形態／NBLAST 再補。
- Zenodo 完整 `proofread_connections_783.feather`／`flywire_synapses_783.feather`：目前不需要。

---

## 分析時建議讀哪些檔

1. 連線：`download/connections_princeton.csv.gz`
2. 傳導物質：`download/neurons.csv.gz` 的 `nt_type`（或註解母表的 `top_nt`）
3. 分類／類型：`download/classification.csv.gz`、`download/consolidated_cell_types.csv.gz`
