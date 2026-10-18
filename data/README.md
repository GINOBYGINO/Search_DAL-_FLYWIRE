# 資料出處紀錄

專案 materialization：**v783**
註解鎖定：**flywire_annotations v2.1.0**（對齊 Nature 2024 公開版）
本 README 最近更新下載日：**2026-10-19**

原始檔放在 `raw/`，**不要手改**。由註解母表拆出的表在 `processed/`。報告引用請依下表。

## 為什麼不用 Codex API 直接抓 neurons？

Codex `download_resource` 現在需要帳號 `api_token`。公開的 `flywire_annotations v2.1.0` 註解表已含 `top_nt`／`top_nt_conf`以及 classification、cell_type 欄位，與論文公開版一致，故以此為準。

| 本機路徑 | 版本 | 下載日期 | 來源／網址 | 用途 |
|---|---|---|---|---|
| `raw/proofread_connections_783.feather` (261,586,944 bytes; md5 `a4f36b87176787f3a909ba9ac52e44f0`) | 783.0 | 2026-10-19 | Zenodo record 10676866 (FlyWire Whole-brain Connectome Connectivity Data)<br>`https://zenodo.org/records/10676866/files/proofread_connections_783.feather?download=1` | 校對後神經元對神經元連線（neuropil、突觸數、傳導物質平均機率） |
| `raw/Supplemental_file1_neuron_annotations.tsv` (27,015,208 bytes; md5 `8527e0a95ed5f112766b13260a91e8e2`) | v2.1.0 | 2026-10-19 | flyconnectome/flywire_annotations v2.1.0<br>`https://raw.githubusercontent.com/flyconnectome/flywire_annotations/v2.1.0/supplemental_files/Supplemental_file1_neuron_annotations.tsv` | 神經元註解母表（含 top_nt、flow／super_class／cell_class、cell_type 等；等同 Codex neurons + classification + cell_types 的公開對應） |
| `processed/neurons_783.csv` (19,566,067 bytes; md5 `46cfab57d040c88eedaca72f03363a56`) | v2.1.0 | 2026-10-19 | flyconnectome/flywire_annotations v2.1.0（本機拆表）<br>`https://raw.githubusercontent.com/flyconnectome/flywire_annotations/v2.1.0/supplemental_files/Supplemental_file1_neuron_annotations.tsv` | 由註解母表匯出（含 top_nt） |
| `processed/classification_783.csv` (8,650,414 bytes; md5 `4aa7e93035147e334198eaecb30fb348`) | v2.1.0 | 2026-10-19 | flyconnectome/flywire_annotations v2.1.0（本機拆表）<br>`https://raw.githubusercontent.com/flyconnectome/flywire_annotations/v2.1.0/supplemental_files/Supplemental_file1_neuron_annotations.tsv` | 由註解母表匯出（flow／super_class／cell_class 等） |
| `processed/cell_types_783.csv` (6,293,301 bytes; md5 `065dd60eaef121a69003a0756983ebf5`) | v2.1.0 | 2026-10-19 | flyconnectome/flywire_annotations v2.1.0（本機拆表）<br>`https://raw.githubusercontent.com/flyconnectome/flywire_annotations/v2.1.0/supplemental_files/Supplemental_file1_neuron_annotations.tsv` | 由註解母表匯出（cell_type／hemibrain_type 等） |

## 重跑下載

```bash
conda activate flywire-dal
python scripts/download_data.py
# 強制重抓：
python scripts/download_data.py --force
```

## 未下載

- 完整突觸點表 `flywire_synapses_783.feather`（Zenodo，約 9.5 GB）：第 3 步若需要再補。
