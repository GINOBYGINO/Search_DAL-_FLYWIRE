# 延伸分析日誌 20261101

- 資料：v783；分塊過濾後 `_relevant_conn.parquet`
- 只寫結構候選；`nt_type` 為預測值
- 圖：本機 matplotlib savefig 崩潰，改以 PIL 輸出 PNG

## 1. 左腦補齊
- Codex：CB3198/CB2122 各 4 顆（左右各 2）；CB0339 左右各 1。
- 配對：同型 left 中，以 3×|ΔMB_CA|＋5×|ΔMBON|＋|ΔKC|＋位置差 選最佳。
- **CB3198**：CB3198_R=`720575940633320785`, CB3198_L=`720575940614928234`
```
 side            root_id  mb_ca_syn  mbon_in_syn  kc_bi_syn nt_type
 left 720575940614928234         75            0         14     ACH
right 720575940633320785         89            2         10     ACH
```
- **CB2122**：CB2122_R=`720575940628347654`, CB2122_L=`720575940631305004`
```
 side            root_id  mb_ca_syn  mbon_in_syn  kc_bi_syn nt_type
 left 720575940631305004          0            2         25     ACH
right 720575940628347654         63            5         26     ACH
```
- **CB0339**：CB0339_R=`720575940632389028`, CB0339_L=`720575940632818221`
```
 side            root_id  mb_ca_syn  mbon_in_syn  kc_bi_syn nt_type
 left 720575940632818221         71            0         26     ACH
right 720575940632389028         27            4          4     ACH
```
白話：左右對稱性看表。若左腦 MB_CA 或 MBON 明顯弱於右種子，不支持「每側一顆且同質」的 DAL 故事。

## 2. 時鐘（Yin 2025 動機；Codex 名稱已核）
- ≥1：直接組合=0；2跳內=623
- ≥5：直接組合=0；2跳內=11
- **直接邊：空的**
白話：細胞內時鐘基因≠一定有很多中樞時鐘突觸。

## 3. 巧合偵測（Xia 2005 動機）
```
focus_label  syn_odor_PN_plus_KC  syn_reinf_PAM_plus_PPL1  frac_odor  frac_reinf  has_both_odor_and_reinf
   CB3198_R                   72                        0   0.143713    0.000000                        0
   CB3198_L                   82                        0   0.157692    0.000000                        0
   CB2122_R                  142                        6   0.108979    0.004605                        1
   CB2122_L                   44                        8   0.040516    0.007366                        1
   CB0339_R                  280                        0   0.312500    0.000000                        0
   CB0339_L                  254                        0   0.321519    0.000000                        0
```

## 4. MB 迴路／MBON（Wu 2017 第三種＝γ3,γ3β′1，GRASP 陰性；Codex MBON09 有、MBON08=0）
```
focus_label  syn_KC_to_focus  syn_focus_to_KC  n_KC_reciprocal_same_id  has_KC_in_and_out
   CB3198_R                7                3                        1                  1
   CB3198_L                9                5                        0                  1
   CB2122_R               10               16                        0                  1
   CB2122_L               14               11                        3                  1
   CB0339_R                2                2                        0                  1
   CB0339_L                4               22                        1                  1
```
```
focus_label           mbon_set  direct_syn_CRE_SMP_SIP  direct_syn_any_neuropil  n_twohop_paths  best_twohop_syn_path
   CB3198_R      MBON14_alpha3                       0                        0              22                    35
   CB3198_R     MBON03_beta2mp                       2                        2              32                   244
   CB3198_R MBON09_gamma3beta1                       0                        0              58                   164
   CB3198_L      MBON14_alpha3                       0                        0              22                    30
   CB3198_L     MBON03_beta2mp                       0                        0              50                   217
   CB3198_L MBON09_gamma3beta1                       0                        0              76                   136
   CB2122_R      MBON14_alpha3                       0                        0              65                   202
   CB2122_R     MBON03_beta2mp                       5                        5              78                   203
   CB2122_R MBON09_gamma3beta1                       0                        0             138                   167
   CB2122_L      MBON14_alpha3                       2                        2              88                   122
   CB2122_L     MBON03_beta2mp                       0                        0              76                   203
   CB2122_L MBON09_gamma3beta1                       0                        0             129                   140
   CB0339_R      MBON14_alpha3                       0                        0              22                    35
   CB0339_R     MBON03_beta2mp                       4                        4              29                   274
   CB0339_R MBON09_gamma3beta1                       0                        0              48                   166
   CB0339_L      MBON14_alpha3                       0                        0              23                    26
   CB0339_L     MBON03_beta2mp                       0                        0              25                   245
   CB0339_L MBON09_gamma3beta1                       0                        0              34                    86
```

## 5. 睡眠路徑（無字面 vFB→用 FB1/FB2；ExR 對照）
```
focus_label     source_set  syn_source_to_MBON03  syn_MBON03_to_focus  chain_direct_source_MBON03_focus  n_MBON03_bridging_both  syn_MBON03_to_source_reverse_lit
   CB3198_R FB1_FB2_as_vFB                     3                    2                                 1                       0                                 9
   CB3198_R    ExR_related                     1                    2                                 1                       0                                 8
   CB3198_L FB1_FB2_as_vFB                     3                    0                                 0                       0                                 9
   CB3198_L    ExR_related                     1                    0                                 0                       0                                 8
   CB2122_R FB1_FB2_as_vFB                     3                    5                                 1                       0                                 9
   CB2122_R    ExR_related                     1                    5                                 1                       0                                 8
   CB2122_L FB1_FB2_as_vFB                     3                    0                                 0                       0                                 9
   CB2122_L    ExR_related                     1                    0                                 0                       0                                 8
   CB0339_R FB1_FB2_as_vFB                     3                    4                                 1                       0                                 9
   CB0339_R    ExR_related                     1                    4                                 1                       0                                 8
   CB0339_L FB1_FB2_as_vFB                     3                    0                                 0                       0                                 9
   CB0339_L    ExR_related                     1                    0                                 0                       0                                 8
```

## 總結：做了什麼／有結果／空的

### 做了什麼
1. 左腦配對與左右比較圖
2. 時鐘六類＋5th-LNv 直接／2跳（≥1、≥5）
3. 上游氣味 vs 獎懲占比
4. KC↔候選；MBON14/03/09→候選
5. FB1/2（與 ExR）→MBON03→候選

### 有結果
- 時鐘 ≥1 有 2 跳訊號組合數=623；直接邊組合數=0
- 巧合同時氣味+獎懲（≥1）細胞數=2
- KC 雙向（≥1）細胞數=6
- FB1/2→MBON03→候選直接鏈（≥1）組合數=6

### 空的或弱的
- 時鐘**直接邊：空的**
- CB2122 左腦最佳配對 **MB_CA=0**（相對右種子弱），僅結構配對候選。
- 字面 `vFB` primary_type：空（改用 FB1*/FB2* 操作定義）。
- MBON08：0；Wu 2017 第三種以 γ3,γ3β′1 → Codex **MBON09**。
- ≥5 門檻常砍掉弱邊。
- NBLAST／濕實驗／NT 驗證：未做。

### 文獻
- Yin 2025 https://doi.org/10.3389/fnins.2025.1623251
- Wu 2017 https://doi.org/10.1038/s41598-017-07600-2
- Xia 2005：僅動機
