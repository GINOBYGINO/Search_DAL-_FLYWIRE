# 第 3 步驗證：已知結果

- **日期**：2026-10-19
- **門檻**：syn_count ≥ 5
- **資料**：`data/download/connections_princeton.csv.gz` 等
- **總評**：通過

---

## 1. KCab-p

- primary_type == `KCab-p` 顆數：**128**
- 輸入突觸落在萼部 `MB_CA_L/MB_CA_R` 的比例：**0.130**
- 全部 `MB_*` neuropil 輸入占比：**0.695**
- 萼部 neuropil 進入輸入 top 的名次：rank=3
- 判定：通過

輸入 neuropil top：

```
neuropil  syn_count  fraction
 MB_ML_L       2727  0.212964
 MB_ML_R       1522  0.118860
 MB_CA_R       1351  0.105506
MB_PED_R       1305  0.101913
   SLP_L       1045  0.081609
   SCL_L       1035  0.080828
MB_PED_L        995  0.077704
   SLP_R        967  0.075517
   PLP_L        385  0.030066
   SCL_R        374  0.029207
```

## 2. APL → Kenyon cell

- APL 顆數：**2**
- 門檻以上連到的 KC 顆數（全下游，非僅 top50）：**5111**
- KC 佔 APL 全部下游突觸量比例：**0.896**
- 判定：通過

APL 下游 top10（單一最強邊，常為 DPM／MBON）：

```
        partner_id  syn_count primary_type class  side
720575940620975696       1847          DPM  MBIN  left
720575940625934094       1137          DPM  MBIN right
720575940630496374        639       MBON05  MBON right
720575940621164720        584       MBON05  MBON  left
720575940638028607        222       MBON06  MBON right
720575940616463477        221       MBON22  MBON  left
720575940624117245        208       MBON01  MBON right
720575940628896168        192        LHMB1   NaN  left
720575940635063135        188       MBON22  MBON right
720575940630770042        188      DM1_lPN  ALPN  left
```

## 3. MBON14／MBON03 下游

- MBON14 顆數：**4**；MBON03 顆數：**2**
- 判定：通過（有下游名單；請與 Codex 目視對照）

### MBON14 下游 top15

```
        partner_id  syn_count primary_type class  side
720575940632998028        610       CB0687   NaN right
720575940622284947        605       CB0687   NaN  left
720575940632457514        394       SIP088   NaN right
720575940630668306        213       SIP088   NaN  left
720575940626385182        154       CB3112   NaN right
720575940639917400        104       CB3369   NaN right
720575940622358604        104      LHAD1k1   NaN  left
720575940641237645        100      LHAD1k1   NaN right
720575940638028607         92       MBON06  MBON right
720575940636148325         83       SLP391   NaN right
720575940627680021         79      LHPV2h1   NaN  left
720575940635870831         73       SMP026   NaN  left
720575940646488430         73       SLP391   NaN  left
720575940626220925         71       CB3369   NaN  left
720575940613809071         70       SMP026   NaN right
```

### MBON03 下游 top15

```
        partner_id  syn_count primary_type class  side
720575940606664894        503       CRE042   NaN right
720575940623387786        300       SMP177   NaN  left
720575940619134789        289       SMP177   NaN right
720575940624449318        257       SIP087   NaN  left
720575940621684358        249       SIP087   NaN right
720575940626315010        229       MBON13  MBON right
720575940621828443        227       MBON12  MBON  left
720575940630075703        220       MBON12  MBON  left
720575940620252075        219      LHPV4m1   NaN right
720575940623905719        201       MBON12  MBON right
720575940623001480        198       MBON12  MBON right
720575940645304430        193       MBON13  MBON  left
720575940618551481        160      LHPV4m1   NaN  left
720575940621033808        148       CRE042   NaN  left
720575940622999688        128      LHPV5e1   NaN right
```

## 備註

- Codex 萼部 neuropil 代碼為 `MB_CA_L`／`MB_CA_R`。
- KC 在葉部（ML／PED 等）也收調控輸入，故萼部占比不一定 >50%。
- MBON 下游請到 Codex 同一細胞類型頁核對前幾名是否一致。
