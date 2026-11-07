# 九顆候選：四點對照表

資料：`connections_princeton_no_threshold`；input／output 取該細胞全部突觸依 neuropil 加總後前 5／前 10。

## 總表

| # | 類型 | 側 | 位置像DAL | 形態(SMP/SIP入＋SLP/CRE出) | 萼部MB_CA | NT身份 | root_id |
|---|---|---|---|---|---|---|---|
| 1 | LHPV5e1 | right | 像（背前外側粗框內） | 像 | 有（加分）（2） | 非DA/5HT（不排除）（ACH） | `720575940622999688` |
| 2 | LHPV5e1 | left | 偏離（不夠外側） | 像 | 有（加分）（6） | 非DA/5HT（不排除）（ACH） | `720575940624155672` |
| 3 | AL-MBDL1 | right | 像（背前外側粗框內） | 像 | 有（加分）（19） | 非DA/5HT（不排除）（NA） | `720575940624206472` |
| 4 | AL-MBDL1 | left | 偏離（偏腹／不夠外側） | 像 | 有（加分）（4） | 非DA/5HT（不排除）（NA） | `720575940648669060` |
| 5 | OA-VPM3 | right | 偏離（太前／不夠外側） | 像 | 有（加分）（504） | 非DA/5HT（不排除）（OCT） | `720575940640356149` |
| 6 | SMP194 | right | 像（背前外側粗框內） | 像 | 有（加分）（3） | 非DA/5HT（不排除）（ACH） | `720575940642272456` |
| 7 | CB2122 | right | 像（背前外側粗框內） | 像 | 有（加分）（63） | 非DA/5HT（不排除）（ACH） | `720575940628347654` |
| 8 | CB0339 | right | 像（背前外側粗框內） | 像 | 有（加分）（27） | 非DA/5HT（不排除）（ACH） | `720575940632389028` |
| 9 | CB3198 | right | 像（背前外側粗框內） | 像 | 有（加分）（89） | 非DA/5HT（不排除）（ACH） | `720575940633320785` |

### 位置備註

- DAL 文獻：每側約 1 顆。本表同型若出現左右各一（LHPV5e1、AL-MBDL1），記「不完全符合」。

### 形態細節

| # | input 前5 | SMP/SIP入 | output 前5 | SLP出 | CRE出 |
|---|---|---|---|---|---|
| 1 | SMP_R:2880, SIP_R:1153, SLP_R:895, CRE_R:785, LH_R:173 | 是：SMP_R,SIP_R,SMP_L,SIP_L | SMP_R:1451, CRE_R:1439, SIP_R:1048, CRE_L:944, SIP_L:676 | 是：SLP_R,SLP_L | 是：CRE_R,CRE_L |
| 2 | SMP_L:2200, SIP_L:903, SLP_L:741, CRE_L:690, LH_L:153 | 是：SMP_L,SIP_L,SMP_R,SIP_R | CRE_L:1191, SMP_L:1082, SMP_R:1008, CRE_R:963, SIP_R:721 | 是：SLP_L,SLP_R | 是：CRE_L,CRE_R |
| 3 | SMP_R:830, CRE_R:576, LH_R:493, SIP_R:405, AL_L:370 | 是：SMP_R,SIP_R | AL_R:3797, AL_L:3520, SMP_R:67, SIP_R:59, LH_R:53 | 是：SLP_R | 是：CRE_R |
| 4 | SMP_L:806, CRE_L:665, LH_L:590, SIP_L:482, AL_R:356 | 是：SMP_L,SIP_L | AL_R:4033, AL_L:3693, SMP_L:45, SIP_L:30, LH_L:29 | 是：SLP_L | 是：CRE_L |
| 5 | FB:624, SAD:425, SLP_R:360, GNG:226, SMP_R:151 | 是：SMP_R,SIP_R | FB:1596, SLP_R:1230, SMP_R:506, MB_CA_R:440, CRE_R:412 | 是：SLP_R | 是：CRE_R |
| 6 | SMP_R:183, CRE_R:54, SIP_R:42, MB_VL_R:3, SLP_R:3 | 是：SMP_R,SIP_R | SMP_R:258, SIP_R:56, SLP_R:11, CRE_R:9, MB_VL_R:7 | 是：SLP_R | 是：CRE_R |
| 7 | SLP_R:483, SMP_R:446, CRE_R:165, LH_R:125, MB_CA_R:35 | 是：SMP_R,SIP_R | SMP_R:340, SLP_R:243, CRE_R:176, LH_R:81, MB_CA_R:28 | 是：SLP_R | 是：CRE_R |
| 8 | LH_R:636, SMP_R:145, SIP_R:32, CRE_R:29, SLP_R:28 | 是：SMP_R,SIP_R | LH_R:379, SMP_R:347, SIP_R:56, CRE_R:29, SLP_R:28 | 是：SLP_R | 是：CRE_R |
| 9 | SMP_R:169, SIP_R:111, MB_CA_R:82, LH_R:69, SLP_R:43 | 是：SMP_R,SIP_R | SMP_R:263, SIP_R:190, CRE_R:34, LH_R:21, MB_VL_R:16 | 是：SLP_R | 是：CRE_R |

### 評分說明

- **位置像**：落在背前外側粗框（\|x−512k\|≥60k、y∈[110–240]k、z≥100k）。
- **形態像**：input 前10 含 SMP/SIP，且 output 前10 含 SLP 或 CRE。
- **萼部**：任一方向有 MB_CA_* 突觸 → 與 GRASP 萼部接觸加分。
- **NT**：Codex 預測為 DA／5HT（SER）→ DDC 相關加分；否則不排除。
