# DAL2 全庫篩選（僅 PPL105＋PPL106 雙輸出）20261024

- **資料**：FlyWire v783 `connections_princeton_no_threshold.csv.gz`
- **硬條件**：同時輸出到 PPL105（α2α′2）與 PPL106（α3）；各達門檻
- **排除**：Kenyon_Cell／MBON／DAN／ALPN／LHCENT／visual_projection／APL／DPM／LHMB1／所有 PPL*
- **soma 粗框**：|x−512k|≥60k、y∈[110k,240k]、z≥100k（僅標記，不刪除）
- **聲明**：連結性候選 ≠ 找到 DAL2；NPF 無法用 Codex NT 驗證。

## 門檻 ≥5（對 PPL105、PPL106 各自加總）

- **通過數**：22
- **側別**：left=13, right=9, other/NA=0
- **背前外側粗框內**：13
- **SMP194 排名**：9（syn_sum=48）

### 關注細胞

```
 rank   type            root_id  side  syn_to_PPL105  syn_to_PPL106  syn_sum_both predicted_NT      x      y      z  in_DAL_box  watch class super_class
    9 SMP194 720575940642272456 right             22             26            48          ACH 650236 148992 172880        True SMP194   NaN     central
```

### 前 20 名

```
 rank    type            root_id  side  syn_to_PPL105  syn_to_PPL106  syn_sum_both predicted_NT      x      y      z  in_DAL_box  watch class super_class
    1  SMP116 720575940629852419 right              7            101           108         GLUT 533036 111732  45120       False          NaN     central
    2  CB1226 720575940637762906  left             91             12           103         GLUT 509064  96656  93440       False          NaN     central
    3  CB1489 720575940633516002  left             75             14            89          ACH 379936 131648 206320        True          NaN     central
    4  CB1226 720575940626259458 right             73              9            82         GLUT 522588 121372  53080       False          NaN     central
    5  CB1226 720575940626495109  left             71              5            76         GLUT 512948 113736  57200       False          NaN     central
    6  CB1489 720575940620054430 right             64              6            70          ACH 651080 149292 172920        True          NaN     central
    7  CB1489 720575940643929111 right             47             15            62          ACH 651648 150148 172920        True          NaN     central
    8  CB2444 720575940615223347 right             43             19            62          ACH 614168 123612  95360       False          NaN     central
    9  SMP194 720575940642272456 right             22             26            48          ACH 650236 148992 172880        True SMP194   NaN     central
   10  CB2444 720575940612090866  left             30             12            42          ACH 401536 141440 204280        True          NaN     central
   11  SMP194 720575940643214628  left             24             16            40          ACH 442536 114056 112920        True          NaN     central
   12  SMP194 720575940637439566 right              6             34            40          ACH 647460 151340 172920        True          NaN     central
   13  CB1489 720575940636833392  left             29             10            39          ACH 392724 164712 211760        True          NaN     central
   14 SIP028b 720575940638824245  left             23              6            29         GABA 404732 202968  53520       False          NaN     central
   15  SMP194 720575940624207145  left              5             24            29          ACH 438192 115264  86680       False          NaN     central
   16  CB1727 720575940632278189  left              9             14            23          ACH 395456 129696 202560        True          NaN     central
   17  SIP057 720575940610593234  left             16              6            22          ACH 439800 113852 115120        True          NaN     central
   18 SIP028a 720575940622227453  left             16              5            21         GABA 484360 171396  91560       False          NaN     central
   19  CB2310 720575940614387198  left              6             12            18          ACH 366368 146016 210600        True          NaN     central
   20  CB1357 720575940616596822  left              9              8            17          ACH 390752 132608 204600        True          NaN     central
```

## 門檻 ≥1（對 PPL105、PPL106 各自加總）

- **通過數**：155
- **側別**：left=83, right=72, other/NA=0
- **背前外側粗框內**：58
- **SMP194 排名**：13（syn_sum=48）

### 關注細胞

```
 rank   type            root_id  side  syn_to_PPL105  syn_to_PPL106  syn_sum_both predicted_NT      x      y      z  in_DAL_box  watch class super_class
   13 SMP194 720575940642272456 right             22             26            48          ACH 650236 148992 172880        True SMP194   NaN     central
```

### 前 20 名

```
 rank   type            root_id  side  syn_to_PPL105  syn_to_PPL106  syn_sum_both predicted_NT      x      y      z  in_DAL_box  watch class super_class
    1 SMP116 720575940629852419 right              7            101           108         GLUT 533036 111732  45120       False          NaN     central
    2 CB1226 720575940637762906  left             91             12           103         GLUT 509064  96656  93440       False          NaN     central
    3 CB1489 720575940633516002  left             75             14            89          ACH 379936 131648 206320        True          NaN     central
    4 SMP116 720575940615763573  left              1             85            86         GLUT 503488  97216  55400       False          NaN     central
    5 CB1226 720575940626259458 right             73              9            82         GLUT 522588 121372  53080       False          NaN     central
    6 CB1226 720575940626495109  left             71              5            76         GLUT 512948 113736  57200       False          NaN     central
    7 CB1489 720575940620054430 right             64              6            70          ACH 651080 149292 172920        True          NaN     central
    8 SLP073 720575940605675808 right              1             62            63          ACH 703712 194688  81800       False          NaN     central
    9 CB1489 720575940643929111 right             47             15            62          ACH 651648 150148 172920        True          NaN     central
   10 CB2444 720575940615223347 right             43             19            62          ACH 614168 123612  95360       False          NaN     central
   11 SMP026 720575940635870831  left              3             47            50          ACH 363324 153404 104680        True          NaN     central
   12 SMP026 720575940613809071 right              2             48            50          ACH 662488 165664  46880       False          NaN     central
   13 SMP194 720575940642272456 right             22             26            48          ACH 650236 148992 172880        True SMP194   NaN     central
   14 CB2444 720575940612090866  left             30             12            42          ACH 401536 141440 204280        True          NaN     central
   15 CB1226 720575940642653147 right             39              2            41         GLUT 529376  99692  84320       False          NaN     central
   16 SMP053 720575940611890289  left              2             39            41          ACH 514464 142752  34200       False          NaN     central
   17 SMP194 720575940643214628  left             24             16            40          ACH 442536 114056 112920        True          NaN     central
   18 SMP194 720575940637439566 right              6             34            40          ACH 647460 151340 172920        True          NaN     central
   19 CB1489 720575940636833392  left             29             10            39          ACH 392724 164712 211760        True          NaN     central
   20 CB1871 720575940607802123  left              2             35            37         GLUT 514976  94528  62280       False          NaN     central
```
