from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "_probe_types.txt"
ct = pd.read_csv(ROOT / "data/download/consolidated_cell_types.csv.gz")
cl = pd.read_csv(ROOT / "data/download/classification.csv.gz")
lines = [f"ct={len(ct)} cl={len(cl)}"]

def dump(title, vals):
    lines.append(f"## {title} n={len(vals)}")
    for v in vals:
        n = int((ct["primary_type"] == v).sum())
        lines.append(f"  {v}\t{n}")

for key in ["LNv", "LNd", "DN1", "DN2", "DN3", "PDF", "s-LN", "l-LN", "CPDN", "LPN"]:
    m = ct["primary_type"].astype(str).str.contains(key, case=False, na=False)
    dump(key, sorted(ct.loc[m, "primary_type"].astype(str).unique()))

for key in ["vFB", "^FB", "ExR", "SMP", "sleep"]:
    m = ct["primary_type"].astype(str).str.contains(key, case=False, na=False, regex=True)
    dump("vFBish:"+key, sorted(ct.loc[m, "primary_type"].astype(str).unique())[:80])

m = ct["primary_type"].astype(str).str.startswith("MBON", na=False)
dump("MBON", sorted(ct.loc[m, "primary_type"].astype(str).unique()))

# classification hints for clock
for col in ["class", "sub_class", "super_class"]:
    if col in cl.columns:
        m = cl[col].astype(str).str.contains("clock|circadian|LNv|LNd|DN1", case=False, na=False)
        vals = sorted(cl.loc[m, col].astype(str).unique())
        lines.append(f"## cl.{col} clockish n={len(vals)}")
        lines.extend([f"  {v}" for v in vals[:100]])

for t in ["CB3198", "CB2122", "CB0339"]:
    sub = ct[ct["primary_type"] == t].merge(cl[["root_id", "side"]], on="root_id", how="left")
    lines.append(f"## {t} n={len(sub)}")
    for r in sub.itertuples(index=False):
        lines.append(f"  {r.root_id}\t{r.side}")

# labels file may have more names
lab = ROOT / "data/download/labels.csv.gz"
if lab.exists():
    L = pd.read_csv(lab)
    lines.append(f"## labels cols={list(L.columns)}")
    for key in ["LNv", "LNd", "DN1a", "DN1p", "DN2", "DN3", "vFB", "MBON08", "MBON09"]:
        hit = []
        for c in L.columns:
            if L[c].dtype == object:
                m = L[c].astype(str).str.contains(key, case=False, na=False)
                if m.any():
                    hit.append((c, int(m.sum()), sorted(L.loc[m, c].astype(str).unique())[:20]))
        lines.append(f"labels {key}: {hit[:5]}")

OUT.write_text("\n".join(lines), encoding="utf-8")
print(OUT)
