
from pathlib import Path
import csv, statistics

path=Path("two_hotspot_docking_out/summary.csv")
rows=[]
with open(path,encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["best_vina_kcal_mol"]:
            rows.append((r["model"],float(r["best_vina_kcal_mol"])))
rows.sort(key=lambda x:x[1])

print("TOP 10")
for m,v in rows[:10]:
    print(m,v)

vals=[v for _,v in rows]
print("\nROBUST SUMMARY")
print("n =",len(vals))
print("best =",min(vals))
print("median =",statistics.median(vals))
print("mean =",round(statistics.mean(vals),3))
print("Q1 =",statistics.quantiles(vals,n=4,method="inclusive")[0])
print("Q3 =",statistics.quantiles(vals,n=4,method="inclusive")[2])
