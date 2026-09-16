from pathlib import Path
import csv, statistics
from collections import defaultdict

meta={}
with open("refined_receptors/refinement_metadata.csv",encoding="utf-8") as f:
    for r in csv.DictReader(f):
        meta[r["model"]]=r

rows=[]
with open("refinement_docking_out/summary.csv",encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if not r["best_vina_kcal_mol"]: continue
        x={"model":r["model"],"score":float(r["best_vina_kcal_mol"])}
        x.update(meta.get(r["model"],{}))
        rows.append(x)

rows.sort(key=lambda r:r["score"])
print("TOP 15")
for r in rows[:15]:
    print(r["model"],
          "score=",r["score"],
          "seed=",r.get("seed"),
          "angle=",round(float(r.get("angle_deg",0)),2),
          "trans=",round(float(r.get("translation_A",0)),2),
          "min_inter=",round(float(r.get("min_interfragment_A",0)),2))

by=defaultdict(list)
for r in rows: by[r.get("seed","?")].append(r["score"])

print("\nBY SEED")
for seed,vals in sorted(by.items()):
    print(seed,
          "n=",len(vals),
          "best=",min(vals),
          "median=",round(statistics.median(vals),3),
          "mean=",round(statistics.mean(vals),3))

print("\nINTERPRETATION TARGETS")
print("Stage 1 SLAHS median = -3.643")
print("Stage 2A LVTKL median = -3.995")
print("Stage 2B random two-hotspot median = -4.700; best = -5.748")
print("Published two-fragment reference ≈ -8.2 (not expected to reproduce exactly without original coordinates/flexibility)")
