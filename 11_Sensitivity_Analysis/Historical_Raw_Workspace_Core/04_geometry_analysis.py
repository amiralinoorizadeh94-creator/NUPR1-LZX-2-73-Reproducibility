# 04_geometry_analysis.py
# Relate two-hotspot geometry to docking score and identify the most promising arrangements.

from pathlib import Path
import csv, statistics, math

meta_path = Path("combined_receptors/ensemble_metadata.csv")
score_path = Path("two_hotspot_docking_out/summary.csv")

if not meta_path.exists():
    raise FileNotFoundError(meta_path)
if not score_path.exists():
    raise FileNotFoundError(score_path)

meta = {}
with meta_path.open(encoding="utf-8") as f:
    for r in csv.DictReader(f):
        model = Path(r["model"]).stem
        meta[model] = {
            "slahs": r["slahs"],
            "lvtkl": r["lvtkl"],
            "centroid_sep_A": float(r["centroid_sep_A"]),
            "min_interfragment_A": float(r["min_interfragment_A"]),
            "span_x": float(r["span_x"]),
            "span_y": float(r["span_y"]),
            "span_z": float(r["span_z"]),
        }

rows = []
with score_path.open(encoding="utf-8") as f:
    for r in csv.DictReader(f):
        model = r["model"]
        score = float(r["best_vina_kcal_mol"])
        if model in meta:
            rows.append({"model": model, "score": score, **meta[model]})

rows.sort(key=lambda x: x["score"])

def mean(vals):
    return sum(vals)/len(vals) if vals else float("nan")

print("TOP 10 WITH GEOMETRY")
for r in rows[:10]:
    print(
        r["model"],
        "score=", r["score"],
        "centroid_sep_A=", round(r["centroid_sep_A"], 2),
        "min_interfragment_A=", round(r["min_interfragment_A"], 2),
        "span=", tuple(round(r[k],2) for k in ("span_x","span_y","span_z")),
        "S=", r["slahs"],
        "L=", r["lvtkl"],
    )

# Compare top quartile vs rest
scores = [r["score"] for r in rows]
q1 = statistics.quantiles(scores, n=4, method="inclusive")[0]
topq = [r for r in rows if r["score"] <= q1]
rest = [r for r in rows if r["score"] > q1]

print("\nGEOMETRY SUMMARY")
print("Q1 cutoff:", q1)
print("Top-quartile n:", len(topq), "Rest n:", len(rest))
for key in ["centroid_sep_A","min_interfragment_A","span_x","span_y","span_z"]:
    print(
        key,
        "topQ_mean=", round(mean([r[key] for r in topq]),3),
        "rest_mean=", round(mean([r[key] for r in rest]),3)
    )

# Best 5 models for refinement
print("\nREFINEMENT SEEDS")
for r in rows[:5]:
    print(r["model"], r["score"])

# Write merged CSV
out = Path("two_hotspot_docking_out/geometry_score_merged.csv")
with out.open("w", newline="", encoding="utf-8") as f:
    fields = ["model","score","slahs","lvtkl","centroid_sep_A","min_interfragment_A","span_x","span_y","span_z"]
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print("\nWrote:", out)
