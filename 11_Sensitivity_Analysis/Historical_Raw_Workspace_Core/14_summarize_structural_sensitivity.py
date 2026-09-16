import csv
import statistics
from pathlib import Path
from collections import defaultdict

p = Path("sensitivity_docking/structural_sensitivity_mode1.csv")

with p.open(newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

for r in rows:
    r["seed"] = int(r["seed"])
    r["score"] = float(r["vina_kcal_mol"])
    r["slahs"] = float(r["slahs_min_A"])
    r["lvtkl"] = float(r["lvtkl_min_A"])
    r["dual"] = int(r["dual_contact_le4A"])

print("=" * 70)
print("PRIMARY STRUCTURAL ROBUSTNESS SUMMARY")
print("=" * 70)

print("Total =", len(rows))
print("Dual =", sum(r["dual"] for r in rows))

slahs_contact = sum(r["slahs"] <= 4.0 for r in rows)
lvtkl_contact = sum(r["lvtkl"] <= 4.0 for r in rows)

print(
    f"SLAHS contact <=4 A = "
    f"{slahs_contact}/{len(rows)} "
    f"({100*slahs_contact/len(rows):.1f}%)"
)

print(
    f"LVTKL contact <=4 A = "
    f"{lvtkl_contact}/{len(rows)} "
    f"({100*lvtkl_contact/len(rows):.1f}%)"
)

print()
print("DISTANCE DISTRIBUTIONS")

for name, key in [
    ("SLAHS", "slahs"),
    ("LVTKL", "lvtkl")
]:
    vals = [r[key] for r in rows]

    print(
        f"{name}: "
        f"mean={statistics.mean(vals):.3f} "
        f"median={statistics.median(vals):.3f} "
        f"min={min(vals):.3f} "
        f"max={max(vals):.3f}"
    )

dual = [r for r in rows if r["dual"] == 1]
nondual = [r for r in rows if r["dual"] == 0]

print()
print("DOCKING SCORE BY STRUCTURAL CLASS")

for name, group in [
    ("Dual", dual),
    ("Non-dual", nondual)
]:
    vals = [r["score"] for r in group]

    print(
        f"{name}: n={len(vals)} "
        f"mean={statistics.mean(vals):.3f} "
        f"median={statistics.median(vals):.3f} "
        f"best={min(vals):.3f}"
    )

print()
print("NON-DUAL MODELS")

for r in nondual:
    print(
        r["seed"],
        r["model"],
        f'score={r["score"]:.3f}',
        f'SLAHS={r["slahs"]:.3f}',
        f'LVTKL={r["lvtkl"]:.3f}'
    )
