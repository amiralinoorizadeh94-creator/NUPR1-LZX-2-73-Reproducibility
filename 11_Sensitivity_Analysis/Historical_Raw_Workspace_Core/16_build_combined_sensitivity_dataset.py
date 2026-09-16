import csv
from pathlib import Path

root = Path("sensitivity_docking")

historical_file = root / "historical_20260911_structural_mode1.csv"
sensitivity_file = root / "structural_sensitivity_mode1.csv"
output_file = root / "combined_288_structural_sensitivity.csv"

rows = []

with historical_file.open(newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        r["run_type"] = "historical_original"
        rows.append(r)

with sensitivity_file.open(newline="", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        r["run_type"] = "sensitivity"
        rows.append(r)

if len(rows) != 288:
    raise RuntimeError(
        f"Expected 288 total models, found {len(rows)}"
    )

fields = [
    "run_type",
    "seed",
    "model",
    "vina_kcal_mol",
    "slahs_min_A",
    "lvtkl_min_A",
    "slahs_contact_le4A",
    "lvtkl_contact_le4A",
    "dual_contact_le4A",
]

with output_file.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print("Combined dataset created")
print("Total rows =", len(rows))
print("Historical =", sum(r["run_type"] == "historical_original" for r in rows))
print("Sensitivity =", sum(r["run_type"] == "sensitivity" for r in rows))
print("Dual =", sum(int(r["dual_contact_le4A"]) for r in rows))
print("LVTKL contact =", sum(int(r["lvtkl_contact_le4A"]) for r in rows))
print("SLAHS contact =", sum(int(r["slahs_contact_le4A"]) for r in rows))
print("Output =", output_file)
