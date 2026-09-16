from pathlib import Path
import subprocess
import re
import csv
import statistics
import time

vina = Path("vina.exe")
lig = Path("LZX-2-73.pdbqt")

SEEDS = [
    20260921,
    20260922,
    20260923,
    20260924,
    20260925,
]

ROOT = Path("sensitivity_reconstructions")
OUTROOT = Path("sensitivity_docking")

MAX_RETRIES = 3


def coords(path):
    pts = []

    for line in Path(path).read_text(
        errors="ignore"
    ).splitlines():

        if line.startswith(("ATOM", "HETATM")):
            try:
                pts.append(
                    (
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
            except Exception:
                pass

    if not pts:
        raise RuntimeError(
            f"No coordinates parsed from {path}"
        )

    return pts


def parse_best(text):

    for line in text.splitlines():

        m = re.match(
            r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
            line
        )

        if m:
            return float(m.group(1))

    return None


if not vina.exists():
    raise FileNotFoundError("vina.exe missing")

if not lig.exists():
    raise FileNotFoundError(
        "LZX-2-73.pdbqt missing"
    )

OUTROOT.mkdir(exist_ok=True)

master_rows = []


for seed in SEEDS:

    receptor_dir = ROOT / f"seed_{seed}"
    outdir = OUTROOT / f"seed_{seed}"

    outdir.mkdir(exist_ok=True)

    recs = sorted(
        receptor_dir.glob("COMBO_*.pdbqt")
    )

    if len(recs) != 48:
        raise RuntimeError(
            f"Seed {seed}: expected 48 receptors, "
            f"found {len(recs)}"
        )

    rows = []

    print()
    print("=" * 70)
    print(f"DOCKING / RESUME SEED {seed}")
    print("=" * 70)

    for idx, rec in enumerate(recs, 1):

        logfile = outdir / (
            rec.stem + ".log"
        )

        outfile = outdir / (
            rec.stem + "_out.pdbqt"
        )

        existing_best = None

        if logfile.exists() and outfile.exists():

            oldtxt = logfile.read_text(
                errors="ignore"
            )

            existing_best = parse_best(
                oldtxt
            )

        if existing_best is not None:

            best = existing_best

            print(
                f"[{idx:02d}/48] "
                f"{rec.stem} "
                f"SKIP existing "
                f"best = {best:.3f}"
            )

        else:

            pts = coords(rec)

            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            zs = [p[2] for p in pts]

            center = (
                sum(xs) / len(xs),
                sum(ys) / len(ys),
                sum(zs) / len(zs),
            )

            cfg = outdir / (
                rec.stem + ".txt"
            )

            cfg.write_text(
                f"""receptor = {rec.resolve().as_posix()}
ligand = {lig.resolve().as_posix()}
center_x = {center[0]:.3f}
center_y = {center[1]:.3f}
center_z = {center[2]:.3f}
size_x = 30
size_y = 30
size_z = 30
exhaustiveness = 8
num_modes = 10
energy_range = 6
""",
                encoding="utf-8"
            )

            best = None

            for attempt in range(
                1,
                MAX_RETRIES + 1
            ):

                print(
                    f"[{idx:02d}/48] "
                    f"{rec.stem} "
                    f"attempt {attempt}/{MAX_RETRIES}"
                )

                cp = subprocess.run(
                    [
                        str(vina.resolve()),
                        "--config",
                        str(cfg.resolve()),
                        "--out",
                        str(outfile.resolve()),
                    ],
                    capture_output=True,
                    text=True
                )

                txt = (
                    (cp.stdout or "")
                    + "\n"
                    + (cp.stderr or "")
                )

                logfile.write_text(
                    txt,
                    encoding="utf-8"
                )

                if cp.returncode == 0:

                    best = parse_best(txt)

                    if best is not None:
                        break

                print(
                    f"WARNING: failed attempt "
                    f"{attempt} for {rec.stem}; "
                    f"return code={cp.returncode}"
                )

                time.sleep(1)

            if best is None:
                raise RuntimeError(
                    f"Vina failed after "
                    f"{MAX_RETRIES} attempts: "
                    f"seed {seed}, {rec.stem}"
                )

            print(
                f"          SUCCESS "
                f"best = {best:.3f}"
            )

        rows.append(
            (
                seed,
                rec.stem,
                best,
            )
        )

        master_rows.append(
            (
                seed,
                rec.stem,
                best,
            )
        )

    with open(
        outdir / "summary.csv",
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        w = csv.writer(f)

        w.writerow(
            [
                "seed",
                "model",
                "best_vina_kcal_mol",
            ]
        )

        w.writerows(rows)

    vals = [x[2] for x in rows]

    print()
    print(f"SEED {seed} SUMMARY")
    print("n =", len(vals))
    print("best =", min(vals))
    print(
        "median =",
        statistics.median(vals)
    )
    print(
        "mean =",
        round(
            statistics.mean(vals),
            3
        )
    )


with open(
    OUTROOT / "all_sensitivity_scores.csv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.writer(f)

    w.writerow(
        [
            "seed",
            "model",
            "best_vina_kcal_mol",
        ]
    )

    w.writerows(master_rows)


all_vals = [
    x[2] for x in master_rows
]

print()
print("=" * 70)
print("ALL SENSITIVITY DOCKING COMPLETE")
print("=" * 70)

print(
    "Total models =",
    len(all_vals)
)

print(
    "Overall best =",
    min(all_vals)
)

print(
    "Overall median =",
    statistics.median(all_vals)
)

print(
    "Overall mean =",
    round(
        statistics.mean(all_vals),
        3
    )
)

for cutoff in [-5, -6, -7, -8]:

    n = sum(
        v <= cutoff
        for v in all_vals
    )

    print(
        f"Overall <= {cutoff:.1f} "
        f"kcal/mol = "
        f"{n}/{len(all_vals)} "
        f"({100*n/len(all_vals):.1f}%)"
    )
