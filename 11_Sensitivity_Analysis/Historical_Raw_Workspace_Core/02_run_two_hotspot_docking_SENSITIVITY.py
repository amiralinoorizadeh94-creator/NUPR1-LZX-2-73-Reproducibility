from pathlib import Path
import subprocess, re, csv, statistics

vina = Path("vina.exe")
lig = Path("LZX-2-73.pdbqt")

if not vina.exists():
    raise FileNotFoundError("vina.exe missing")

if not lig.exists():
    raise FileNotFoundError("LZX-2-73.pdbqt missing")


SEEDS = [
    20260921,
    20260922,
    20260923,
    20260924,
    20260925,
]

ROOT = Path("sensitivity_reconstructions")
OUTROOT = Path("sensitivity_docking")
OUTROOT.mkdir(exist_ok=True)


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
            except:
                pass

    if not pts:
        raise RuntimeError(
            f"No coordinates parsed from {path}"
        )

    return pts


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
    print(f"DOCKING SEED {seed}")
    print("=" * 70)

    for idx, rec in enumerate(recs, 1):

        pts = coords(rec)

        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        zs = [p[2] for p in pts]

        center = (
            sum(xs) / len(xs),
            sum(ys) / len(ys),
            sum(zs) / len(zs),
        )

        cfg = outdir / (rec.stem + ".txt")

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

        out = outdir / (
            rec.stem + "_out.pdbqt"
        )

        cp = subprocess.run(
            [
                str(vina.resolve()),
                "--config",
                str(cfg.resolve()),
                "--out",
                str(out.resolve()),
            ],
            capture_output=True,
            text=True
        )

        txt = (
            (cp.stdout or "")
            + "\n"
            + (cp.stderr or "")
        )

        (
            outdir / (rec.stem + ".log")
        ).write_text(
            txt,
            encoding="utf-8"
        )

        if cp.returncode != 0:
            print(txt)
            raise SystemExit(cp.returncode)

        best = None

        for line in txt.splitlines():

            m = re.match(
                r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
                line
            )

            if m:
                best = float(m.group(1))
                break

        if best is None:
            raise RuntimeError(
                f"Could not parse best score for "
                f"seed {seed} {rec.stem}"
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

        print(
            f"[{idx:02d}/48] "
            f"{rec.stem} "
            f"best = {best:.3f}"
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

    vals = [
        v for _, _, v in rows
        if v is not None
    ]

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
        round(statistics.mean(vals), 3)
    )

    for cutoff in [-5, -6, -7, -8]:

        n = sum(
            v <= cutoff
            for v in vals
        )

        print(
            f"fraction <= {cutoff:.1f} "
            f"kcal/mol = "
            f"{n}/{len(vals)} "
            f"({100*n/len(vals):.1f}%)"
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
    v for _, _, v in master_rows
    if v is not None
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
