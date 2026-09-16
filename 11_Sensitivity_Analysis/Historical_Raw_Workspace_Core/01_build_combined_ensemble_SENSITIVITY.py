from pathlib import Path
import math
import random
import numpy as np
import csv

SEEDS = [
    20260921,
    20260922,
    20260923,
    20260924,
    20260925,
]

S_DIR = Path("slahs_pdbqt")
L_DIR = Path("lvtkl_pdbqt")
ROOT_OUT = Path("sensitivity_reconstructions")
ROOT_OUT.mkdir(exist_ok=True)


def parse_atoms(path):
    atoms = []
    for line in Path(path).read_text(
        encoding="utf-8",
        errors="ignore"
    ).splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                xyz = np.array(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ],
                    dtype=float,
                )
            except Exception:
                continue

            atoms.append([line, xyz])

    if not atoms:
        raise RuntimeError(f"No atoms parsed from {path}")

    return atoms


def random_rotation():
    u1, u2, u3 = np.random.rand(3)

    q = np.array(
        [
            math.sqrt(1 - u1) * math.sin(2 * math.pi * u2),
            math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
            math.sqrt(u1) * math.sin(2 * math.pi * u3),
            math.sqrt(u1) * math.cos(2 * math.pi * u3),
        ]
    )

    x, y, z, w = q

    return np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ]
    )


def transform(atoms, R, shift):
    coords = np.array([a[1] for a in atoms])
    center = coords.mean(axis=0)
    new = (coords - center) @ R.T + shift
    return new


def min_dist(A, B):
    d = A[:, None, :] - B[None, :, :]
    return np.sqrt((d * d).sum(axis=2)).min()


def rewrite(line, xyz, serial, chain):
    s = list(line)

    while len(s) < 80:
        s.append(" ")

    serial_txt = f"{serial:5d}"
    s[6:11] = list(serial_txt)
    s[21] = chain

    x, y, z = xyz

    s[30:38] = list(f"{x:8.3f}")
    s[38:46] = list(f"{y:8.3f}")
    s[46:54] = list(f"{z:8.3f}")

    return "".join(s).rstrip()


sfiles = sorted(S_DIR.glob("SLAHS_*.pdbqt"))
lfiles = sorted(L_DIR.glob("LVTKL_*.pdbqt"))

if len(sfiles) != 12:
    raise RuntimeError(f"Expected 12 SLAHS files, found {len(sfiles)}")

if len(lfiles) != 12:
    raise RuntimeError(f"Expected 12 LVTKL files, found {len(lfiles)}")


all_summary = []

for seed in SEEDS:

    random.seed(seed)
    np.random.seed(seed)

    out = ROOT_OUT / f"seed_{seed}"
    out.mkdir(exist_ok=True)

    target_models = 48
    attempts = 0
    made = 0
    meta = []

    print()
    print("=" * 70)
    print(f"SEED {seed}")
    print("=" * 70)

    while made < target_models and attempts < 10000:

        attempts += 1

        sf = random.choice(sfiles)
        lf = random.choice(lfiles)

        sa = parse_atoms(sf)
        la = parse_atoms(lf)

        sep = random.uniform(8.0, 14.0)

        axis = np.random.normal(size=3)
        axis = axis / np.linalg.norm(axis)

        Rs = random_rotation()
        Rl = random_rotation()

        Sc = transform(sa, Rs, -axis * sep / 2)
        Lc = transform(la, Rl, axis * sep / 2)

        dmin = min_dist(Sc, Lc)

        if not (2.2 <= dmin <= 7.0):
            continue

        allc = np.vstack([Sc, Lc])
        span = allc.max(axis=0) - allc.min(axis=0)

        if span.max() > 26:
            continue

        made += 1

        outfile = out / f"COMBO_{made:03d}.pdbqt"

        lines = []
        serial = 1

        for (orig, _), xyz in zip(sa, Sc):
            lines.append(rewrite(orig, xyz, serial, "A"))
            serial += 1

        for (orig, _), xyz in zip(la, Lc):
            lines.append(rewrite(orig, xyz, serial, "B"))
            serial += 1

        outfile.write_text(
            "\n".join(lines) + "\n",
            encoding="utf-8"
        )

        meta.append(
            (
                outfile.name,
                sf.name,
                lf.name,
                sep,
                dmin,
                *span,
            )
        )

        print(
            f"{outfile.name}: "
            f"{sf.name} + {lf.name} "
            f"sep={sep:.3f} "
            f"min_inter={dmin:.3f}"
        )

    with open(
        out / "ensemble_metadata.csv",
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        w = csv.writer(f)

        w.writerow(
            [
                "model",
                "slahs",
                "lvtkl",
                "centroid_sep_A",
                "min_interfragment_A",
                "span_x",
                "span_y",
                "span_z",
            ]
        )

        w.writerows(meta)

    if made != 48:
        raise RuntimeError(
            f"Seed {seed} produced only {made} accepted models"
        )

    all_summary.append(
        [
            seed,
            made,
            attempts,
        ]
    )

    print(
        f"\nSeed {seed}: "
        f"built {made} accepted models "
        f"after {attempts} attempts."
    )


with open(
    ROOT_OUT / "seed_generation_summary.csv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    w = csv.writer(f)
    w.writerow(
        [
            "seed",
            "accepted_models",
            "attempts",
        ]
    )

    w.writerows(all_summary)


print()
print("=" * 70)
print("Sensitivity reconstruction complete")
print("=" * 70)

for row in all_summary:
    print(
        f"seed={row[0]} "
        f"accepted={row[1]} "
        f"attempts={row[2]}"
    )

print(
    f"\nTotal new models: "
    f"{sum(r[1] for r in all_summary)}"
)