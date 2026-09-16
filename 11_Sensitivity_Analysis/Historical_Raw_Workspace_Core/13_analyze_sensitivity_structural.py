from pathlib import Path
import csv
import math
from collections import defaultdict

SEEDS = [
    20260921,
    20260922,
    20260923,
    20260924,
    20260925,
]

RECON_ROOT = Path("sensitivity_reconstructions")
DOCK_ROOT = Path("sensitivity_docking")

SCORE_FILE = DOCK_ROOT / "all_sensitivity_scores.csv"
OUT_FILE = DOCK_ROOT / "structural_sensitivity_mode1.csv"

CONTACT_CUTOFF = 4.0


def atom_xyz(line):
    try:
        return (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        )
    except Exception:
        return None


def is_hydrogen(line):
    """
    Conservative PDB/PDBQT hydrogen identification.
    First prefer the atom name field.
    """
    atom_name = line[12:16].strip().upper()

    # Strip leading digits, e.g. 1H, 2H
    atom_name_no_digits = atom_name.lstrip("0123456789")

    return atom_name_no_digits.startswith("H")


def read_receptor_hotspots(path):
    """
    Combined receptor contains SLAHS followed by LVTKL.

    We classify peptide atoms using residue sequence order,
    based on residue number + chain identifier encountered
    in the receptor PDBQT.

    First five unique residues = SLAHS
    Next five unique residues = LVTKL
    """

    residues = []
    residue_atoms = defaultdict(list)

    for line in Path(path).read_text(
        errors="ignore"
    ).splitlines():

        if not line.startswith(("ATOM", "HETATM")):
            continue

        if is_hydrogen(line):
            continue

        xyz = atom_xyz(line)

        if xyz is None:
            continue

        chain = line[21:22]
        resnum = line[22:26].strip()
        resname = line[17:20].strip()

        key = (chain, resnum, resname)

        if key not in residue_atoms:
            residues.append(key)

        residue_atoms[key].append(xyz)

    if len(residues) != 10:
        raise RuntimeError(
            f"{path}: expected 10 peptide residues, "
            f"found {len(residues)}: {residues}"
        )

    slahs_res = residues[:5]
    lvtkl_res = residues[5:]

    slahs_atoms = [
        xyz
        for res in slahs_res
        for xyz in residue_atoms[res]
    ]

    lvtkl_atoms = [
        xyz
        for res in lvtkl_res
        for xyz in residue_atoms[res]
    ]

    return slahs_atoms, lvtkl_atoms


def read_mode1_ligand(path):
    """
    Read only MODEL 1 from Vina output.
    If MODEL records are absent, read atoms until ENDMDL/end.
    Heavy atoms only.
    """

    atoms = []

    in_model1 = False
    saw_model = False

    for line in Path(path).read_text(
        errors="ignore"
    ).splitlines():

        if line.startswith("MODEL"):

            saw_model = True

            parts = line.split()

            if len(parts) >= 2 and parts[1] == "1":
                in_model1 = True
                continue

            if in_model1:
                break

            continue

        if line.startswith("ENDMDL"):

            if in_model1:
                break

            continue

        if saw_model and not in_model1:
            continue

        if not line.startswith(("ATOM", "HETATM")):
            continue

        if is_hydrogen(line):
            continue

        xyz = atom_xyz(line)

        if xyz is not None:
            atoms.append(xyz)

    if not atoms:
        raise RuntimeError(
            f"No ligand heavy atoms parsed from mode 1: {path}"
        )

    return atoms


def min_distance(group_a, group_b):

    dmin = float("inf")

    for ax, ay, az in group_a:
        for bx, by, bz in group_b:

            dx = ax - bx
            dy = ay - by
            dz = az - bz

            d = math.sqrt(
                dx*dx + dy*dy + dz*dz
            )

            if d < dmin:
                dmin = d

    return dmin


scores = {}

with SCORE_FILE.open(
    newline="",
    encoding="utf-8"
) as f:

    reader = csv.DictReader(f)

    for row in reader:

        key = (
            int(row["seed"]),
            row["model"]
        )

        scores[key] = float(
            row["best_vina_kcal_mol"]
        )


rows = []


for seed in SEEDS:

    print()
    print("=" * 70)
    print(f"STRUCTURAL ANALYSIS SEED {seed}")
    print("=" * 70)

    for i in range(1, 49):

        model = f"COMBO_{i:03d}"

        receptor = (
            RECON_ROOT
            / f"seed_{seed}"
            / f"{model}.pdbqt"
        )

        docked = (
            DOCK_ROOT
            / f"seed_{seed}"
            / f"{model}_out.pdbqt"
        )

        if not receptor.exists():
            raise FileNotFoundError(receptor)

        if not docked.exists():
            raise FileNotFoundError(docked)

        slahs, lvtkl = read_receptor_hotspots(
            receptor
        )

        ligand = read_mode1_ligand(
            docked
        )

        d_s = min_distance(
            ligand,
            slahs
        )

        d_l = min_distance(
            ligand,
            lvtkl
        )

        slahs_contact = (
            d_s <= CONTACT_CUTOFF
        )

        lvtkl_contact = (
            d_l <= CONTACT_CUTOFF
        )

        dual = (
            slahs_contact
            and lvtkl_contact
        )

        score = scores[
            (seed, model)
        ]

        rows.append(
            {
                "seed": seed,
                "model": model,
                "vina_kcal_mol": score,
                "slahs_min_A": round(d_s, 4),
                "lvtkl_min_A": round(d_l, 4),
                "slahs_contact_le4A": int(
                    slahs_contact
                ),
                "lvtkl_contact_le4A": int(
                    lvtkl_contact
                ),
                "dual_contact_le4A": int(
                    dual
                ),
            }
        )

        print(
            f"{model} "
            f"score={score:7.3f} "
            f"SLAHS={d_s:6.3f} Å "
            f"LVTKL={d_l:6.3f} Å "
            f"dual={dual}"
        )


with OUT_FILE.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    fields = [
        "seed",
        "model",
        "vina_kcal_mol",
        "slahs_min_A",
        "lvtkl_min_A",
        "slahs_contact_le4A",
        "lvtkl_contact_le4A",
        "dual_contact_le4A",
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields
    )

    writer.writeheader()
    writer.writerows(rows)


print()
print("=" * 70)
print("STRUCTURAL SENSITIVITY COMPLETE")
print("=" * 70)

print("Total models =", len(rows))

dual_all = [
    r for r in rows
    if r["dual_contact_le4A"] == 1
]

print(
    "Dual-contact models =",
    len(dual_all),
    "/",
    len(rows),
    f"({100*len(dual_all)/len(rows):.1f}%)"
)


for seed in SEEDS:

    sub = [
        r for r in rows
        if r["seed"] == seed
    ]

    dual = [
        r for r in sub
        if r["dual_contact_le4A"] == 1
    ]

    print(
        f"Seed {seed}: "
        f"{len(dual)}/48 dual "
        f"({100*len(dual)/48:.1f}%)"
    )


if dual_all:

    ranked = sorted(
        dual_all,
        key=lambda x: x["vina_kcal_mol"]
    )

    print()
    print("TOP 10 DUAL-CONTACT MODELS")

    for r in ranked[:10]:

        print(
            f'{r["seed"]} '
            f'{r["model"]} '
            f'score={r["vina_kcal_mol"]:.3f} '
            f'SLAHS={r["slahs_min_A"]:.3f} '
            f'LVTKL={r["lvtkl_min_A"]:.3f}'
        )


print()
print("Output:")
print(OUT_FILE)
