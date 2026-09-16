from pathlib import Path
import csv
import math
from collections import defaultdict

RECEPTOR_DIR = Path("combined_receptors")
DOCK_DIR = Path("two_hotspot_docking_out")

OUT_FILE = Path(
    "sensitivity_docking/"
    "historical_20260911_structural_mode1.csv"
)

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
    atom_name = line[12:16].strip().upper()
    atom_name = atom_name.lstrip("0123456789")
    return atom_name.startswith("H")


def read_receptor_hotspots(path):

    residues = []
    residue_atoms = defaultdict(list)

    for line in Path(path).read_text(
        errors="ignore"
    ).splitlines():

        if not line.startswith(
            ("ATOM", "HETATM")
        ):
            continue

        if is_hydrogen(line):
            continue

        xyz = atom_xyz(line)

        if xyz is None:
            continue

        key = (
            line[21:22],
            line[22:26].strip(),
            line[17:20].strip(),
        )

        if key not in residue_atoms:
            residues.append(key)

        residue_atoms[key].append(xyz)

    if len(residues) != 10:
        raise RuntimeError(
            f"{path}: expected 10 residues, "
            f"found {len(residues)}"
        )

    slahs = [
        xyz
        for res in residues[:5]
        for xyz in residue_atoms[res]
    ]

    lvtkl = [
        xyz
        for res in residues[5:]
        for xyz in residue_atoms[res]
    ]

    return slahs, lvtkl


def read_mode1(path):

    atoms = []
    saw_model = False
    in_model1 = False

    for line in Path(path).read_text(
        errors="ignore"
    ).splitlines():

        if line.startswith("MODEL"):

            saw_model = True
            parts = line.split()

            if (
                len(parts) >= 2
                and parts[1] == "1"
            ):
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

        if not line.startswith(
            ("ATOM", "HETATM")
        ):
            continue

        if is_hydrogen(line):
            continue

        xyz = atom_xyz(line)

        if xyz is not None:
            atoms.append(xyz)

    if not atoms:
        raise RuntimeError(
            f"No mode-1 atoms parsed: {path}"
        )

    return atoms


def min_distance(a, b):

    best = float("inf")

    for ax, ay, az in a:
        for bx, by, bz in b:

            d = math.sqrt(
                (ax-bx)**2
                + (ay-by)**2
                + (az-bz)**2
            )

            if d < best:
                best = d

    return best


scores = {}

summary = DOCK_DIR / "summary.csv"

with summary.open(
    newline="",
    encoding="utf-8"
) as f:

    for row in csv.DictReader(f):

        scores[row["model"]] = float(
            row["best_vina_kcal_mol"]
        )


rows = []


for i in range(1, 49):

    model = f"COMBO_{i:03d}"

    receptor = (
        RECEPTOR_DIR
        / f"{model}.pdbqt"
    )

    docked = (
        DOCK_DIR
        / f"{model}_out.pdbqt"
    )

    slahs, lvtkl = (
        read_receptor_hotspots(receptor)
    )

    ligand = read_mode1(docked)

    ds = min_distance(
        ligand,
        slahs
    )

    dl = min_distance(
        ligand,
        lvtkl
    )

    s_contact = ds <= CONTACT_CUTOFF
    l_contact = dl <= CONTACT_CUTOFF

    dual = (
        s_contact
        and l_contact
    )

    score = scores[model]

    rows.append(
        {
            "seed": 20260911,
            "model": model,
            "vina_kcal_mol": score,
            "slahs_min_A": round(ds, 4),
            "lvtkl_min_A": round(dl, 4),
            "slahs_contact_le4A": int(
                s_contact
            ),
            "lvtkl_contact_le4A": int(
                l_contact
            ),
            "dual_contact_le4A": int(
                dual
            ),
        }
    )

    print(
        f"{model} "
        f"score={score:7.3f} "
        f"SLAHS={ds:6.3f} "
        f"LVTKL={dl:6.3f} "
        f"dual={dual}"
    )


with OUT_FILE.open(
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=rows[0].keys()
    )

    writer.writeheader()
    writer.writerows(rows)


dual_n = sum(
    r["dual_contact_le4A"]
    for r in rows
)

s_n = sum(
    r["slahs_contact_le4A"]
    for r in rows
)

l_n = sum(
    r["lvtkl_contact_le4A"]
    for r in rows
)


print()
print("=" * 70)
print("HISTORICAL SEED 20260911")
print("=" * 70)

print("Total =", len(rows))

print(
    f"SLAHS contact = "
    f"{s_n}/48 "
    f"({100*s_n/48:.1f}%)"
)

print(
    f"LVTKL contact = "
    f"{l_n}/48 "
    f"({100*l_n/48:.1f}%)"
)

print(
    f"Dual contact = "
    f"{dual_n}/48 "
    f"({100*dual_n/48:.1f}%)"
)

print()
print(
    "Historical best score =",
    min(
        r["vina_kcal_mol"]
        for r in rows
    )
)

print("Output =", OUT_FILE)
