from pathlib import Path
import math
from collections import Counter

BASE = Path(__file__).resolve().parent

RECEPTOR_DIR = BASE / "refined_receptors"
DOCK_DIR = BASE / "refinement_docking_out"

MODELS = [
    "COMBO_011_R28",
    "COMBO_011_R01",
    "COMBO_011_R22",
    "COMBO_044_R28",
    "COMBO_030_R11",
]

CUTOFF = 4.5

# Protein side-chain carbon atoms only.
# Backbone carbon atoms are excluded to avoid inflating hydrophobic counts.
BACKBONE_NAMES = {"C", "CA"}

# Ligand atoms treated as hydrophobic/aromatic contact partners.
LIGAND_HYDROPHOBIC_TYPES = {"A", "C", "S", "Br"}


def parse_atom(line):
    try:
        return {
            "serial": int(line[6:11]),
            "name": line[12:16].strip(),
            "res": line[17:20].strip(),
            "chain": line[21:22].strip(),
            "resnum": line[22:26].strip(),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "type": line.split()[-1],
        }
    except Exception:
        return None


def read_receptor(path):
    atoms = []

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom = parse_atom(line)
                if atom:
                    atoms.append(atom)

    return atoms


def read_first_vina_model(path):
    atoms = []
    inside = False

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:

            if line.startswith("MODEL"):
                if not inside:
                    inside = True
                    continue

            if line.startswith("ENDMDL") and inside:
                break

            if inside and line.startswith(("ATOM", "HETATM")):
                atom = parse_atom(line)
                if atom:
                    atoms.append(atom)

    return atoms


def distance(a, b):
    return math.sqrt(
        (a["x"] - b["x"]) ** 2 +
        (a["y"] - b["y"]) ** 2 +
        (a["z"] - b["z"]) ** 2
    )


def peptide_name(chain):
    if chain == "A":
        return "SLAHS"
    if chain == "B":
        return "LVTKL"
    return chain


def analyze_model(model):

    receptor = read_receptor(
        RECEPTOR_DIR / f"{model}.pdbqt"
    )

    ligand = read_first_vina_model(
        DOCK_DIR / f"{model}_out.pdbqt"
    )

    protein_sidechain_c = [
        a for a in receptor
        if a["type"] == "C"
        and a["name"] not in BACKBONE_NAMES
    ]

    ligand_hydrophobic = [
        a for a in ligand
        if a["type"] in LIGAND_HYDROPHOBIC_TYPES
    ]

    atom_contacts = []
    residue_min = {}

    for p in protein_sidechain_c:

        for lig in ligand_hydrophobic:

            d = distance(p, lig)

            if d > CUTOFF:
                continue

            pep = peptide_name(p["chain"])

            atom_contacts.append({
                "peptide": pep,
                "res": p["res"],
                "resnum": p["resnum"],
                "protein_atom": p["name"],
                "lig_serial": lig["serial"],
                "lig_name": lig["name"],
                "lig_type": lig["type"],
                "distance": d,
            })

            key = (
                pep,
                p["res"],
                p["resnum"],
            )

            if key not in residue_min:
                residue_min[key] = d
            else:
                residue_min[key] = min(
                    residue_min[key],
                    d
                )

    atom_contacts.sort(
        key=lambda x: x["distance"]
    )

    return atom_contacts, residue_min


all_results = {}
residue_presence = Counter()

print()
print("STAGE 3B-4 HYDROPHOBIC CONTACT COMPARISON")
print("=" * 90)
print(f"Side-chain heavy-atom contact cutoff = {CUTOFF:.1f} A")
print()


for model in MODELS:

    contacts, residue_min = analyze_model(model)

    all_results[model] = (
        contacts,
        residue_min,
    )

    print(model)
    print("-" * 90)

    slahs_contacts = [
        x for x in contacts
        if x["peptide"] == "SLAHS"
    ]

    lvtkl_contacts = [
        x for x in contacts
        if x["peptide"] == "LVTKL"
    ]

    print(
        f"  atom-pair contacts: "
        f"SLAHS={len(slahs_contacts)} "
        f"LVTKL={len(lvtkl_contacts)} "
        f"total={len(contacts)}"
    )

    print("  contacting residues:")

    if not residue_min:
        print("    none")

    for key, d in sorted(
        residue_min.items(),
        key=lambda x: (
            x[0][0],
            int(x[0][2]),
        )
    ):
        pep, res, resnum = key

        print(
            f"    {pep:5s} "
            f"{res}{resnum} "
            f"min={d:.2f} A"
        )

        residue_presence[key] += 1

    print()

    print("  closest atom-pair contacts:")

    for x in contacts[:12]:
        print(
            f"    {x['peptide']:5s} "
            f"{x['res']}{x['resnum']}:{x['protein_atom']} "
            f"-> Lig#{x['lig_serial']}:{x['lig_name']} "
            f"({x['lig_type']}) "
            f"d={x['distance']:.2f} A"
        )

    print()


print()
print("CONSERVED CONTACTING RESIDUES")
print("=" * 90)

for key, count in residue_presence.most_common():

    pep, res, resnum = key

    print(
        f"{count}/5 poses  "
        f"{pep:5s} "
        f"{res}{resnum}"
    )


print()
print("MODEL-LEVEL HYDROPHOBIC SUPPORT TO BOTH HOTSPOTS")
print("=" * 90)

dual_count = 0

for model, (contacts, residue_min) in all_results.items():

    has_slahs = any(
        key[0] == "SLAHS"
        for key in residue_min
    )

    has_lvtkl = any(
        key[0] == "LVTKL"
        for key in residue_min
    )

    dual = has_slahs and has_lvtkl

    if dual:
        dual_count += 1

    print(
        f"{model}: "
        f"SLAHS_hydrophobic={has_slahs} "
        f"LVTKL_hydrophobic={has_lvtkl} "
        f"DUAL_HYDROPHOBIC_SUPPORT={dual}"
    )


print()
print(
    f"Representative poses with hydrophobic "
    f"support to both hotspots: "
    f"{dual_count}/{len(MODELS)}"
)