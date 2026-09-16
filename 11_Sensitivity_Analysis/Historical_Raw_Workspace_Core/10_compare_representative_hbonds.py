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

H_ACCEPTOR_MAX = 2.5
D_ACCEPTOR_MAX = 3.5
D_H_A_MIN_ANGLE = 135.0
DONOR_H_BOND_MAX = 1.25

DONOR_HEAVY_TYPES = {"N", "NA", "OA"}
ACCEPTOR_TYPES = {"OA", "NA"}


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


def angle(a, b, c):
    ba = (
        a["x"] - b["x"],
        a["y"] - b["y"],
        a["z"] - b["z"],
    )

    bc = (
        c["x"] - b["x"],
        c["y"] - b["y"],
        c["z"] - b["z"],
    )

    dot = sum(x * y for x, y in zip(ba, bc))

    nba = math.sqrt(sum(x * x for x in ba))
    nbc = math.sqrt(sum(x * x for x in bc))

    if nba == 0 or nbc == 0:
        return 0.0

    cosine = dot / (nba * nbc)
    cosine = max(-1.0, min(1.0, cosine))

    return math.degrees(math.acos(cosine))


def assign_hydrogens(hydrogens, donors):
    assignments = []

    for h in hydrogens:
        candidates = []

        for donor in donors:
            d = distance(h, donor)

            if d <= DONOR_H_BOND_MAX:
                candidates.append((d, donor))

        if candidates:
            candidates.sort(key=lambda x: x[0])
            d, donor = candidates[0]
            assignments.append((h, donor, d))

    return assignments


def analyze_model(model):

    receptor_file = RECEPTOR_DIR / f"{model}.pdbqt"
    ligand_file = DOCK_DIR / f"{model}_out.pdbqt"

    receptor = read_receptor(receptor_file)
    ligand = read_first_vina_model(ligand_file)

    receptor_h = [a for a in receptor if a["type"] == "HD"]
    ligand_h = [a for a in ligand if a["type"] == "HD"]

    receptor_heavy = [a for a in receptor if a["type"] != "HD"]
    ligand_heavy = [a for a in ligand if a["type"] != "HD"]

    receptor_donors = [
        a for a in receptor_heavy
        if a["type"] in DONOR_HEAVY_TYPES
    ]

    ligand_donors = [
        a for a in ligand_heavy
        if a["type"] in DONOR_HEAVY_TYPES
    ]

    receptor_acceptors = [
        a for a in receptor_heavy
        if a["type"] in ACCEPTOR_TYPES
    ]

    ligand_acceptors = [
        a for a in ligand_heavy
        if a["type"] in ACCEPTOR_TYPES
    ]

    receptor_assignments = assign_hydrogens(
        receptor_h,
        receptor_donors
    )

    ligand_assignments = assign_hydrogens(
        ligand_h,
        ligand_donors
    )

    hbonds = []

    # Ligand donor -> receptor acceptor
    for h, donor, dh_dist in ligand_assignments:

        for acceptor in receptor_acceptors:

            ha = distance(h, acceptor)
            da = distance(donor, acceptor)

            if ha > H_ACCEPTOR_MAX:
                continue

            if da > D_ACCEPTOR_MAX:
                continue

            dha = angle(donor, h, acceptor)

            if dha < D_H_A_MIN_ANGLE:
                continue

            peptide = (
                "SLAHS" if acceptor["chain"] == "A"
                else "LVTKL" if acceptor["chain"] == "B"
                else acceptor["chain"]
            )

            hbonds.append({
                "direction": "LIGAND->PEPTIDE",
                "peptide": peptide,
                "res": acceptor["res"],
                "resnum": acceptor["resnum"],
                "atom": acceptor["name"],
                "lig_atom": donor["serial"],
                "HA": ha,
                "DA": da,
                "angle": dha,
            })

    # Receptor donor -> ligand acceptor
    for h, donor, dh_dist in receptor_assignments:

        for acceptor in ligand_acceptors:

            ha = distance(h, acceptor)
            da = distance(donor, acceptor)

            if ha > H_ACCEPTOR_MAX:
                continue

            if da > D_ACCEPTOR_MAX:
                continue

            dha = angle(donor, h, acceptor)

            if dha < D_H_A_MIN_ANGLE:
                continue

            peptide = (
                "SLAHS" if donor["chain"] == "A"
                else "LVTKL" if donor["chain"] == "B"
                else donor["chain"]
            )

            hbonds.append({
                "direction": "PEPTIDE->LIGAND",
                "peptide": peptide,
                "res": donor["res"],
                "resnum": donor["resnum"],
                "atom": donor["name"],
                "lig_atom": acceptor["serial"],
                "HA": ha,
                "DA": da,
                "angle": dha,
            })

    hbonds.sort(key=lambda x: (x["peptide"], x["HA"]))

    return hbonds


all_results = {}
interaction_counter = Counter()

print()
print("STAGE 3B-3 REPRESENTATIVE H-BOND COMPARISON")
print("=" * 88)
print(
    f"Criteria: H...A <= {H_ACCEPTOR_MAX:.2f} A; "
    f"D...A <= {D_ACCEPTOR_MAX:.2f} A; "
    f"D-H-A >= {D_H_A_MIN_ANGLE:.1f} deg"
)
print()


for model in MODELS:

    hbonds = analyze_model(model)
    all_results[model] = hbonds

    print(model)
    print("-" * 88)

    if not hbonds:
        print("  No validated geometric H-bonds")

    for hb in hbonds:

        key = (
            hb["direction"],
            hb["peptide"],
            hb["res"],
            hb["resnum"],
            hb["atom"],
        )

        interaction_counter[key] += 1

        print(
            f"  {hb['direction']:16s} "
            f"{hb['peptide']:5s} "
            f"{hb['res']}{hb['resnum']}:{hb['atom']} "
            f"LigAtom#{hb['lig_atom']} "
            f"H...A={hb['HA']:.2f} A "
            f"D...A={hb['DA']:.2f} A "
            f"angle={hb['angle']:.1f} deg"
        )

    n_slahs = sum(1 for x in hbonds if x["peptide"] == "SLAHS")
    n_lvtkl = sum(1 for x in hbonds if x["peptide"] == "LVTKL")

    print(
        f"  SUMMARY: total={len(hbonds)} "
        f"SLAHS={n_slahs} "
        f"LVTKL={n_lvtkl}"
    )

    print()


print()
print("CONSERVED INTERACTIONS ACROSS REPRESENTATIVE POSES")
print("=" * 88)

for key, count in interaction_counter.most_common():

    direction, peptide, res, resnum, atom = key

    print(
        f"{count}/5 poses  "
        f"{direction:16s} "
        f"{peptide:5s} "
        f"{res}{resnum}:{atom}"
    )


print()
print("MODEL-LEVEL DUAL H-BOND SUPPORT")
print("=" * 88)

dual_count = 0

for model, hbonds in all_results.items():

    has_slahs = any(x["peptide"] == "SLAHS" for x in hbonds)
    has_lvtkl = any(x["peptide"] == "LVTKL" for x in hbonds)

    dual = has_slahs and has_lvtkl

    if dual:
        dual_count += 1

    print(
        f"{model}: "
        f"SLAHS_HB={has_slahs} "
        f"LVTKL_HB={has_lvtkl} "
        f"DUAL_HB_SUPPORT={dual}"
    )

print()
print(
    f"Representative poses with H-bond support to both hotspots: "
    f"{dual_count}/{len(MODELS)}"
)