from pathlib import Path
import math

BASE = Path(__file__).resolve().parent

RECEPTOR_FILE = BASE / "refined_receptors" / "COMBO_011_R28.pdbqt"
LIGAND_FILE = BASE / "refinement_docking_out" / "COMBO_011_R28_out.pdbqt"

# Geometric criteria
H_ACCEPTOR_MAX = 2.5      # H...A distance in Angstrom
D_ACCEPTOR_MAX = 3.5      # D...A distance in Angstrom
D_H_A_MIN_ANGLE = 135.0   # degrees
DONOR_H_BOND_MAX = 1.25   # donor-heavy-atom to H assignment


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


def xyz(a):
    return (a["x"], a["y"], a["z"])


def distance(a, b):
    ax, ay, az = xyz(a)
    bx, by, bz = xyz(b)

    return math.sqrt(
        (ax - bx) ** 2 +
        (ay - by) ** 2 +
        (az - bz) ** 2
    )


def angle(a, b, c):
    """
    Angle ABC in degrees.
    For hydrogen bonds we use D-H-A:
    a = donor
    b = hydrogen
    c = acceptor
    """

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

    norm_ba = math.sqrt(sum(x * x for x in ba))
    norm_bc = math.sqrt(sum(x * x for x in bc))

    if norm_ba == 0 or norm_bc == 0:
        return 0.0

    cosine = dot / (norm_ba * norm_bc)

    cosine = max(-1.0, min(1.0, cosine))

    return math.degrees(math.acos(cosine))


receptor = read_receptor(RECEPTOR_FILE)
ligand = read_first_vina_model(LIGAND_FILE)

# Heavy donor candidates.
# We restrict these to common N/O donor-capable atoms.
DONOR_HEAVY_TYPES = {"N", "NA", "OA"}

# Acceptor types available in these PDBQT files.
ACCEPTOR_TYPES = {"OA", "NA"}

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


def assign_hydrogens(hydrogens, donors):
    """
    Assign each explicit donor H to its nearest plausible donor heavy atom.
    """
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


receptor_assignments = assign_hydrogens(
    receptor_h,
    receptor_donors
)

ligand_assignments = assign_hydrogens(
    ligand_h,
    ligand_donors
)


def find_hbonds(assignments, acceptors, donor_source):
    results = []

    for h, donor, dh_dist in assignments:

        for acceptor in acceptors:

            # Prevent intramolecular self-contact if ever applicable
            if donor_source == "receptor":
                if acceptor["chain"] == donor["chain"]:
                    continue

            ha = distance(h, acceptor)
            da = distance(donor, acceptor)

            if ha > H_ACCEPTOR_MAX:
                continue

            if da > D_ACCEPTOR_MAX:
                continue

            dha = angle(donor, h, acceptor)

            if dha < D_H_A_MIN_ANGLE:
                continue

            results.append(
                {
                    "donor": donor,
                    "H": h,
                    "acceptor": acceptor,
                    "DH": dh_dist,
                    "HA": ha,
                    "DA": da,
                    "angle": dha,
                }
            )

    return sorted(results, key=lambda x: x["HA"])


# Ligand donor -> receptor acceptor
lig_to_rec = find_hbonds(
    ligand_assignments,
    receptor_acceptors,
    "ligand",
)

# Receptor donor -> ligand acceptor
rec_to_lig = []

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

        rec_to_lig.append(
            {
                "donor": donor,
                "H": h,
                "acceptor": acceptor,
                "DH": dh_dist,
                "HA": ha,
                "DA": da,
                "angle": dha,
            }
        )

rec_to_lig = sorted(rec_to_lig, key=lambda x: x["HA"])


print()
print("STAGE 3B-2 ANGLE-AWARE HYDROGEN-BOND VALIDATION")
print("=" * 82)

print(f"H...A cutoff       = {H_ACCEPTOR_MAX:.2f} A")
print(f"D...A cutoff       = {D_ACCEPTOR_MAX:.2f} A")
print(f"D-H-A minimum angle= {D_H_A_MIN_ANGLE:.1f} deg")
print()

print("HYDROGEN ASSIGNMENT SANITY CHECK")
print("-" * 82)

print(
    f"Receptor HD atoms: {len(receptor_h)}; "
    f"assigned to donor heavy atom: {len(receptor_assignments)}"
)

print(
    f"Ligand HD atoms: {len(ligand_h)}; "
    f"assigned to donor heavy atom: {len(ligand_assignments)}"
)

print()
print("Ligand donor-H assignments:")

for h, donor, d in ligand_assignments:
    print(
        f"  Lig H#{h['serial']} -> "
        f"Lig {donor['name']}#{donor['serial']} "
        f"({donor['type']}), D-H={d:.2f} A"
    )


print()
print("VALIDATED: LIGAND DONOR -> PEPTIDE ACCEPTOR")
print("-" * 82)

if lig_to_rec:

    for hb in lig_to_rec:

        d = hb["donor"]
        h = hb["H"]
        a = hb["acceptor"]

        peptide = (
            "SLAHS" if a["chain"] == "A"
            else "LVTKL" if a["chain"] == "B"
            else a["chain"]
        )

        print(
            f"{peptide:5s}  "
            f"Lig {d['name']}#{d['serial']}-H#{h['serial']} "
            f"-> {a['res']}{a['resnum']}:{a['name']}  "
            f"H...A={hb['HA']:.2f} A  "
            f"D...A={hb['DA']:.2f} A  "
            f"angle={hb['angle']:.1f} deg"
        )

else:
    print("None")


print()
print("VALIDATED: PEPTIDE DONOR -> LIGAND ACCEPTOR")
print("-" * 82)

if rec_to_lig:

    for hb in rec_to_lig:

        d = hb["donor"]
        h = hb["H"]
        a = hb["acceptor"]

        peptide = (
            "SLAHS" if d["chain"] == "A"
            else "LVTKL" if d["chain"] == "B"
            else d["chain"]
        )

        print(
            f"{peptide:5s}  "
            f"{d['res']}{d['resnum']}:{d['name']}-H#{h['serial']} "
            f"-> Lig {a['name']}#{a['serial']}  "
            f"H...A={hb['HA']:.2f} A  "
            f"D...A={hb['DA']:.2f} A  "
            f"angle={hb['angle']:.1f} deg"
        )

else:
    print("None")


print()
print("SUMMARY")
print("-" * 82)

sla_ligdon = sum(
    1 for hb in lig_to_rec
    if hb["acceptor"]["chain"] == "A"
)

lvt_ligdon = sum(
    1 for hb in lig_to_rec
    if hb["acceptor"]["chain"] == "B"
)

sla_protdon = sum(
    1 for hb in rec_to_lig
    if hb["donor"]["chain"] == "A"
)

lvt_protdon = sum(
    1 for hb in rec_to_lig
    if hb["donor"]["chain"] == "B"
)

print(
    f"Ligand-donor validated H-bonds: "
    f"SLAHS={sla_ligdon}, LVTKL={lvt_ligdon}"
)

print(
    f"Peptide-donor validated H-bonds: "
    f"SLAHS={sla_protdon}, LVTKL={lvt_protdon}"
)

print(
    f"Total validated geometric H-bonds = "
    f"{len(lig_to_rec) + len(rec_to_lig)}"
)