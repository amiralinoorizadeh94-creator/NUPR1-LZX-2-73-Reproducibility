from pathlib import Path
import math

BASE = Path(__file__).resolve().parent

RECEPTOR_FILE = BASE / "refined_receptors" / "COMBO_011_R28.pdbqt"
LIGAND_FILE = BASE / "refinement_docking_out" / "COMBO_011_R28_out.pdbqt"

HBOND_CUTOFF = 3.5
HYDROPHOBIC_CUTOFF = 4.5
BR_CUTOFF = 4.5


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


def dist(a, b):
    return math.sqrt(
        (a["x"] - b["x"]) ** 2
        + (a["y"] - b["y"]) ** 2
        + (a["z"] - b["z"]) ** 2
    )


receptor = read_receptor(RECEPTOR_FILE)
ligand = read_first_vina_model(LIGAND_FILE)

slahs = [a for a in receptor if a["chain"] == "A"]
lvtkl = [a for a in receptor if a["chain"] == "B"]

# AutoDock-style acceptor types
ACCEPTOR_TYPES = {"OA", "NA"}

# Protein donor nitrogens/oxygens carrying explicit donor H nearby
DONOR_PARENT_TYPES = {"N", "NA", "OA"}

# Ligand donor hydrogens
ligand_H = [a for a in ligand if a["type"] == "HD"]

# Ligand acceptors
ligand_acceptors = [a for a in ligand if a["type"] in ACCEPTOR_TYPES]

# Ligand hydrophobic / aromatic heavy atoms
ligand_hydrophobic = [
    a for a in ligand
    if a["type"] in {"A", "C", "S", "Br"}
]

br_atoms = [a for a in ligand if a["type"] == "Br"]


def protein_acceptor_candidates(peptide_atoms):
    return [a for a in peptide_atoms if a["type"] in ACCEPTOR_TYPES]


def protein_hydrophobic_atoms(peptide_atoms):
    return [
        a for a in peptide_atoms
        if a["type"] == "C"
    ]


def ligand_H_to_protein_acceptors(peptide_atoms, peptide_name):
    results = []

    acceptors = protein_acceptor_candidates(peptide_atoms)

    for h in ligand_H:
        for acc in acceptors:
            d = dist(h, acc)

            if d <= HBOND_CUTOFF:
                results.append(
                    (
                        d,
                        peptide_name,
                        h["serial"],
                        acc["res"],
                        acc["resnum"],
                        acc["name"],
                        acc["type"],
                    )
                )

    return sorted(results)


def protein_to_ligand_acceptor_candidates(peptide_atoms, peptide_name):
    """
    Distance-based candidate screen only.
    Since receptor donor-H topology is incomplete for rigorous angular
    validation, these are labeled candidates rather than confirmed H-bonds.
    """
    results = []

    protein_polar = [
        a for a in peptide_atoms
        if a["type"] in DONOR_PARENT_TYPES
    ]

    for p in protein_polar:
        for lig in ligand_acceptors:
            d = dist(p, lig)

            if d <= HBOND_CUTOFF:
                results.append(
                    (
                        d,
                        peptide_name,
                        p["res"],
                        p["resnum"],
                        p["name"],
                        p["type"],
                        lig["serial"],
                        lig["name"],
                        lig["type"],
                    )
                )

    return sorted(results)


def hydrophobic_contacts(peptide_atoms, peptide_name):
    results = []

    prot_c = protein_hydrophobic_atoms(peptide_atoms)

    for p in prot_c:
        for lig in ligand_hydrophobic:
            d = dist(p, lig)

            if d <= HYDROPHOBIC_CUTOFF:
                results.append(
                    (
                        d,
                        peptide_name,
                        p["res"],
                        p["resnum"],
                        p["name"],
                        lig["serial"],
                        lig["name"],
                        lig["type"],
                    )
                )

    return sorted(results)


def br_environment(peptide_atoms, peptide_name):
    results = []

    for br in br_atoms:
        for p in peptide_atoms:

            if p["type"] == "HD":
                continue

            d = dist(br, p)

            if d <= BR_CUTOFF:
                results.append(
                    (
                        d,
                        peptide_name,
                        p["res"],
                        p["resnum"],
                        p["name"],
                        p["type"],
                    )
                )

    return sorted(results)


print()
print("STAGE 3B SPECIFIC INTERACTION ANALYSIS")
print("=" * 78)
print("Model: COMBO_011_R28")
print("Vina mode: 1 only")
print()

print("Ligand atom inventory:")
for atom in ligand:
    print(
        f"  {atom['serial']:>2} {atom['name']:<4} "
        f"type={atom['type']}"
    )

print()
print("H-BOND CANDIDATES: ligand donor H -> peptide acceptor")
print("-" * 78)

hb1 = ligand_H_to_protein_acceptors(slahs, "SLAHS")
hb2 = ligand_H_to_protein_acceptors(lvtkl, "LVTKL")

for x in hb1 + hb2:
    print(
        f"{x[1]:5s}  d={x[0]:.2f} A  "
        f"LigH#{x[2]} -> {x[3]}{x[4]}:{x[5]} ({x[6]})"
    )

if not hb1 and not hb2:
    print("None within cutoff")

print()
print("POLAR CONTACT / H-BOND CANDIDATES: peptide polar atom -> ligand acceptor")
print("-" * 78)

ph1 = protein_to_ligand_acceptor_candidates(slahs, "SLAHS")
ph2 = protein_to_ligand_acceptor_candidates(lvtkl, "LVTKL")

for x in ph1 + ph2:
    print(
        f"{x[1]:5s}  d={x[0]:.2f} A  "
        f"{x[2]}{x[3]}:{x[4]} ({x[5]}) -> "
        f"Lig#{x[6]}:{x[7]} ({x[8]})"
    )

if not ph1 and not ph2:
    print("None within cutoff")

print()
print("HYDROPHOBIC CONTACT CANDIDATES")
print("-" * 78)

hyd1 = hydrophobic_contacts(slahs, "SLAHS")
hyd2 = hydrophobic_contacts(lvtkl, "LVTKL")

for x in (hyd1 + hyd2)[:40]:
    print(
        f"{x[1]:5s}  d={x[0]:.2f} A  "
        f"{x[2]}{x[3]}:{x[4]} -> "
        f"Lig#{x[5]}:{x[6]} ({x[7]})"
    )

print()
print(f"SLAHS hydrophobic atom-pair contacts = {len(hyd1)}")
print(f"LVTKL hydrophobic atom-pair contacts = {len(hyd2)}")

print()
print("BROMINE ENVIRONMENT")
print("-" * 78)

br1 = br_environment(slahs, "SLAHS")
br2 = br_environment(lvtkl, "LVTKL")

for x in br1 + br2:
    print(
        f"{x[1]:5s}  d={x[0]:.2f} A  "
        f"{x[2]}{x[3]}:{x[4]} ({x[5]})"
    )

if not br1 and not br2:
    print("No peptide heavy atom within cutoff")

print()
print("SUMMARY")
print("-" * 78)

print(
    f"Ligand-donor H-bond candidates: "
    f"SLAHS={len(hb1)}, LVTKL={len(hb2)}"
)

print(
    f"Peptide-polar/ligand-acceptor candidates: "
    f"SLAHS={len(ph1)}, LVTKL={len(ph2)}"
)

print(
    f"Hydrophobic atom-pair contacts: "
    f"SLAHS={len(hyd1)}, LVTKL={len(hyd2)}"
)

print(
    f"Br-neighbor atoms: "
    f"SLAHS={len(br1)}, LVTKL={len(br2)}"
)