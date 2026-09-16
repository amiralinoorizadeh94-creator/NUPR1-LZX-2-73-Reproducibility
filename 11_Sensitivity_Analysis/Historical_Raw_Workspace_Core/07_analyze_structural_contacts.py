from pathlib import Path
import math

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

CONTACT_CUTOFF = 4.0  # Angstrom


def read_receptor_atoms(path):
    atoms = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom = parse_atom_line(line)
            if atom is not None:
                atoms.append(atom)

    return atoms


def read_first_vina_model(path):
    """
    Read only the first Vina docking pose (MODEL 1).
    If MODEL records are absent, read all atom records as one pose.
    """
    atoms = []
    inside_first_model = False
    model_records_seen = False

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:

            if line.startswith("MODEL"):
                model_records_seen = True

                if not inside_first_model:
                    inside_first_model = True
                    continue

            if line.startswith("ENDMDL"):
                if inside_first_model:
                    break

            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            if model_records_seen and not inside_first_model:
                continue

            atom = parse_atom_line(line)
            if atom is not None:
                atoms.append(atom)

    return atoms


def parse_atom_line(line):
    try:
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21:22].strip()
        res_num = line[22:26].strip()

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

    except Exception:
        return None

    # Ignore hydrogens where identifiable.
    element_guess = atom_name[0].upper() if atom_name else ""

    if element_guess == "H":
        return None

    return {
        "atom": atom_name,
        "res": res_name,
        "chain": chain,
        "resnum": res_num,
        "xyz": (x, y, z),
    }


def distance(a, b):
    return math.sqrt(
        (a[0] - b[0]) ** 2
        + (a[1] - b[1]) ** 2
        + (a[2] - b[2]) ** 2
    )


def split_receptor_by_chain(atoms):
    chains = {}

    for atom in atoms:
        chains.setdefault(atom["chain"], []).append(atom)

    if "A" not in chains or "B" not in chains:
        raise RuntimeError(
            f"Expected receptor chains A and B, found: {list(chains.keys())}"
        )

    # Confirmed from peptide sequence:
    # Chain A = SLAHS
    # Chain B = LVTKL
    return chains["A"], chains["B"]


def min_distance(group_a, group_b):
    return min(
        distance(a["xyz"], b["xyz"])
        for a in group_a
        for b in group_b
    )


def count_contacts(group_a, group_b, cutoff):
    count = 0

    for a in group_a:
        for b in group_b:
            if distance(a["xyz"], b["xyz"]) <= cutoff:
                count += 1

    return count


def residue_contacts(peptide_atoms, ligand_atoms, cutoff):
    contacts = {}

    for p in peptide_atoms:
        key = (p["res"], p["resnum"])

        for lig in ligand_atoms:
            d = distance(p["xyz"], lig["xyz"])

            if d <= cutoff:
                if key not in contacts or d < contacts[key]:
                    contacts[key] = d

    return contacts


print()
print("STAGE 3A STRUCTURAL CONTACT ANALYSIS - VINA MODE 1 ONLY")
print("=" * 76)
print(f"Heavy-atom contact cutoff = {CONTACT_CUTOFF:.1f} A")
print("Chain A = SLAHS")
print("Chain B = LVTKL")
print()


for model in MODELS:

    receptor_file = RECEPTOR_DIR / f"{model}.pdbqt"
    ligand_file = DOCK_DIR / f"{model}_out.pdbqt"

    receptor_atoms = read_receptor_atoms(receptor_file)
    ligand_atoms = read_first_vina_model(ligand_file)

    slahs, lvtkl = split_receptor_by_chain(receptor_atoms)

    min_slahs = min_distance(slahs, ligand_atoms)
    min_lvtkl = min_distance(lvtkl, ligand_atoms)

    contacts_slahs = count_contacts(
        slahs, ligand_atoms, CONTACT_CUTOFF
    )

    contacts_lvtkl = count_contacts(
        lvtkl, ligand_atoms, CONTACT_CUTOFF
    )

    residues_slahs = residue_contacts(
        slahs, ligand_atoms, CONTACT_CUTOFF
    )

    residues_lvtkl = residue_contacts(
        lvtkl, ligand_atoms, CONTACT_CUTOFF
    )

    close_slahs = contacts_slahs > 0
    close_lvtkl = contacts_lvtkl > 0

    if close_slahs and close_lvtkl:
        classification = "DUAL_CONTACT"
    elif close_slahs:
        classification = "SLAHS_ONLY"
    elif close_lvtkl:
        classification = "LVTKL_ONLY"
    else:
        classification = "NO_CLOSE_CONTACT"

    print(model)
    print(f"  receptor heavy atoms : {len(receptor_atoms)}")
    print(f"  ligand heavy atoms   : {len(ligand_atoms)}")
    print(f"  SLAHS heavy atoms    : {len(slahs)}")
    print(f"  LVTKL heavy atoms    : {len(lvtkl)}")

    print(f"  min ligand-SLAHS     : {min_slahs:.3f} A")
    print(f"  min ligand-LVTKL     : {min_lvtkl:.3f} A")

    print(f"  contacts SLAHS       : {contacts_slahs}")
    print(f"  contacts LVTKL       : {contacts_lvtkl}")

    print(f"  classification       : {classification}")

    print("  SLAHS residues       :", end=" ")
    if residues_slahs:
        print(
            ", ".join(
                f"{res}{num}({d:.2f}A)"
                for (res, num), d in residues_slahs.items()
            )
        )
    else:
        print("none")

    print("  LVTKL residues       :", end=" ")
    if residues_lvtkl:
        print(
            ", ".join(
                f"{res}{num}({d:.2f}A)"
                for (res, num), d in residues_lvtkl.items()
            )
        )
    else:
        print("none")

    print("-" * 76)