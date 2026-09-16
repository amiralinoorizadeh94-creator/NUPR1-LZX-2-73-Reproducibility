from pathlib import Path

BASE = Path(__file__).resolve().parent

RECEPTOR = (
    BASE / "refined_receptors" /
    "COMBO_011_R28.pdbqt"
)

DOCKED = (
    BASE / "refinement_docking_out" /
    "COMBO_011_R28_out.pdbqt"
)

OUTDIR = BASE / "md_preparation"
OUTDIR.mkdir(exist_ok=True)

OUTPUT = OUTDIR / "COMBO_011_R28_complex_raw.pdb"


def pdb_atom_line(
    serial,
    name,
    resname,
    chain,
    resnum,
    x,
    y,
    z,
    element
):
    return (
        f"HETATM{serial:5d} "
        f"{name:<4s} "
        f"{resname:>3s} "
        f"{chain:1s}"
        f"{resnum:4d}    "
        f"{x:8.3f}"
        f"{y:8.3f}"
        f"{z:8.3f}"
        f"{1.00:6.2f}"
        f"{0.00:6.2f}"
        f"          "
        f"{element:>2s}\n"
    )


def parse_pdbqt_atom(line):
    return {
        "name": line[12:16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21:22].strip(),
        "resnum": int(line[22:26]),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "type": line.split()[-1],
    }


def element_from_type(atom_type):
    mapping = {
        "C": "C",
        "A": "C",
        "N": "N",
        "NA": "N",
        "OA": "O",
        "O": "O",
        "S": "S",
        "SA": "S",
        "HD": "H",
        "H": "H",
        "Br": "Br",
    }

    return mapping.get(atom_type, atom_type[:2])


# ==========================================================
# READ RECEPTOR
# ==========================================================

receptor_atoms = []

with open(
    RECEPTOR,
    encoding="utf-8",
    errors="ignore"
) as f:

    for line in f:

        if not line.startswith(
            ("ATOM", "HETATM")
        ):
            continue

        atom = parse_pdbqt_atom(line)

        # Remove docking/preparation hydrogens.
        # Hydrogens will later be rebuilt during MD preparation.
        if atom["type"] in {"HD", "H"}:
            continue

        receptor_atoms.append(atom)


# ==========================================================
# READ ONLY VINA MODEL 1
# ==========================================================

ligand_atoms = []

inside_model1 = False

with open(
    DOCKED,
    encoding="utf-8",
    errors="ignore"
) as f:

    for line in f:

        if line.startswith("MODEL"):

            model_number = int(
                line.split()[1]
            )

            if model_number == 1:
                inside_model1 = True
                continue

            elif inside_model1:
                break

        if (
            inside_model1
            and line.startswith(
                ("ATOM", "HETATM")
            )
        ):

            atom = parse_pdbqt_atom(line)

            # Exclude explicit docking hydrogens here.
            if atom["type"] in {"HD", "H"}:
                continue

            ligand_atoms.append(atom)


# ==========================================================
# WRITE RAW PDB COMPLEX
# ==========================================================

serial = 1

with open(OUTPUT, "w") as out:

    out.write(
        "REMARK Raw coordinate-preserving complex for MD preparation\n"
    )

    out.write(
        "REMARK Receptor geometry: COMBO_011_R28\n"
    )

    out.write(
        "REMARK Ligand: LZX-2-73, AutoDock Vina mode 1\n"
    )

    out.write(
        "REMARK Hydrogens removed intentionally; "
        "they will be rebuilt during force-field preparation\n"
    )

    # ------------------------------------------------------
    # Receptor
    # ------------------------------------------------------

    for atom in receptor_atoms:

        element = element_from_type(
            atom["type"]
        )

        line = pdb_atom_line(
            serial=serial,
            name=atom["name"],
            resname=atom["resname"],
            chain=atom["chain"],
            resnum=atom["resnum"],
            x=atom["x"],
            y=atom["y"],
            z=atom["z"],
            element=element,
        )

        # Protein/peptide records should be ATOM
        line = "ATOM  " + line[6:]

        out.write(line)

        serial += 1

    out.write("TER\n")

    # ------------------------------------------------------
    # Ligand
    # ------------------------------------------------------

    for atom in ligand_atoms:

        element = element_from_type(
            atom["type"]
        )

        out.write(
            pdb_atom_line(
                serial=serial,
                name=atom["name"],
                resname="LZX",
                chain="L",
                resnum=1,
                x=atom["x"],
                y=atom["y"],
                z=atom["z"],
                element=element,
            )
        )

        serial += 1

    out.write("TER\n")
    out.write("END\n")


# ==========================================================
# SANITY CHECKS
# ==========================================================

EXPECTED_RECEPTOR_HEAVY = 76
EXPECTED_LIGAND_HEAVY = 24

print()
print(
    "STAGE 4A-1 RAW MD COMPLEX BUILD"
)
print("=" * 70)

print(
    f"Receptor heavy atoms = "
    f"{len(receptor_atoms)}"
)

print(
    f"Ligand heavy atoms   = "
    f"{len(ligand_atoms)}"
)

print(
    f"Total heavy atoms    = "
    f"{len(receptor_atoms) + len(ligand_atoms)}"
)

print()
print(
    f"Output: {OUTPUT}"
)

print()
print(
    f"Expected receptor heavy atoms = "
    f"{EXPECTED_RECEPTOR_HEAVY}"
)

if len(receptor_atoms) == EXPECTED_RECEPTOR_HEAVY:
    print(
        "Receptor atom-count check = PASS"
    )
else:
    print(
        "Receptor atom-count check = FAIL"
    )

print()
print(
    f"Expected ligand heavy atoms = "
    f"{EXPECTED_LIGAND_HEAVY}"
)

if len(ligand_atoms) == EXPECTED_LIGAND_HEAVY:
    print(
        "Ligand atom-count check = PASS"
    )
else:
    print(
        "Ligand atom-count check = FAIL"
    )

chains = sorted(
    set(
        atom["chain"]
        for atom in receptor_atoms
    )
)

print()
print(
    f"Receptor chains = {chains}"
)

if chains == ["A", "B"]:
    print(
        "Receptor chain check = PASS"
    )
else:
    print(
        "Receptor chain check = WARNING"
    )

print()
print(
    "IMPORTANT: This file preserves the docking coordinates "
    "but is NOT yet an MD-ready topology."
)