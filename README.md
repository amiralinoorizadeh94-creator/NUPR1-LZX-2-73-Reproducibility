# NUPR1-LZX-2-73 Reproducibility

Computational scripts, selected inputs, processed source data, representative structures, and analysis materials supporting the manuscript **“Modeling Dual-Hotspot Contact Geometry of LZX-2-73 with NUPR1 Peptide Fragments: Stochastic Fragment-Conformer Docking, Structural Refinement, and Molecular Dynamics.”**

**Author:** Amirali Noorizadehsalout, Department of Chemistry, Shahid Beheshti University, Tehran, Iran  
**ORCID:** 0009-0008-7015-4386  
**Archived release:** Zenodo DOI `10.5281/zenodo.22774318`

## Scope and interpretation

The study uses two untethered five-residue NUPR1 hotspot fragments, SLAHS (S31–S35) and LVTKL (L66–L70). The repository supports reproducibility of the reduced-model computational workflow. It does not establish a unique or native full-length NUPR1 binding mode, binding free energy, equilibrium occupancy, residence time, dissociation kinetics, or biological efficacy.

## Repository contents

- `02_Historical_Stage2B_Inputs/` — recovered historical receptor/ligand docking inputs and input manifest.
- `03_Historical_Stage2B_Stage2C_Scripts/` — historical reconstruction, docking, refinement, and analysis scripts.
- `05_Structural_Interactions/` — selected COMBO_011_R28 structural source data and verification material.
- `06_MD_System/` — reduced-system HMR topology/restart and ion-concentration check.
- `08_Trajectory_Analysis/` — processed three-replicate trajectory-analysis outputs.
- `10_Supplementary_Data/` — verified Supplementary Tables S1–S6 source data and representative structures.
- `11_Sensitivity_Analysis/` — frozen 288-model sensitivity summary products plus recovered sensitivity-analysis scripts.
- `12_v2.6_Revision_Analysis/` — Replicate 1 MIC analysis, refined-pose clustering, Figure S2 source data, and associated scripts.

## Core workflow

The historical workflow reconstructed 48 accepted dual-fragment configurations from fixed 12-member SLAHS and LVTKL conformer libraries, performed AutoDock Vina screening, locally refined the five highest-ranking initial configurations (150 refined models), and selected COMBO_011_R28 as the highest-ranking refined configuration. Five additional reconstruction seeds generated 240 models for within-protocol sensitivity analysis, giving 288 configurations in the combined frozen dataset.

The molecular-dynamics analysis comprises three 50 ns explicit-solvent replicate trajectories. Processed outputs for all three replicates are retained here. Raw-coordinate retrospective MIC analysis is available for Replicate 1 only. Equivalent raw-coordinate reanalysis could not be performed for Replicates 2 and 3 because their raw trajectories were unavailable.

## Important reproducibility boundary

This GitHub repository is a curated companion to the archived research release, not a replacement for Zenodo. Large raw production trajectories are intentionally excluded from Git history. The authoritative archived release and deposited source materials are identified by Zenodo DOI `10.5281/zenodo.22774318`.

Supplementary Table S4 is restricted to verified residue-level values. Unavailable interval-residue values were not reconstructed or imputed.

## Software recorded in the project

- AutoDock Vina 1.2.3
- OpenMM 8.6
- CPPTRAJ 7.6.2
- Amber ff19SB for peptide fragments
- GAFF2 with AM1-BCC charges for LZX-2-73
- OPC water model
- NVIDIA RTX 4090 for production MD

Ancillary software versions are not asserted where they were not independently preserved in project provenance.

## Running the code

The scripts are preserved primarily for provenance and reproducibility. Several historical scripts retain the paths, executable assumptions, or workspace conventions used during the original calculations. Inspect and adapt input/output paths before execution on a new machine. The historical Vina workflow did not specify a fixed internal Vina random seed, so rerunning it should not be expected to reproduce every docking score bit-for-bit.

A useful starting sequence for the historical reconstruction/docking workflow is:

```text
00_collect_inputs.py
01_build_combined_ensemble.py
02_run_two_hotspot_docking.py
03_analyze_results.py
04_geometry_analysis.py
04_build_local_refinement.py
05_run_local_refinement_docking.py
06_analyze_local_refinement.py
```

Revision analyses are documented separately in `12_v2.6_Revision_Analysis/README_v2.6_NEW_ANALYSES.md`.

## Licenses

Code is released under the MIT License. Research data/source-data files are provided under CC BY 4.0 unless otherwise noted.

## Citation

For an immutable record of the research materials, cite the Zenodo release:

**Noorizadehsalout, A. (2026). NUPR1 / LZX-2-73 computational reproducibility materials. Zenodo. DOI: 10.5281/zenodo.22774318.**

The manuscript should also be cited after publication.
