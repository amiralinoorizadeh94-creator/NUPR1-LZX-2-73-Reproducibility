
STAGE 2C — local geometric refinement around the five best two-hotspot seeds

What it does
------------
Uses the top 5 Stage-2B models:
COMBO_011, COMBO_031, COMBO_030, COMBO_044, COMBO_025

For each seed it generates 30 local perturbations by:
- rotating LVTKL by up to ±18 degrees
- translating it locally (typically <2 A)
- rejecting hard peptide-peptide clashes

Total target = 150 refined two-hotspot receptors.

Then it docks LZX-2-73 with Vina 1.2.3 using exhaustiveness 16 and reports
the score distribution.

Scientific meaning
------------------
This asks whether the stronger Stage-2B geometries lie near even better local
arrangements. It is still a rigid-receptor ensemble approximation, not an exact
reproduction of the paper's fully flexible peptide calculation.

Required in the working Stage-2B folder:
combined_receptors/
vina.exe
LZX-2-73.pdbqt

Run:
python 04_build_local_refinement.py
python 05_run_local_refinement_docking.py
python 06_analyze_local_refinement.py
