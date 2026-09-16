
STAGE 2B — Combined SLAHS + LVTKL geometric ensemble docking

Goal
----
Test whether a combined two-hotspot environment can generate substantially
stronger LZX-2-73 Vina scores than the matched single-fragment controls.

Important scientific limitation
-------------------------------
The paper did not publish the exact starting coordinates of its two flexible
peptide fragments. Therefore this is NOT an exact reproduction of the paper's
-8.2 kcal/mol pose.

Instead, this stage performs a transparent geometric ensemble search:
- samples independent SLAHS/LVTKL conformers and relative orientations,
- rejects peptide-peptide steric clashes,
- retains 48 plausible local two-fragment arrangements,
- docks LZX-2-73 into every combined receptor with the same 30 A box.

This lets us ask a defensible question:
Does adding a second hotspot create a lower-scoring regime in an independent
ensemble, compared with the SLAHS (~-3.64 median) and LVTKL (~-4.00 median)
controls?

Run in active (nupr1-md):
  python 00_collect_inputs.py
  python 01_build_combined_ensemble.py
  python 02_run_two_hotspot_docking.py
  python 03_analyze_results.py

The first script automatically searches Downloads for the required previous
files/folders, so you do not need to copy them manually.

Main outputs:
  combined_receptors/ensemble_metadata.csv
  two_hotspot_docking_out/summary.csv
  two_hotspot_docking_out/*_out.pdbqt
