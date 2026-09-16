# v2.6 revision-analysis additions

This directory contains source data and scripts added after the v2.5 reproducibility package.

## Refined-pose clustering
`refined_pose_clustering.py` is the exact analysis script used in this v2.6 build to reproduce the reported 150-pose clustering from archived refined receptors and top-ranked docking outputs. It aligns the 30 peptide-backbone N/CA/C atoms to COMBO_011_R28, applies the same rigid transform to the ligand, computes pairwise ligand heavy-atom RMSD without ligand-only refitting, evaluates average-linkage solutions for k=2–10 by silhouette score, and writes cluster membership and the RMSD matrix.

Expected result: k=8; cluster populations 39, 28, 26, 25, 23, 6, 2, 1. COMBO_011_R28 is in cluster 1; cluster-1 medoid is COMBO_025_R09.

## Replicate 1 MIC analysis / Figure S2
`rep1_mic_analysis.py` is the exact analysis script used in this v2.6 build to derive the MIC interfragment source series from the supplied Replicate 1 topology and two DCD segments. In the compact package, the topology and DCD checksums/provenance are deposited; the large DCD payloads are not embedded. To rerun the script, provide a directory containing the topology and both DCD files as the first command-line argument. The minimum interfragment metric is all-atom atom–atom distance, including hydrogens. This is distinct from the ligand–fragment contact definition elsewhere in the study, which is heavy-atom distance <=4.0 Å.

`plot_figure_s2.py` regenerates a publication-style plot from the deposited processed CSV files. `Figure_S2_Rep1_MIC_interfragment_v2.6.png` is the canonical final audited Figure S2 asset embedded in Supplement v1.11.

## Final figure assets
The final audited assets embedded in Main v1.17 are deposited as Figure1, Figure2, Figure3, and Figure5 PNGs. Their underlying structural/source files remain in the inherited v2.5 directories. The original prior-session plotting code that produced the exact final pixel styling of Figures 1/2/3/5 was not available in the current build environment, so this package does not falsely label a reconstructed plotting script as the historical exact script. The canonical deposited PNGs and source coordinates/data are retained instead.

No MD simulation was rerun. No raw Replicate 2 or 3 trajectory is included or implied.
