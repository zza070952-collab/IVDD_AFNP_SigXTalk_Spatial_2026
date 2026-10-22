Figure 4 was regenerated locally from the real completed V3 depth tables and the
unchanged V2 all-bin score and regional candidate summaries. The PNG was opened
for visual inspection. All six panels, tick labels, scale bars, legend entries
and annotation text are legible without material overlap. The final SVG uses
the same layout and data as the inspected PNG; top-row point clouds are
rasterized inside vector containers, while text, axes and ROI outlines remain
editable. High-resolution original H&E overlays accompany the compact figure.

Scientific checks:

- The exact frozen ROI SHA256 is
  `ef8cf1b835a0dbb9e938f5a8be3ca55757bc8c1ebcda4bb218620d2cac72fd31`.
  Both original image hashes and image dimensions are checked before plotting.
- The original map source has 12,464 CTR and 12,508 cKO QC bins. All original
  QC values determine the common color scale; explicit view-window flags retain
  points outside the displayed crop in the full source TSV. Coordinates were
  not transformed and values were not smoothed, interpolated or imputed.
- Original disc contrasts are recomputed from the baseline regional table and
  matched to the completed depth table within absolute tolerance 1e-9.
- Panel e retains eight original discs, connects the seven adequately covered
  discs through all three conditions, and displays the insufficient-retention
  original disc as an open point. Black diamonds always use those same seven
  discs. The original-eligible mean is 0.179486702860283; the 100-UMI mean is
  0.0316979897236584. Six of seven adequate pairs remain NP greater than AF.
- cKO_BL has only three retained AF bins and therefore fails the predeclared
  minimum of five. Its sampled estimate is retained in the full source table
  rather than included in the adequate-pair summary. cKO_TR is adequately
  covered but reverses direction; that negative contrast remains visible.
- Candidate dots are explicitly unadjusted and give each disc equal weight.
  PRSS23 is a very small nonzero CTR-NP dot and zero-detection hollow dots in
  the other capture-region combinations. This follows one real detected UMI
  in one CTR NP-core bin; it is not a spatial-validation claim.
- No animal n=8, genotype significance, bin-level P value, CEP ROI, ligand-to-
  receptor causality, or independent human spatial replication is asserted.

Material scientific implication: the plot now weakens the original spatial
NP-enrichment interpretation because sequencing depth substantially affects the
score. It cannot establish that the residual contrast is robust to all depth
or tissue-composition effects, or that the biological program is absent.
