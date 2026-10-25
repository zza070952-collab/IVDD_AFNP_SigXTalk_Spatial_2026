# Human figure authoring from frozen caches

`build_human_figures.py` rebuilds the human Figure 1, Figure 2 and the original V4 supplementary Figures S1/S2. The latter two are the full internal authoring views, not the later three-file compact supplement. It also exports their source tables, a source manifest, descriptive fact registry and figure contracts.

The code does not refit cNMF, rerun QC or identity assignment, calculate cell scores, perform clinical-grade tests, train a model or generate a cell embedding. Processed expression remains processed expression. Donors are the descriptive sampling unit. Missing comparator scores and unavailable mappings remain distinct from zero cell counts or measured zero detection. AUCell and UCell use separate scales. The locked K=5/top25 definition, original illustrative marker panels and figure geometry are method settings. Donor identifiers, totals, comparison counts, mapping ranges and review classifications are read or computed from inputs.

## Run

From the repository root:

```sh
python scripts/authoring/v4/human/build_human_figures.py \
  --source-root /path/to/frozen-review-cache \
  --v4-root /path/to/read-only-v4-cache \
  --output-root analysis/human-authoring-v4
```

`IVDD_REVIEW_ROOT` and `IVDD_V4_ROOT` may replace the corresponding arguments. `--eligibility-table` and `--design-table` may replace the V4-root requirements individually. Output must be outside both input roots. Use a new output directory; the script writes files beneath it. `analysis/` is ignored by this repository. Do not force-add generated output.

The exact tested plotting dependencies are in `requirements.txt`. Python 3.8.5 was used for the local cache smoke run. This does not establish compatibility with every environment or a clean-machine reconstruction of all analysis steps.

## Required external cache contract

The source root is the **review** directory, not the repository or the V4 output directory. It must contain these existing frozen files; this code-only publication does not bundle or download them:

```text
01_human_programs/
  COHORT_AND_DONOR_MANIFEST.tsv
  FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv
  FIGURE2_DONOR_PROGRAM_SOURCE.tsv
  FROZEN_K5_PROGRAMS.tsv
  PROGRAM_GENE_COVERAGE.tsv
  PROGRAM_GENE_MAPPING.tsv
  AUCELL_EFFECTIVE_PROGRAM_COVERAGE.tsv
  TU_MARKER_DOTPLOT_SOURCE.tsv
  TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv
  TU_CANDIDATE_GENE_BY_DONOR_IDENTITY.tsv
  ALL_PROGRAMS_BRIEF.tsv
  GEP2_DONOR_SUPPORT.tsv
  TU_MATRIX_AND_QC_MANIFEST.tsv
  AUCELL_RUN_METADATA.json
  PROVENANCE_AND_SCORE_SCALE.md
05_figures/Figure1_Figure2/Figure1a_donor_tissue_counts.tsv
02_sigxtalk_official_new/formal_new/tables/DONOR_GROUP_INPUT_COUNTS.tsv
07_public_repo/resources/frozen_k5/K_selection_and_parameters.tsv
07_public_repo/resources/frozen_k5/PROGRAM_LOCK.json
```

The V4 cache supplies:

```text
00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv
03_tables/human/Figure1A_source.tsv
```

The last file is a small frozen **input**, with columns `layer`, `label`, `scope`, describing the four study-design boxes. Reading this source avoids embedding observed cross-module counts in public code and avoids reading experimental wells. The script does not reconstruct this input from a manuscript. The schematic retains any source-status wording supplied in that table; a successful run cannot authenticate experiments.

Tables are tab-delimited with their frozen column names. Cohort and donor manifests supply `cohort`, `donor_id`, `gsm`, `submitted_cells` and `QC_cells`. Identity counts use `cohort`, `donor_id`, `identity`, `n_cells`. Scores use `cohort`, `donor_id`, `identity`, `program`, `score_method`, `score_median`. Gene tables retain original and mapped symbols, identity, donor, detection fraction and mean log1p published expression. Program membership uses `program`, `rank`, `gene`. The exact cached schemas are required; these are not generic single-cell importers. Eligibility preserves donor ID, actual tail counts, quartile cutoffs and eligibility. No raw cell-expression matrix or wet-lab record is required by this entry point.

## Outputs and verification boundary

Generated figures go to `04_figures/` and `06_supplementary/figures/`; generated tables go to `03_tables/human/`; structural QA, figure contracts and the input-hash manifest go to `08_review/human/`. Exports are editable PDF/SVG, 300-dpi PNG and 600-dpi TIFF. No Methods, Results, full legends, manuscript documents or original data are copied into this source directory or emitted as prose.

`PASS_AUTOMATED` reports only bounding-box and structural checks. Every run retains `PENDING_ACTUAL_INSPECTION` for manual image review. Rendering is not a substitute for that review. The source manifest records the supplied cache bytes, not independent source authentication. Syntax/help checks, cache smoke runs and exact numerical-source comparisons can validate this authoring step; they do not certify a complete raw-to-final analysis, an experimental source, or a publication release.
