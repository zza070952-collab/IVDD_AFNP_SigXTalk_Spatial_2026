# Communication and spatial figure authoring

`plot_communication_spatial_v4.py` renders the V4 communication/spatial main panels and the earlier supplementary panels from frozen private aggregate caches. It does not train a model, generate a depth draw, change ROI geometry, or write a manuscript. Numeric study results are derived from the inputs. The repository contains code, not the required private tables, ROI coordinates, H&E pixels, or generated figures.

Python 3.9 or later is required. Install the libraries in `requirements.txt` in an existing analysis environment. The code uses the noninteractive matplotlib Agg backend.

## Required caches

`IVDD_V3_ROOT` is the read-only V3 review-cache directory. Required relative paths are:

- `02_sigxtalk_official_new/formal_new/tables/`: `ALL_TESTED_NEW_AXES_WITH_PRS.tsv`, `AF_TO_NP_CANDIDATE_LR.tsv`, `CELLCHAT_ALL_SOURCE_TARGET_LR.tsv`, `LR_COMPLEX_COMPONENTS.tsv`, `TARGET_GENE_SELECTION_ALL.tsv`, `DONOR_GROUP_INPUT_COUNTS.tsv`, `CANDIDATE_TESTABILITY.tsv`, and `NEW_REPRESENTATIVE_OFFICIAL_FID_SPE.tsv`.
- `02_sigxtalk_official_new/formal_new/prs_outputs/`: `FID_SPE_RAW_UNIVERSE_OFFICIAL_NEW.tsv` and `FID_SPE_FILTERED_UNIVERSE_OFFICIAL_NEW.tsv`.
- `01_human_programs/FROZEN_K5_PROGRAMS.tsv` and `04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv`.
- `03_spatial_gse287499/`: `DEPTH_DISC_PAIR_COMPARISON.tsv`, `DEPTH_SENSITIVITY_VALUES.tsv`, `DISC_LEVEL_PROGRAM_AND_DEPTH.tsv`, and `CANDIDATE_GENES_PER_DISC_REGION.tsv`.
- `05_figures/source_data/Figure4ab_frozen_HE_vertices.tsv`.

`IVDD_V4_ROOT` is the read-only V4 cache directory. Required relative paths are:

- `01_sources_readonly/spatial_reference/`: `SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz`, `NP_PROGRAM_BY_REGION.tsv`, and `PROGRAM_GENE_SPATIAL_COVERAGE.tsv`.
- `01_sources_readonly/spatial_reference/roi_validation/`: `INDEPENDENT_MARKERS_BY_DISC_REGION.tsv` and `INDEPENDENT_MARKER_VALIDATION_BY_DISC.tsv`.
- `02_analysis/spatial_allprograms/NP_PROGRAM_BY_REGION.tsv`: the existing cached projection of every frozen program, including the previous summaries used for identity checks. This script does not produce that projection.
- `01_sources_readonly/spatial_reference/prepared48/<capture_id>/source_images/<source_image_filename>`: authorized local original H&E images. Filenames, pixel calibration and expected image SHA256 values come from the private ROI JSON. Pixels are unchanged; no contrast adjustment is made.

`IVDD_ROI_JSON` or `--roi-json` selects the private frozen ROI JSON. If omitted, the script expects `07_public_repo/configs/frozen_spatial_rois.json` inside `IVDD_V3_ROOT`. Its schema contains capture IDs, calibration, source-image metadata and disc/core polygons. Supply it from a private authorized cache. No ROI JSON is shipped in this folder. These authoring scripts do not use `IVDD_V2_ROOT`.

## Execution

Set the environment variables to your own local cache paths. `IVDD_OUTPUT_ROOT` must be separate from both input cache trees. CLI arguments override the equivalent environment variables.

```text
python scripts/authoring/v4/communication_spatial/plot_communication_spatial_v4.py --review-root <private-v3-review> --project-root <private-v4-cache> --roi-json <private-roi-json> --output-root <new-v4-render-directory>
```

The equivalent variables are `IVDD_V3_ROOT`, `IVDD_V4_ROOT`, `IVDD_ROI_JSON` and `IVDD_OUTPUT_ROOT`. Relative cache subpaths are intentionally stable; cache roots are configurable. Output is rejected if it overlaps an input cache. All generated source TSVs, figures and audit JSON remain local. Put local smoke output under this folder's ignored `.smoke/` directory, or in another private output directory. Never stage generated outputs.

This step requires completed frozen human/SigXTalk outputs, ROI assignment, the original cached fixed-depth draw, marker summaries and cached all-program spatial projection. Run it after those caches exist. To render the compact supplement next, point its `--main-tables-root` to `<new-v4-render-directory>/03_tables/communication_spatial`, while retaining the same read-only V4 cache root.

## Definitions and boundaries

LR complex keys remain separate from receptor component membership. Displayed receiver paths retain their receptor/SSC/target key. Path widths or PRS values do not establish a ligand-specific causal chain. Official raw-positive and filtered normalization universes remain separate. Missing PRS remains unavailable rather than zero. Main displayed path selection uses descending fixed-filter PRS and one path per SSC/target; LR ordering uses probability with lexical tie-breaking.

Spatial panels reuse actual cached bin centers, frozen polygons and the existing depth draw. Attenuation is `100 * (1 - mean_after / mean_before)` on the same eligible structure set. Structures are nested in captures; no animal- or genotype-level P value is generated. No model fitting, resampling or experimental raw values are loaded. Runtime source hashes and canvas/font checks are written under the selected private output root.
