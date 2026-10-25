# Compact communication and spatial supplement code

`build_compact_communication_spatial.py` renders compact supplementary prediction-boundary and spatial-coverage figures and writes static workbook payloads from private frozen caches. It does not create an XLSX workbook or manuscript. It performs no training, ROI fitting or additional depth draw. Study counts, coverage, depth results and counterevidence are read or computed at execution. No input tables, result payloads, manuscript text, ROI vertices or H&E images are bundled in this folder. The compact spatial drawing uses ROI geometry and actual bin centers and does not load H&E pixels.

Use Python 3.9 or later and the libraries listed in `requirements.txt`.

## Required caches and order

First complete the frozen input caches described in `../../v4/communication_spatial/README.md`, then run the V4 plotter. This compact script needs:

- The same V3 official axis, LR, complex-component, donor, frozen-program, target-selection, official raw/filtered metric and integrated-evidence tables.
- `02_sigxtalk_official_new/SIGXTALK_EXECUTION_PROOF.json` in `IVDD_V3_ROOT`. Training metadata and the upstream commit are read from this private receipt, not inserted as result literals.
- V3 `03_spatial_gse287499/DEPTH_DISC_PAIR_COMPARISON.tsv`, `DISC_LEVEL_PROGRAM_AND_DEPTH.tsv` and `CANDIDATE_GENES_PER_DISC_REGION.tsv`.
- V4 `01_sources_readonly/spatial_reference/SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz`, `roi_validation/INDEPENDENT_MARKERS_BY_DISC_REGION.tsv`, and `02_analysis/spatial_allprograms/NP_PROGRAM_BY_REGION.tsv`.
- The private ROI JSON selected by `--roi-json` or `IVDD_ROI_JSON`. It defaults to `IVDD_V3_ROOT/07_public_repo/configs/frozen_spatial_rois.json`; nothing in this folder supplies it.
- Main-panel source tables `Figure3_D_representative_paths.tsv`, `Figure3_E_official_fidelity_specificity.tsv` and `Figure3_C_receptor_target_max_PRS.tsv`. Select their directory with `--main-tables-root` or `IVDD_MAIN_TABLES_ROOT`. The default is `IVDD_V4_ROOT/03_tables/communication_spatial`; a fresh V4 run can instead supply its separate output table directory.

No V2 cache is used by this script. Input schemas match the named upstream cached tables; missing files cause a normal file error rather than synthetic substitution.

```text
python scripts/authoring/compact/communication_spatial/build_compact_communication_spatial.py --v3-review <private-v3-review> --v4 <private-v4-cache> --roi-json <private-roi-json> --main-tables-root <new-v4-render-directory>/03_tables/communication_spatial --output-root <new-compact-render-directory>
```

`IVDD_V3_ROOT`, `IVDD_V4_ROOT`, `IVDD_ROI_JSON`, `IVDD_MAIN_TABLES_ROOT` and `IVDD_OUTPUT_ROOT` provide equivalent defaults. The output directory must be separate from the input cache trees and main table directory. Source files are read only. Local smoke output may use the ignored `.smoke/` directory; generated payloads, tables, ROI drawings and figures must not be staged for publication.

## Fixed representative-path rule

Select the union of all displayed main path keys, all main metric keys, the winning max-positive-PRS axis for each displayed receptor/target pair, the top three fixed-filter PRS paths per discussed target, the highest-PA active axis without a returned PRS per discussed target, the top three ITGB1 receptor-role fixed-filter paths, and the prespecified ITGB1/HIF1A/PLOD2 key. Deduplicate exact receptor/SSC/target keys, retain every reason, then fill to **80** by descending fixed-filter PRS with receptor/SSC/target lexical ties. Stop with an assertion if mandatory keys already exceed that fixed display limit. These are representative paths, never described as all paths. The selection rule and number of mandatory keys are emitted at execution.

LR display retains the main top-18 probability keys plus all observed FN1 complexes. The original complex and component keys are kept separate. Candidate receptor, target and ligand roles are recorded separately. Raw-positive and filtered official fidelity/specificity values retain their distinct normalization universes. PA-active axes without PRS keep unavailable metrics and explicit reasons. Thresholds and display limits are analysis parameters, not biological significance tests.

The spatial output retains every supplied structure's ROI/QC, candidate and independent-marker coverage, eligibility, original/same-eligible/fixed-draw contrasts and the discordant comparator program. It calculates paired-set attenuation directly from cached depth contrasts and never redraws counts. Counterevidence is read from the integrated table, without a combined candidate score. No animal/genotype P value is generated. Experimental role strings describe design/release boundaries only; this script reads no experimental raw data and does not validate an experimental claim.

The outputs are private working artifacts: source TSVs, figure PDF/SVG/PNG, a typed JSON payload, input hashes and figure export checks. Only this code and its documentation belong in the public repository.
