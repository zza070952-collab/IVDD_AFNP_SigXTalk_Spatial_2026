# Compact human supplement selection

`build_compact_human.py` reuses frozen human summaries to rebuild one compact Supplementary Figure S1 and select two static long-format tables: `S1_Cohorts_QC` and `S2_Fixed_NP_Programs`. It writes TSV and typed JSON payloads for a separate workbook assembler; it does not author an XLSX file. No result payload is included in this code directory.

Panel A displays donor QC retention. Panel B displays actual NP/immune coverage and comparator availability, preserving absent comparators as NA. Panel C displays the predeclared illustrative GEP5 matrix/stress-member subset. This subset is not an enrichment test or a new program selection. All frozen programs and ranked members are retained in the table selection. Program status, observed donor counts, score comparisons, coverage and missing features come from the supplied frozen inputs. No new statistical test, score, fit, embedding or grade analysis is performed; AUCell and UCell magnitudes are not compared.

## Run

From the repository root:

```sh
python scripts/authoring/compact/human/build_compact_human.py \
  --source-root /path/to/frozen-review-cache \
  --v4-root /path/to/read-only-v4-cache \
  --output-root analysis/human-authoring-compact
```

`IVDD_REVIEW_ROOT` and `IVDD_V4_ROOT` may replace the corresponding arguments. `--eligibility-table` can replace the V4-root requirement. The output directory must be outside both read-only input roots. Use a new output directory. `analysis/` is ignored by this repository; do not force-add output. The exact tested plotting dependencies are in `requirements.txt`; Python 3.8.5 was used locally.

## Required external caches

The source root must contain these frozen review files:

```text
01_human_programs/
  COHORT_AND_DONOR_MANIFEST.tsv
  FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv
  FIGURE2_DONOR_PROGRAM_SOURCE.tsv
  FROZEN_K5_PROGRAMS.tsv
  PROGRAM_GENE_COVERAGE.tsv
  TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv
  ALL_PROGRAMS_BRIEF.tsv
  GEP2_DONOR_SUPPORT.tsv
05_figures/Figure1_Figure2/Figure1a_donor_tissue_counts.tsv
02_sigxtalk_official_new/formal_new/tables/DONOR_GROUP_INPUT_COUNTS.tsv
07_public_repo/resources/metadata/SWAHN_VERIFIED_DONOR_PAIRING.tsv
```

The V4 root must contain `00_inventory/RECOMPUTED_DONOR_ELIGIBILITY.tsv`. An explicit eligibility-table argument can provide that same schema without a V4 root. The input schemas are the original cached TSV schemas: cohort/donor/GSM and QC counts; donor/identity cell counts; method-labelled program medians; original gene membership/rank, mapped genes and coverage; donor/identity gene expression/detection; frozen qualitative program review; actual comparator medians; verified donor-pair metadata; and actual donor-quartile eligibility. These inputs are not bundled and are not fetched automatically. No raw cell matrix, experiment Ct/OD well or manuscript is read.

## Outputs and limits

The script emits `figures/` PDF/SVG/300-dpi PNG and panel-source TSVs, `tables/` S1/S2 TSVs, a static `HUMAN_WORKBOOK_PAYLOAD.json`, data dictionary, input-hash manifest, figure contract, calculated facts and structural QA. `value` and `frozen_rank` payload fields are numeric or null; null is not a measured zero. Full fixed member lists and GEP2 donor summaries are retained, together with source-derived missing/mixed/reversed evidence. Operational identity support does not certify histological purity.

Full manuscript paragraphs, the original supplementary figure legend and manuscript reference-remapping text were excluded from this public entry point. The figure labels are concise plotting annotations, not a distributed manuscript. This code does not change the original main figures or original cached outputs.

`PASS_AUTOMATED_AWAITING_VISUAL_QA` checks only structural bounding boxes. `actual_visual_review` remains `PENDING` for every new run. A cached-input smoke run and source-table numeric comparison do not establish a new manual visual review, experimental provenance, clean-machine reconstruction of the upstream analysis or formal release approval. The exact dependency and source-manifest receipts are local output, not proof that every required cache has been publicly deposited.
