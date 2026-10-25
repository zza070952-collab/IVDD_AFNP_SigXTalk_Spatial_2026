# Compact supplementary workbook export

These scripts export private prepared source tables into six supplementary-table
sheets and three experimental-source sheets. No real payload, Ct/OD observation,
oligo sequence, experimental protocol document or workbook is included in Git.
All experimental exports retain `INTERNAL_REVIEW_SOURCE_UNRESOLVED`.

Python dependencies are in `envs/authoring_requirements.txt`. The optional JS
builder requires an authorized `@oai/artifact-tool` runtime. The library is not
vendored; no public npm availability is assumed.

Run the compact human and communication/spatial builders first, using the same
private revision output root. They create `human/*PAYLOAD*.json` and
`communication_spatial/*PAYLOAD*.json`. Then prepare the external experiment
tables and assemble the private workbook payload:

```bash
python scripts/authoring/compact/workbooks/prepare_experiment_payload.py \
  --project /private/v4-cache --revision ./analysis/compact \
  --protocol /private/supplied-protocol.docx
python scripts/authoring/compact/workbooks/assemble_workbook_payload.py \
  --revision ./analysis/compact --work-dir ./analysis/compact/workbook-build
node scripts/authoring/compact/workbooks/build_compact_workbooks.mjs \
  --payload=./analysis/compact/workbook-build/COMPACT_WORKBOOK_PAYLOAD.json
python scripts/authoring/compact/workbooks/validate_compact_workbooks.py \
  --project /private/v4-cache --revision ./analysis/compact \
  --payload ./analysis/compact/workbook-build/COMPACT_WORKBOOK_PAYLOAD.json \
  --original-workbook /private/supplied-observations.xlsx
```

The experiment preparer reads `03_tables/experimental/internal/source_unresolved`
under `--project`. It requires `RAW_EXTRACT.json`, batch qPCR/CCK8 summaries,
matched blank values, contrasts, experimental design and optional external
`SUPPLIED_PROTOCOL_CONDITIONS.tsv`. It extracts both oligo strands from the same
original protocol paragraph. No unpublished protocol sentences are embedded in
the source. `--experimental` on the assembler omits the five computational
sheets for a private experiment-only check; this is not a full compact package.

The JS output is private `INTERNAL_REVIEW/`, with rendered previews under
`audit/workbook_render/`. The validator requires every book listed in the payload,
compares every payload data cell, rejects duplicate record keys, and compares
every original raw observation back to the external supplied workbook. It checks
filters, panes, formulas, merged cells, hidden rows/columns and private paths.
It does not authenticate source provenance, validate instrument runs, or certify
manual visual inspection. Never stage its payloads, logs, outputs or source data.
