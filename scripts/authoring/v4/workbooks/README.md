# V4 static workbook export

These optional exporters consume existing V4 source tables. No input tables,
experimental observations, sequences, payload JSON, workbook or manuscript is
distributed here. Keep inputs outside the checkout and outputs under ignored
`analysis/`. The canonical computational runner does not call these scripts.

Install the Python dependencies from `envs/authoring_requirements.txt`. The JS
builder additionally requires an authorized environment providing
`@oai/artifact-tool`; that library is not vendored or installed by pip. JSON
preparation and native XLSX validation can run without that optional JS tool.

```bash
python scripts/authoring/v4/workbooks/prepare_workbook_payload.py \
  --project-root /private/v4-cache --v3-review-root /private/v3-review \
  --output-root ./analysis/v4-workbooks
node scripts/authoring/v4/workbooks/build_workbooks.mjs \
  --payload=./analysis/v4-workbooks/02_analysis/workbook_build/workbook_payload.json
python scripts/authoring/v4/workbooks/validate_workbooks.py \
  --payload ./analysis/v4-workbooks/02_analysis/workbook_build/workbook_payload.json
```

The preparer requires `03_tables/human`,
`03_tables/communication_spatial`, `08_review/communication_spatial`, and
`03_tables/experimental/internal/source_unresolved` in the private V4 cache,
plus `04_integrated/MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv` in the V3 review cache.
It uses the panel source maps to resolve the actual table paths. Missing tables
fail; this script does not generate scientific results or infer missing protocol
conditions. In particular, protocol/reagent tables must be separately supplied.
Use the corrected `../experiments/extract_inputs.py` for sequence pairs rather
than an obsolete derivative based on adjacent paragraphs.

The independent output root receives four static typed workbooks: computational
supplementary tables, main-figure source values, and two internal experimental
books. `--books=1,2` selects a subset in the JS builder; to validate a subset,
provide a private payload containing exactly those selected books. The default
validator requires every expected book. Full candidate evidence is also copied
to a private gzip; all generated files remain local and ignored.

The native validator compares every expected data cell and type, checks tables,
panes, formulas, merged regions, hidden sheets and public-book private paths.
Rendering and structural validation do not establish manual visual approval,
experimental authenticity or a clean-machine raw-to-final reproduction.
