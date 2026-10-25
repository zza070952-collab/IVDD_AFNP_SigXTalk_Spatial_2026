# V4 and compact-supplement authoring

These optional tools assemble figures, static table selections and documents from supplied frozen inputs. They are separate from the canonical computational dispatcher in `run_all.sh` and `scripts/run_project.py`. That dispatcher has no wet-lab stage. Publishing authoring code does not add experimental records or certify a full upstream analysis, an experiment or a submission.

This is a code-only addition. Experimental Ct/OD measurements, instrument exports, supplied oligo sequences, manuscript paragraphs/legends, reference-replacement plans containing manuscript text, generated table payloads, figures, workbooks and review packages are not distributed. The small fixed-program definitions already in `resources/frozen_k5/` are method resources, not model weights or a replacement for the required caches. Scripts read authorized private inputs at runtime; they do not download a complete V4 authoring dataset from this repository.

## Entry points

| Directory | Purpose and entry points | Required external inputs |
|---|---|---|
| [V4 human](../scripts/authoring/v4/human/README.md) | `build_human_figures.py`: human main figures, earlier full supplementary figures, source tables and structural QA | Frozen V3 human summaries, donor eligibility, frozen four-box design-source table |
| [V4 communication/spatial](../scripts/authoring/v4/communication_spatial/README.md) | `plot_communication_spatial_v4.py`: main and earlier supplementary plots, source tables and audit records | Completed official SigXTalk tables, frozen ROI/depth/marker summaries, cached spatial program projections and authorized original H&E pixels |
| [V4 experiments](../scripts/authoring/v4/experiments/README.md) | `extract_inputs.py`, `analyze.py`, `plot.py`; `sequence_parser.py` supports exact protocol transcription | Separately supplied private workbook and protocol; no measurements or sequences are embedded in code |
| [V4 workbook tools](../scripts/authoring/v4/workbooks/) | `prepare_workbook_payload.py`, `build_workbooks.mjs`, `validate_workbooks.py` | Complete private human/communication/spatial authoring tables and any requested unresolved experiment records; optional authorized spreadsheet runtime |
| [V4 spatial-cache bridge](../scripts/authoring/v4/cache/README.md) | `project_frozen_spatial_programs.py`: existing fixed-member spatial scorer exposed for all frozen programs | Prepared sparse spatial counts, native-feature coverage, ortholog mapping and frozen ROI metadata |
| [Compact human](../scripts/authoring/compact/human/README.md) | `build_compact_human.py`: compact human figure and two typed static table selections | Frozen human summary tables, donor-pair metadata and eligibility |
| [Compact communication/spatial](../scripts/authoring/compact/communication_spatial/README.md) | `build_compact_communication_spatial.py`: compact prediction-boundary/spatial-coverage figures and table selections | Frozen official/spatial inputs plus the V4 plotter's displayed-path source tables and actual run receipt |
| [Compact workbook tools](../scripts/authoring/compact/workbooks/) | `prepare_experiment_payload.py`, `assemble_workbook_payload.py`, `build_compact_workbooks.mjs`, `validate_compact_workbooks.py` | Private compact module payloads and the original experimental workbook/protocol where required; optional authorized spreadsheet runtime |
| [Common authoring utilities](../scripts/authoring/common/) | PDF assembly, PDF render/audit, explicit reference replacement, DOCX formatting/export, shared checks and generic synthetic smoke | External layout/style/reference plans, authorized text/figure inputs and actual manual-review records |

Module READMEs and `--help` describe their precise flags and fixed relative schemas. These are cache-specific authoring entry points, not generic importers for arbitrary datasets. Missing files stop execution; the code does not fabricate observations or fill unavailable scientific values with zeros.

## Input roots and execution order

Use separate read-only input roots, a private working directory and a new independent output directory. Public examples below are placeholders for paths selected by the operator. Do not write into an original cache or force-add generated output. The repository ignores `analysis/`, while some utilities also provide a local ignored `.smoke/` directory.

The human scripts accept `IVDD_REVIEW_ROOT` and `IVDD_V4_ROOT`. Communication/spatial scripts use `IVDD_V3_ROOT`, `IVDD_V4_ROOT`, `IVDD_ROI_JSON` and `IVDD_OUTPUT_ROOT`; their explicit CLI flags override those defaults. The V3 review root is the frozen analysis **review** directory, not this source repository. ROI JSON, coordinates and images are external private inputs. The compact communication/spatial script additionally needs `IVDD_MAIN_TABLES_ROOT` or `--main-tables-root`, pointing to the actual V4 plotter's source tables.

1. Complete or obtain authorized upstream frozen caches using the documented computational pipeline. The spatial-cache bridge can generate the required all-program projection from its prepared sparse inputs. It preserves fixed membership, coordinates and measured zeros; it does not refit ROIs or repeat the frozen depth draw.
2. Run the V4 human and communication/spatial figure builders against those caches, choosing a separate private output root. They can populate their distinct subdirectories beneath the same new output root. The human design-source input is explicit so cross-module counts do not become source-code constants.
3. Invoke the experiment module only when the private workbook and protocol are available and analysis of them is authorized. It is optional and is not called by the canonical computational dispatcher.
4. Run the compact human and communication/spatial selectors. Preserve their expected `human/` and `communication_spatial/` subdirectories beneath a private compact-revision root. The compact communication/spatial selector consumes actual V4 main-panel sources, not guessed representative paths.
5. Prepare any requested workbook payloads, build them only in an authorized spreadsheet runtime, and compare saved values/types against their source payloads. Payload and workbook products remain private. A compact three-file review set can include unresolved experimental material; its file count is not a provenance decision.
6. Assemble or format documents only from externally supplied authorized content, then render the actual exported PDF and complete manual visual review. Code does not ship manuscript prose or a pre-completed review record.

## Optional experiment authoring

The optional experiment entry points are independent commands:

```sh
python scripts/authoring/v4/experiments/extract_inputs.py \
  --workbook PRIVATE_INPUT.xlsx --protocol PRIVATE_PROTOCOL.docx \
  --project analysis/private-experiment-review
python scripts/authoring/v4/experiments/analyze.py \
  --project analysis/private-experiment-review
python scripts/authoring/v4/experiments/plot.py \
  --project analysis/private-experiment-review
```

They preserve source row identifiers and missing/undetermined states, apply the declared batch-level calculations, and can generate conditional Figure 5 and related internal source tables. They do not create or recover missing laboratory observations. All such outputs retain `INTERNAL_REVIEW_SOURCE_UNRESOLVED`; arithmetic agreement does not verify instrument provenance, assay QC, protein-level effects or a causal mechanism. Original measurements, protocol sequences and generated results are not committed. The module README defines its input worksheets and statistical-method contract without distributing observations.

## Workbook tools and runtime

Workbook payload preparation is Python code. V4 preparation accepts `--project-root`, `--v3-review-root`, `--output-root` and optional `--work-dir`; it expects the complete documented private authoring cache. Compact preparation reads module payloads below `--revision` and writes its build specification under `--work-dir`. Experiment-specific preparation additionally requires an external `--protocol` and the private project tables.

The `.mjs` builders use `@oai/artifact-tool`. This is an **optional, separately authorized runtime dependency**. The repository does not distribute it or claim that public pip/NPM installation is available. Use an existing runtime that is actually entitled to load it and resolve its modules. If that runtime is unavailable, Python payload preparation and source inspection remain possible; the JavaScript workbook build is unavailable. No alternative package is silently substituted or installed.

After preparing a private specification, the builder interfaces are:

```sh
node scripts/authoring/v4/workbooks/build_workbooks.mjs \
  --payload=PRIVATE_WORK/workbook_payload.json
python scripts/authoring/v4/workbooks/validate_workbooks.py \
  --payload PRIVATE_WORK/workbook_payload.json

node scripts/authoring/compact/workbooks/build_compact_workbooks.mjs \
  --payload=PRIVATE_WORK/COMPACT_WORKBOOK_PAYLOAD.json
python scripts/authoring/compact/workbooks/validate_compact_workbooks.py \
  --project PRIVATE_V4_CACHE --revision PRIVATE_COMPACT_ROOT \
  --payload PRIVATE_WORK/COMPACT_WORKBOOK_PAYLOAD.json \
  --original-workbook PRIVATE_INPUT.xlsx
```

These commands are examples for an available authorized runtime, not proof that inputs or generated workbooks are public. Validate static cell types, actual missing-value states, source completeness and final saved values. Render the saved workbook ranges and inspect them. Successful serialization or cell comparison does not replace visual QA or certify native application behavior that has not been checked.

## Common document and visual-QA utilities

| Utility | Interface and evidence boundary |
|---|---|
| `assemble_pdf.py` | `--layout --output --receipt`; assembles external Markdown/text and existing figure PDFs according to an external layout. It does not supply scientific prose. |
| `pdf_qa.py audit-render` | `--pdf --output-dir`, optional minimum font, expected page count, DPIs and review template. Checks structure and renders pages; the review template remains pending. |
| `pdf_qa.py validate-review` | Binds an externally completed manual review to the actual PDF and automated report hashes. It rejects stale or incomplete records rather than inventing a PASS. |
| `replace_references.py` | `--input --output --plan --receipt`, optional `--dry-run`; applies an explicit private reference plan and checks preservation boundaries. No manuscript-specific plan is shipped. |
| `format_docx.py format` | `--input --output --styles --receipt`; applies an external style specification while preserving document content and media boundaries. |
| `format_docx.py word-pdf` | `--input --output --receipt`; optional native Microsoft Word export on a supported Windows host with an available Word installation and `pywin32`. Export success still requires actual PDF QA. |
| `authoring_checks.py` | Shared hashing, package/text inspection and plan-validation helpers; no manuscript payload. |
| `smoke_synthetic.py` | Generic local fixtures for assembly, formatting, reference preservation and rejection of stale/invalid review claims. These fixtures are not study data. |

Python plotting/document dependencies are recorded in the per-module requirements. The original upstream analysis has separate R/Python requirements and is not replaced by a plotting environment. No additional Word installation, account access or third-party entitlement is implied by publishing these utilities.

## Validation recorded on 25 October 2026

Local authoring checks used the available authorized caches and runtime, not a newly provisioned clean machine. Human builders were syntax/help checked and executed against frozen cached summaries; generated panel-source values and selected-table numeric/NA fields were compared with their frozen authoring baselines. Communication/spatial, experiment and workbook checks use their available private inputs and retain local receipts. The spatial-cache bridge was syntax/CLI checked without replaying the large sparse input. Common document utilities use generic synthetic fixtures. These are distinct scopes; a syntax check is not reported as a real-input execution.

The repository's current technical suite is run with:

```sh
python scripts/11_quality_audit.py --source-only
python scripts/authoring/audit_public_release.py
python -m unittest discover -s tests -v
python scripts/authoring/common/smoke_synthetic.py
```

The existing CI runs the source audit and unittest discovery, including the new publication-scope regressions. The separate publication-source gate inspects the actual Git index; the common smoke uses generic synthetic document/figure fixtures locally with its optional dependencies. The suite may gain tests; no fixed test count is a release criterion. CI fixtures do not download GEO inputs, train biological models, replace research observations or authenticate experiments. A local cache smoke run does not demonstrate a complete independent GEO-to-all-final reconstruction. Structural figure QA remains separate from actual manual viewing, which every new export requires. Private runtime receipts and generated previews are not committed as distributed study results.

No future commit identifier, new archived DOI or unperformed anonymous-access verification is claimed here. The study owner's own-code reuse licence remains unselected; third-party licences and access conditions continue to apply.
