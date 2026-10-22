# Reproducible human NP-program module

This is the human-data component of the study's own-code release. It covers public GEO retrieval, original discovery QC/identity/cNMF invocation, fixed-program projection to the eight Tu adult NP samples, authoritative gene-ID mapping, actual AUCell scoring, donor summaries, evidence tables and Figures1-2. Official cNMF, UCell and decoupler implementations are installed separately; no third-party source tree is presented as our own code.

## Inputs and environments

Use a new project root and a private copy of `config/human.example.yaml`. Paths are relative by default; executable locations and cache locations are configurable. Do not edit an old analysis directory in place. Raw-download receipts record URLs, sizes and SHA256; per-cell runtime data stay outside version control.

The Tu projection environment needs Python>=3.11 and `envs/human_projection_requirements.txt`; decoupler is pinned to2.2.0. Historical discovery reconstruction uses the separate versions in `envs/historical_discovery_requirements.txt` and the historical R/UCell dependencies described in `docs/FROZEN_K5_PROVENANCE.md`. Do not replace a shared base environment. A cached run can use genuine previously computed discovery metadata/scores; a raw reconstruction additionally needs the compatible R environment.

## Raw and cached modes

```bash
python scripts/02_crosscohort_np.py --config config/human.example.yaml --project-root ./analysis --from-raw
python scripts/02_crosscohort_np.py --config config/private-human.yaml --project-root ./analysis --from-cached-inputs
python scripts/10_make_human_figures.py --config config/private-human.yaml --project-root ./analysis
```

Raw mode downloads the exact public disease NP and Tu files, executes raw discovery QC/identity/UCell reconstruction, parses the Tu matrices, applies only unique authoritative symbol mappings, scores Tu with official AUCell, then generates donor/evidence tables and figures. It reuses the fixed distributed K5 rather than implicitly repeating60 NMF initialisations. An explicit historical cNMF replay is supplied separately in the provenance document.

Cached mode requires actual discovery per-cell metadata/scores. For Tu, either provide its published matrix files and let the scripts parse/rescore them, or set `tu.cached_tables_dir` to genuine existing Tu cell metadata, cell scores, measured-expression and coverage tables. Summary-only tables do not satisfy this requirement. Reuse is recorded as SKIPPED_REUSED_ACTUAL_CELL_LEVEL_INPUTS, not a new model run.

Each actual subprocess produces a log and exit code in `00_audit/HUMAN_MODULE_EXECUTION.json`. Failure leaves a PARTIAL_FAILED proof and preserved error/log rather than claiming completion. Frozen cNMF reuse is recorded separately from executed stages. Figures are exported as PNG/PDF/SVG, with plot-source TSVs and captions; visual QA remains pending until the actual images are inspected.

## Principal outputs

`01_human_programs/` contains FROZEN_K5_PROGRAMS.tsv, COHORT_AND_DONOR_MANIFEST.tsv, GEP2_DONOR_SUPPORT.tsv, PROGRAM_GENE_COVERAGE.tsv, OTHER_PROGRAMS_BRIEF.tsv, STAGE_A_METHODS_RESULTS.md and PROVENANCE_AND_SCORE_SCALE.md. Actual full cell scores/matrices are private runtime outputs, not public release files. Figures1-2 and panel-source TSVs are in `05_figures/Figure1_Figure2/`.

The original125 ranked gene memberships, lock and provenance are preserved in `resources/frozen_k5/`. Full gene spectra, consensus loadings and per-gene model coefficients remain private and are not redistributed. Unweighted UCell/AUCell projection needs the member lists rather than fitted coefficients. HGNC snapshots and public sample-level GEO metadata are small reference resources. No individual-cell barcode table, raw expression matrix, Seurat/RDS object, trained model, wet-lab data, internal review ZIP or manuscript is bundled.

## Interpretation and reproducibility boundaries

No grade association is claimed for Tu because the GEO-to-paper-patient crosswalk is unresolved. AUCell and UCell absolute scales are not directly compared. Working marker-supported identities are not definitive atlas labels; the GEP5-marker overlap and small populations remain limitations. Program projection does not establish communication or causality.

The raw scripts and historical scientific parameters are supplied, but a clean-machine historical reconstruction has not been rerun merely for this release. Cached-table/figure regeneration and technical tests must be distinguished from a validated full raw-to-final run. No100% reproducibility claim is made without an actual documented run. See the provenance and limitations documents.

Run lightweight technical checks with `python -m unittest discover -s tests`. These use tiny structural fixtures and frozen public gene resources, not mock research results or large GEO downloads. A project-wide workflow may invoke these module scripts from its own run_all/Makefile.

Licensing for the complete own-code repository is a maintainer decision; this component does not assign a licence to third-party code or data.

