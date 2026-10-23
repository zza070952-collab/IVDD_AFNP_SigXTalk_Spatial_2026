# IVDD fixed-program and official SigXTalk computation

This repository contains the study's own code for adult human fixed-program projection, genuine raw AF/NP model-input preparation, a new official SigXTalk run, frozen adult mouse spatial context, evidence integration and four computational figures. It does not contain raw matrices, individual-cell barcode tables, full HGNN pathway outputs, trained model weights, manuscripts, review archives or vendored third-party implementations. Only the fixed Top25 gene memberships, original lock and provenance are distributed; fitted cNMF spectra, loadings and other model weights remain private.

The computational source code is publicly available at https://github.com/zza070952-collab/IVDD_AFNP_SigXTalk_Spatial_2026. The initial public source snapshot is commit `01ab5ecfe71dabb833a9e03987891d7fde02aa76`; anonymous public access and all 162 source-file hashes were verified on 23 October 2026. This publication check is separate from the historical computation receipt. No licence has been granted for this own-code release yet; all rights are retained pending the study owner's choice. Third-party projects retain their own licences. Citation author names and a release identifier remain pending; no study author identity is invented.

Start with [environment setup](docs/ENVIRONMENT_SETUP.md), then copy the public configuration examples into private files and set real cache/environment locations. Public examples use only portable paths. The human module's [provenance](docs/FROZEN_K5_PROVENANCE.md), [human guide](README_HUMAN.md), [raw-selection guide](README_RAW_SELECTION.md), [official SigXTalk input guide](docs/SIGXTALK_INPUT_AND_EXECUTION.md) and [spatial guide](docs/SPATIAL_PUBLIC_REPRODUCTION.md) explain the component contracts.

```bash
bash run_all.sh --config config/private-pipeline.yaml --project-root ./analysis --from-raw
bash run_all.sh --config config/private-pipeline.yaml --project-root ./analysis --from-cached-inputs
python scripts/run_project.py --config config/private-pipeline.yaml --project-root ./analysis --from-cached-inputs --stage spatial
python scripts/run_project.py --config config/private-pipeline.yaml --project-root ./analysis --stage figures
python scripts/11_quality_audit.py --source-only
```

The stages are `all`, `human`, `raw-selection`, `sigxtalk`, `spatial`, `integrate`, `figures` and `audit`. The complete run retrieves public inputs, prepares the donor-resolved raw model input, generates adult human donor summaries, completes or verifies the single spatial depth check **before** the main formal model, completes the real technical R→Python→R closure and formal official SigXTalk run, integrates measured evidence, exports figures and audits actual outputs. Missing inputs, incomplete official closure and changed frozen inputs fail explicitly; a list of empty directories is never accepted as completion.

Raw mode supplies the full own-code route from the released matrices. It reuses the original locked K5 rather than silently refitting cNMF. An explicit historical fitting command is documented separately. Tu's released expression matrices contain published normalized values and are not called raw UMI counts. The SigXTalk input uses genuine GSE230808 integer count matrices and a fixed, independently recomputed donor/state selection.

Cached mode reuses actual cell-level QC/identity/scores and raw matrix files, with hashes and explicit skip receipts. It does not relabel historical SigXTalk predictions as a new run. Existing current model inputs are independently reselected and verified read-only, preserving their byte hashes once the new model has trained. Existing completed official models are verified against their frozen configuration, trained inputs, prediction hashes and pinned source. A completed HGNN can resume downstream PRS; it is not silently retrained.

The spatial sensitivity has one frozen100-UMI target, one seed, fixed ROI geometry and one draw. A valid completed check with matching configuration, ROI and protected-input hashes is reused without drawing again. Changed or partial depth outputs stop execution for review. No threshold search is performed. The [output contract](docs/EXPECTED_OUTPUTS.md) and [execution boundaries](docs/PIPELINE_BOUNDARIES.md) list what a completion claim does and does not establish.

The assembled workflow has technical tests and component checks. Unless an actual clean-machine raw-run receipt accompanies a release, it has not been independently validated from raw data to every final output. CI runs small structural fixtures and source audits only; it performs no GEO downloads, model training or research-data substitution. Figure export records actual source hashes; visual inspection remains a separate requirement.

Human donor expression is descriptive, Tu clinical-grade mapping remains unresolved, and AUCell/UCell scales are separate. Receptor-component joins do not fit ligand-specific downstream causality. SSC labels do not establish TF activity. Mouse genotype is confounded with capture, fixed cores are provisional and depth sensitivity weakens spatial specificity. No wet-lab outcome, Figure5 result or combined LR×PA×PRS score is generated.
