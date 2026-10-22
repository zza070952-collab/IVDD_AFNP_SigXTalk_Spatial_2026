# Actual output contract

| Folder | Required measured result/proof |
|---|---|
| `00_audit` | actual public-download hashes, per-step logs/exit receipts, `PROJECT_EXECUTION_STATUS.json`, `QUALITY_AUDIT_STATUS.json`, `PIPELINE_OUTPUT_MANIFEST.tsv` |
| `01_human_programs` | fixed125 genes, coverage and mapping, donor scores, `GEP2_DONOR_SUPPORT.tsv`, other-program context and actual Tu candidate donor/identity expression |
| `02_sigxtalk_official_new/input` | genuine integer gene×cell counts, genes and exactly aligned metadata at the configured selected-input subdirectory; independent selection/barcode/source audits |
| `02_sigxtalk_official_new/<technical directory>` | separate real technical input, official HGNN and R PRS closure receipts, without formal biological interpretation |
| `02_sigxtalk_official_new/formal_new` | new official R preparation,100epoch HGNN, complete prediction hash, official PRS/missingness, fidelity/specificity status and `NEW_OFFICIAL_FINAL_STAGE_STATUS.json` |
| `03_spatial_gse287499` | one frozen depth check, per-disc retention/contrasts, score–UMI description, candidate context and protected-input hashes |
| `04_integrated` | candidate evidence matrix, exact source context, current wet-lab nomination decision, limitations and `INTEGRATION_STATUS.json` |
| `05_figures` | Figures1–4 as PNG/PDF/SVG, measured panel-source TSVs, legends and render receipts; no experimental Figure5 |
| `06_manuscript` | computation Methods/Results, evidence limitations, data/code availability and source-structure drafts derived from actual completed inputs |

Only `DONE_NEW_FORMAL_INTEGRATION` establishes complete new-formal integration. `PARTIAL_NEW_B_PENDING` is a draft state and cannot satisfy top-level `all`. Official core completion with missing fidelity/specificity remains labelled by its actual B `run_class`; the dispatcher does not fabricate those metrics.

The output manifest contains hashes of scientific tables, figures and documents; model weights and raw count caches remain separate runtime objects. Files containing barcodes, full pathways, private machine paths or large expression objects are not committed wholesale. The public release uses explicit source allowlists.

GitHub availability is independent of computational completion. Until a real public URL and commit are verified, the code-availability draft must retain `BLOCKED_UNVERIFIED`; no placeholder URL or fake release is supplied.
