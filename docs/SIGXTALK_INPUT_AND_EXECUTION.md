# Configurable new official SigXTalk module

These are our orchestration and reporting scripts. They are ready to be run against real inputs but do not themselves establish computational success. Only newly produced completion statuses, nonempty outputs, actual logs and hashes support a new-run claim. The project owner manages downloads and the isolated environments.

## Real IVDD input

The configured `input_dir` contains:

- `counts.mtx.gz`: nonnegative integer sparse MatrixMarket, genes×selectedcells.
- `genes.tsv`: one unique original human gene symbol per matrix row, no header.
- `meta.tsv`: one row per matrix column, in the exact same order, with unique `cell_id` and `donor_id`, `gsm`, `geo_series`, `condition`, `bridge_group`. Expected groups are AF_stromal, NP_GEP2_high and NP_GEP2_low, all GSE230808/diseased. Do not create counts from saved scores or predicted outputs.
- An explicit configured pairing TSV with `donor_id,AF_gsm,NP_gsm`. Every selected cell is checked against its actual donor/tissue GSM.
- Optional frozen threshold TSV with `donor_id,q25,q75` and an optional `GEP2_top25_UCell` metadata column allow a direct saved-score state check. If these optional score fields are not supplied, state derivation must be documented by the upstream raw-selection script and its saved score/label manifest, not inferred by this module.

The original selected design uses six adequately supported paired donors, seed20261017 and cap150/group/donor. It historically selected2,520 cells, but the new R wrapper does not force old count/target/pathway totals. It checks the real current cells and records their actual dimensions and hashes. A larger all-common raw GEX gene universe is legitimate when the original objects had filtered/collapsed features; record the actual feature-identity and duplicate-handling audit separately. Preserve raw genes and original IDs. Authoritatively justified global aliases can be handled by the upstream preprocessing audit, but target genes must not be inserted or renamed merely to obtain nominated paths.

Raw input inspection confirmed CCN2 is present in the current selected raw feature universe and CTGF is absent, so no CCN2 alias rename is needed for this main input. The official human TFT prior includes CCN2. Candidate testability exports preserve raw presence, sender/receiver count-detection fractions, target-selection status and official-prior filtering. Seurat's own feature-character normalization is recorded as a complete identity map.

## Configuration and commands

Copy `sigxtalk.config.example.json` into the project's `configs/` directory. Its relative `project_root=..` is resolved against the config file's parent by both Python and R. All other relative data/source/library paths are resolved against that root. Set actual private Rscript/Python paths in the private execution config; the public example uses generic executables. Put the wrapper scripts and `sigxtalk_upstream_lock.json` together under `scripts/`.

```
bash scripts/04_run_sigxtalk_hgnn.sh configs/sigxtalk.json smoke all
bash scripts/04_run_sigxtalk_hgnn.sh configs/sigxtalk.json formal all
python scripts/06_sigxtalk_reports.py --config configs/sigxtalk.json --mode formal --plot-official-metrics
python scripts/10_plot_new_sigxtalk_figure3.py --config configs/sigxtalk.json
```

Use `smoke.input_kind=pbmc_rds` when the official real PBMC RDS is available. Otherwise the declared `ivdd_raw` mode is a real-IVDD technical smoke and must not be called a completed official PBMC demo. PBMC uses its actual original identities and CD14_Mono receiver, the official positive DEG rule, then a predefined first50 p_adj/logFC/gene-ordered targets. The formal IVDD rule remains200 variance-selected receiver targets, >=.1 detection, with MT-/RPL/RPS excluded and no candidate panel insertion. Parameters are chosen once before results inspection. PBMC input file hash and actual target/cell counts are recorded.

The smoke must complete real R preparation→official HGNN→official R PRS with nonempty PRS before formal training is allowed. The driver verifies the same frozen config hash. It copies the unmodified official main.py into the private mode directory, captures actual stdout/stderr and exit codes, checks all epoch-end evidence and hashes the complete prediction file. It refuses silently overwriting an existing new prediction CSV. Resume a failed downstream PRS with `formal prs` rather than retraining completed HGNN.

The initial two-epoch PBMC subset completed new R preparation and official HGNN but produced zero PA>.75 paths, so PRS correctly did not complete. Preserve that undertrained attempt and its outputs. A separate `smoke.directory_name=technical_smoke_100epochs` and `smoke.epochs=100` applies the official standard training length to the same frozen PBMC context, genes and PA threshold; it changes neither the formal IVDD design nor biological evidence. The wrapper accepts predeclared integer smoke lengths from1 to100 and resolves the formal prerequisite from the configured smoke directory. A100-epoch technical subset is still a technical subset, not a full published PBMC reproduction, and its result must be reported from actual completion artifacts.

## Dependencies and provenance

Official package source: <https://github.com/LithiumHou/SigXTalk>, commit `4d8055c6c9539ff9e6875d0e731d933cb2dd69f9`. Critical source files are verified against original Git blob IDs; installed R function formals/bodies are compared with the pinned checkout. Human priors are imported from SigXTalkR and object hashes are recorded. Actual R sessions and new input/output hashes are private run audit outputs.

Author paper-reproduction scripts: <https://github.com/LithiumHou/SigXTalk_scripts>, current verified commit `e5a15be86a254dda56e0f45ae082910325cca4cf`. Its older HNSCC script uses `Prepare_Input_New`/`Rec_To_TFTG` and author-specific paths, so those calls are not blindly mixed with the packaged current API.

R requires SigXTalkR, Seurat>=5.1, SeuratObject, CellChat, ranger, tidymodels, tidyfst, Matrix and listed DESCRIPTION imports. CellChat v2.1.2 is verified tag commit `a0d3b2d231d46c8787177fffeac908270c253747`; presto v1.1.0 is `b5df6ee6097eb62522f2e557aa93ac21bab2d05f`. tidyfst tagv1.7.5 points to `46df63f303c1325085964e5241174782211b2aa6`, whose actual DESCRIPTION reports1.7.3; record the installed value rather than relabelling it. Its imports are data.table, fst and stringr. Install only in the new isolated environment/project library and do not use upstream `upgrade=always` to mutate shared libraries.

Python package constraints are >=3.8,<3.9, pandas2.0.3, scipy1.10.1, scikit-learn1.3.0, numpy1.24.3, dhg0.9.4 and torch1.13.1. DHG requires torch<2. CPU is officially supported. The existing compatible Python environment can be referenced read-only. The separate plotting/report environment needs pandas/numpy/matplotlib. See [resource patch audit](SIGXTALK_RESOURCE_PATCH.md) for the explicit four-worker RNG/thread wrapper, unchanged mathematical definition and actual runtime diff.

## Output truth and public-code boundary

Private `formal_new/` and `technical_smoke/` outputs are separate. Their new Seurat/CellChat RDS, receiver expression RDS, complete original HGNN CSV/weights, raw+filtered official PRS, two declared official fidelity/specificity universes, row-level missingness and metadata are never inferred from the October17 tables. Full original LR context is preserved before AF-only scoping. These objects/full results remain private project artifacts and are excluded from the public code repository.

Publish only the explicit own-code allowlist supplied with this module. Do not publish `upstream_inspection/`, original downloaded repo snapshots, third-party source copies, runtime deparsed PRS functions, private execution configs, new logs containing machine paths, raw matrices, RDS/h5ad, complete pathway tables, model weights or cell-barcode metadata. The source lock and small attributed patch are public-compatible provenance, not a vendored implementation.

Interpretation remains bounded: a receptor-complex/component-key join links AF-source LR context to receiver predictions, not a ligand-specific fitted chain. The raw TF column is SSC/intermediate; PA, PRS, fidelity and specificity are not clinical p-values or mechanistic proof. GEP2-defined receiver/target coverage is not independent program validation. Report unsupported candidates, low coverage and current missing metrics without inserting old positive rows.
