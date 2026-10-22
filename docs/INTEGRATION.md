# Human module integration

Copy `scripts/`, `resources/`, `config/`, `envs/`, `docs/`, `tests/`, `README_HUMAN.md` and `.gitignore` into the new own-code repository. Merge shared directory contents deliberately. Do not publish `validation/`, Python bytecode, private configuration, cell-level caches, raw matrices or execution logs containing private filesystem locations. `PUBLIC_MODULE_MANIFEST.tsv` enumerates the intended files and SHA256 values.

Run from the new repository root. Set the new analysis root separately from the immutable old cache root. A private configuration may use absolute paths to read actual caches, but the distributed configuration contains portable paths only. The new scripts do not write into the cache inputs.

```bash
python scripts/02_crosscohort_np.py --config config/private-human.json --project-root NEW_ANALYSIS_ROOT --from-cached-inputs
```

Minimal private JSON configuration, with the uppercase placeholders replaced by actual paths:

```json
{
  "executables": {"human_python": "PYTHON_EXECUTABLE", "python_addons": "PYTHON_ADDONS_DIR"},
  "discovery": {
    "cached_metadata": "OLD_WORK/data_cache/discovery_cell_metadata.tsv.gz",
    "cached_scores": "OLD_WORK/data_cache/discovery_programs_UCell.tsv.gz",
    "af_donor_counts": "OLD_WORK/data_cache/AF_NP_DONOR_CELL_COUNTS.tsv"
  },
  "tu": {"cached_tables_dir": "OLD_WORK/01_crosscohort_NP"}
}
```

Remove `python_addons` if dependencies are installed normally. Discovery cache paths are independent; point them at the existing actual files rather than moving files to the illustrated layout. `tu.cached_tables_dir` is the actual completed old source-table directory, **not** the `Tu_sparse` directory alone. It must contain the two real per-cell tables `TU_CELL_METADATA.tsv.gz` and `TU_PROGRAM_AUCELL_CELL_SCORES.tsv.gz` plus actual coverage, gene mapping, expression sources, QC and AUCell provenance listed in `scripts/02_crosscohort_np.py`. The complete source-table cache is copied read-only as inputs, then donor summaries, supporting tables and both figures are newly generated from those measured data. Summary-only reuse is rejected.

The new output paths are `01_human_programs/` and `05_figures/Figure1_Figure2/`; the execution receipt is `00_audit/HUMAN_MODULE_EXECUTION.json`. Each subprocess exit code and log SHA256 is recorded. Frozen K5 is audited and reused, and cached cell-level computations are explicitly marked as reused. This cached run does not claim a fresh cNMF fit, raw QC run or AUCell fit.

`interpretation.review_decisions`, if supplied, must name an explicitly reviewed JSON bound to the exact new `PROGRAM_SCORE_DONOR_SUMMARY.tsv` SHA256. The supplied V2 decisions bind the genuine prior summary and will intentionally become `REVIEW_REQUIRED_INPUT_HASH_CHANGED` when method/provenance fields change in V3. Inspect the regenerated actual tables and rebind a review only after that inspection. The code does not infer CONSISTENT/MIXED labels from an arbitrary score cutoff.

For a new raw human reconstruction, use `--from-raw` and configure separate compatible discovery Python and historical R/UCell environments. The fixed K5 is still reused. Optional historical cNMF replay is explicitly requested through `scripts/01_discovery_np_cnmf.py --step cnmf --execute-historical-cnmf`; it cannot overwrite the distributed frozen lock. Raw discovery reconstruction and full model replay have not been rerun locally for this bundle.

The local private validation replays genuine completed small source tables through the new support-table and plotting code, rather than copying prior figures. Both PNGs were visually inspected for readable text, clipping, legends, grade caveats and separate UCell/AUCell axes. This validates table/figure regeneration, not the entire raw-to-final workflow. Technical fixtures cover sparse parsing, conservative aliases, frozen provenance and method-aware summaries.

One historical QC nuance matters when comparing newly reconstructed raw SigXTalk cell inputs: original discovery QC canonicalized symbols and summed duplicate symbols before counting nonzero genes. Thus equal total UMI can coexist with different feature counts if a new raw matrix retains duplicate feature symbols. The new raw-count universe should be documented; it should not be silently altered to force agreement with historical feature-count metadata.
