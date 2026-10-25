# Main-figure renumbering and source-table composition

These authoring scripts reuse frozen outputs. They do not rerun cNMF, UCell,
AUCell, SigXTalk, CellChat, spatial inference, or experimental statistical models.
Plotting summaries and error bars are evaluated with the original formulas from
the supplied tables. No data tables, experimental measurements, manuscripts,
finished figures, reference images, or private result hashes are distributed here.

| New figure | Operation |
| --- | --- |
| Figure 2 | Replot old Figure 1 data panels B–E as A–D; omit the old schematic. |
| Figures 3–5 | Copy old Figures 2–4 in all four formats without changing bytes. |
| Figure 6 | Replot old Figure 5 data panels B–G as A–F; omit the old schematic. |

The fixed cohort labels, group definitions, marker order, axis limits, colors,
time points, annotations, and layout are inherited from the frozen authoring
code. These are reproduction scripts for that design, not a general-purpose
analysis interface. Do not silently substitute a different cohort or design.

## Dependencies and commands

Use Python 3.12 with NumPy, pandas, SciPy, Matplotlib, Pillow, and PyMuPDF. Exact
historical environment information belongs in the parent reproducibility
documentation. Set `MPLCONFIGDIR` to a writable local cache if needed.

From this directory, with external inputs stored outside the repository:

```sh
python replot_figures2_6.py --source-root /path/to/read_only_v4 --output-root /path/to/local_replot --figures both
python copy_retained_figures.py --source-root /path/to/read_only_v4 --output-root /path/to/local_delivery
python validate_legacy_data_panels.py --source-root /path/to/read_only_v4 --legacy-code-root /path/to/trusted_legacy_code --output /path/to/local_audits/legacy_oracle.json --trust-legacy-code
python -m unittest discover -s . -p 'test_*.py'
```

Source and output must be disjoint trees. Existing output symlinks/junctions are
rejected. Inputs are hashed before and after the operation; output hash receipts
are generated locally and must not be added to the public repository when they
refer to private data. The copy command preflights all twelve required files,
then confirms each input and output digest plus a final input recheck. Existing
same-named files inside the chosen output directory can be replaced.

## External input contract

`--source-root` has the following frozen table paths. Only the listed columns
are consumed; no example rows or real values are bundled.

| Relative path | Required columns |
| --- | --- |
| `03_tables/human/Figure1B_source.tsv` | `donor_id`, `AF_retained_cells`, `NP_like_cells`, `NP_high_input`, `NP_low_input`, `communication_included` (boolean) |
| `03_tables/human/Figure1C_source.tsv` | `donor_id`, `identity`, `n_cells` |
| `03_tables/human/Figure1D_source.tsv` | `identity`, `gene`, `fraction_detected`, `mean_log1p_published_expression` |
| `03_tables/human/Figure1E_source.tsv` | `donor_id`, `identity`, `n_cells` |
| `03_tables/experimental/internal/source_unresolved/QPCR_BATCH_DELTACt.tsv` | `Gene`, `Group_ID`, `Batch_ID`, `fold_vs_same_batch_G1` |
| `03_tables/experimental/internal/source_unresolved/QPCR_CONTRASTS_AND_CI.tsv` | `gene`, `contrast`, `fold_ratio`, `fold_CI95_low`, `fold_CI95_high`, `p_adjusted`; `SOURCE_STATUS` when available must retain the unresolved marker |
| `03_tables/experimental/internal/source_unresolved/QPCR_BATCH_CONTRASTS.tsv` | `gene`, `contrast`, `Batch_ID`, `fold_ratio` |
| `03_tables/experimental/internal/source_unresolved/CCK8_BATCH_SUMMARY.tsv` | `Group_ID`, `Batch_ID`, `actual_read_hours_after_reseeding`, `ratio_vs_same_batch_time_G1` |

Numeric columns must retain their frozen scale, reference group, and missingness.
The plotting code preserves published-expression and donor order, logged fold
ratios, confidence intervals, and adjusted probabilities. It does not estimate
new experimental contrasts. Experimental files remain privately supplied inputs;
their provenance boundary is `INTERNAL_REVIEW_SOURCE_UNRESOLVED`. A successful
plot or validation does not authenticate the underlying experimental source.

The byte-preserving copy requires these stems under `04_figures/`, each with
`.pdf`, `.svg`, `.png`, and `.tiff`:

- `Figure2/Figure2_Frozen_Programs_CrossCohort`
- `Figure3/Figure3_Predicted_Communication`
- `Figure4/Figure4_Anatomy_and_Depth`

## Legacy oracle and trust boundary

The optional oracle accepts `build_human_figures.py` and `plot_experiments.py`
from `--legacy-code-root` (default `SOURCE_ROOT/07_public_code_delta`). The two
legacy scripts are not bundled in this directory. It extracts selected helper
functions and data-panel statements from the exact frozen layout, executes them
in memory, and compares artists with the replacement plots. The whole legacy
`main`, saving routines, manuscript generation, and supplementary exports are
not invoked. Different versions of those scripts may fail the extraction.

**Extraction is not a security sandbox.** The selected Python statements execute
with the caller's privileges. Inspect and trust both local files before passing
the required `--trust-legacy-code` opt-in. Never use downloaded or unreviewed
code merely to obtain an oracle result.

The comparator checks data coordinates, matrix values, bars, error segments,
axis limits/labels, colors, markers, and line styles, while deliberately allowing
the requested panel geometry and lettering changes. Automated checks accompany,
but do not replace, actual-size visual review of the exported PDFs and TIFFs.
