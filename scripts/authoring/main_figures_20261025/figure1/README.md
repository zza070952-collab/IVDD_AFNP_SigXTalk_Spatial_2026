# Original study-design schematic

`build_figure1.py` preserves the completed original vector drawing: 180 x 220 mm,
editable SVG/PDF text and PNG/TIFF exports. Mini matrices, spatial bins and cells
are geometry-only illustrations. The example JSON records aggregate design facts;
it contains no measured Ct/OD, source workbook, protocol, manuscript or legend.

Run from the repository root, with an output directory outside the checkout:

```bash
python scripts/authoring/main_figures_20261025/figure1/build_figure1.py --facts scripts/authoring/main_figures_20261025/figure1/figure1_facts.example.json --output-root ../private-figure-revision/figures
python scripts/authoring/main_figures_20261025/figure1/audit_figure1_exports.py --figure-root ../private-figure-revision/figures --audit-root ../private-figure-revision/qa/figure1
python scripts/authoring/main_figures_20261025/figure1/draw_layout_sketches.py --output-root ../private-figure-revision/layout-sketches
```

Dependencies: NumPy, matplotlib, Pillow and PyMuPDF for the export audit. Font
resolution and byte hashes are recorded at runtime. Exports can vary with library
and font versions; the original drawing functions are retained. The source
schematic reflects candidate associations and complementary evidence, not a
causal chain. PA prediction and downstream PRS quantification remain separate.
The optional experimental input provenance remains
`INTERNAL_REVIEW_SOURCE_UNRESOLVED`.

Design grammar was informed by the viewed official Figure 1 images of
[scDrugPrio](https://doi.org/10.1186/s13073-024-01314-7) and
[LIANA+](https://doi.org/10.1038/s41556-024-01469-w). All drawing primitives are
original; no publisher figure pixels, vector objects or reference-image files
are distributed. No full manuscript legend is bundled. Automated layout checks
leave actual visual review pending.
