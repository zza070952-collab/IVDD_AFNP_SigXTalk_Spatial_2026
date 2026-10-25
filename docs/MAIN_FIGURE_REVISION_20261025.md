# Main-figure revision code, 25 October 2026

The optional [main-figure module](../scripts/authoring/main_figures_20261025/README.md)
adds the original study-overview Figure 1, relays old Figure 1B–E to new Figure
2A–D and old Figure 5B–G to new Figure 6A–F, and shifts old Figures 2–4 to new
Figures 3–5 without changing their artwork bytes. Its source-based replots
remove the two former schematic panels and redistribute space to data panels.

This version's figure numbers apply to the authoring revision. Earlier pipeline
scripts and V4 authoring scripts retain their historical numbering. They are
not silently rewritten or rerun. Supplementary figures and source data are not
renumbered by this module.

External inputs are necessary: frozen V4 human panel tables, private derived
experimental summary tables for the optional experimental replot, original
retained artwork, and private manuscript/plan/legend inputs for DOCX revision.
The repository supplies code and a small aggregate design configuration only.
The independent comparison utility can execute selected drawing statements
from explicitly trusted local legacy code. It never runs the legacy full main
or scientific export workflow; inspect that local source before opting in.

Technical checks cover semantic figure references and protected non-figure
numbers, read-only input hashing, unchanged retained artwork, data-coordinate
artist comparisons, editable SVG/PDF text, minimum font size and page bounds.
Private local component checks use existing authorized frozen inputs; their
research content and full logs remain outside this repository. Actual visual
review is a separate recorded action, never inferred from file existence.

No raw experiment, protocol/oligo content, manuscript prose or full legend is
published. The DOCX tool's plan is private because moving prose into public
JSON would still publish the manuscript. Private outputs must not be committed.
The experimental source status remains unresolved, and this code release does
not constitute a new biological analysis, source authentication or independent
raw-to-final reproduction.
