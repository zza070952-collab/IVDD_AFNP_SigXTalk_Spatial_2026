# Figure-reference and private DOCX editing code

The mapping is simultaneous and semantic: old Figure 1B–E becomes 2A–D, old
Figures 2–4 become 3–5, and old Figure 5B–G becomes 6A–F. Deleted schematic
references map to the new overview or its experimental module. Supplementary
S references are outside the mapping grammar. Numerical-token guards exclude
explicit figure and bibliography labels, then compare the remaining tokens.

`figure_references.py` and its unit tests require only Python's standard library.
The remaining tools need `python-docx`, `lxml` and Pillow. Tests use only small
synthetic documents and coloured rectangles. The DOCX integration tests skip
explicitly when those optional packages are absent.

## Private source inputs

No manuscript, legend, original experimental data, complete study-specific
revision plan, source-document hash or exact research-result guard is included.
`plan.example.json` contains placeholders and illustrative synthetic indices;
it is intentionally not an executable study plan. The schema documents every
field. The CLI validates critical fields and source/hash/structural assertions;
it does not require or automatically run a JSON Schema validator.

Keep the source DOCX, author-approved overview, complete legend, exact paragraph
anchors, protected phrases, replacement text, reference exceptions, figure PNG
manifest and filled plan in a private directory outside this publication
checkout. Inspect the source with `inspect_docx.py` before setting its zero-based
paragraph indices. Ranges are half-open. `old_images_start` is valid only when
all later original paragraphs belong to the replaceable figure appendix.
The source retains its section break immediately before that appendix.

```bash
python inspect_docx.py /private/input/source.docx \
  --output /private/audit/source_inventory.json

python revise_manuscript.py \
  --source /private/input/source.docx \
  --plan /private/input/approved_plan.json \
  --figure-manifest /private/input/figures.json \
  --figure1-legend /private/input/figure1_legend.txt \
  --output-dir /private/new-edit
```

The private manifest maps string keys `"1"` through `"6"` to absolute paths of
the six final PNG files. Exactly one private legend insertion is required; its
panels must occur once each in A–E order. Generic headings may be overridden by
private `figure_headings` and `figure_alt_text` values to retain source-status
qualifiers. All study-specific prose and result checks come from the private
plan or legend file. The editor preserves unedited paragraph XML, selected
reference XML ranges, all protected tokens and supplementary reference spans;
it checks original image relationships and rejects specified obsolete image
payloads from the final package. The selected source is never written.

The editor writes a new DOCX and manuscript-bearing audit tables, deleted
legend spans, numeric-token reports and complete legends. `inspect_docx.py`
also writes full paragraph text. All tools require explicit output locations;
`output_guard.py` rejects paths within the nearest publication Git checkout.
Use a fresh private destination. Do not commit generated output or private plans.
A failed validation may leave a private draft DOCX but never a success receipt.
No success state certifies visual approval.

## Layout-only correction after visual inspection

```bash
python fix_legend_pagination.py \
  --source /private/new-edit/manuscript_renumbered.docx \
  --expected-sha256 SELECTED_INPUT_SHA256 \
  --output /private/layout/manuscript_renumbered.docx \
  --audit /private/layout/layout_audit.json
```

Optional `--legend-heading`, `--end-heading` and `--expected-groups` identify
generic legend boundaries. The output changes only `w:keepNext` and
`w:keepLines` inside those groups. Full text, images, all other package parts
and all other document XML are checked unchanged. Input and output are separate;
this command does not alter a previous edit receipt. Render the resulting DOCX
with a word processor and inspect every page before recording visual approval.

## Technical verification

```bash
python -B -m unittest discover -s . -p 'test_*.py' -v
```

The synthetic tests exercise the full edit, source-hash mismatch, scientific
numeric drift, missing protected phrases, obsolete-image removal, six-image
embedding, private inventory, output boundary protection and layout-only XML
comparison. No biological inference is run by these tools or tests. The fixed
mapping intentionally supports this six-main-figure transition; applying it to
an already-renumbered manuscript requires a new reviewed plan, not a second pass.
