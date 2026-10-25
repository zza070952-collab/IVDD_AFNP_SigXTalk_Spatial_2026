# External content authoring tools

These tools contain no manuscript paragraphs, study values, real figure legends, replacement maps, ROI data or machine paths. All content and layout decisions are supplied through private files at execution. Keep those inputs and generated documents outside Git, or under this folder's ignored `private_inputs/` and `.smoke/` directories. Never stage real layout/replacement JSON or generated artifacts.

Use Python 3.9 or later with the dependencies in `requirements.txt`. Native Word export is optional and requires Windows, Microsoft Word and `pywin32`. The tools do not download software, contact a service, train a model or redraw a scientific figure.

## PDF assembly

`assemble_pdf.py` uses ReportLab for text and PyMuPDF `show_pdf_page` for existing vector artwork. It keeps figure aspect ratio and validates effective source-figure font sizes. Pages are explicit; overflow causes an error instead of silent shrinking or omitted text. The standard geometry is A4, 15 mm side margins, 11 mm top clearance and 17 mm bottom clearance. Text styles default to black, with a declared minimum of 8 pt. The layout JSON can change these parameters.

```text
python scripts/authoring/common/assemble_pdf.py --layout <private-layout.json> --output <private-output.pdf> --receipt <private-build-audit.json>
```

A minimal external layout has this schema. The strings below are generic example filenames, not supplied content:

```json
{
  "metadata": {"title": "Example document", "author": ""},
  "minimum_font_pt": 8,
  "page": {"width_mm": 210, "height_mm": 297, "margin_mm": 15, "top_mm": 11, "bottom_mm": 17},
  "font_files": {"BodyFont": "fonts/body.ttf"},
  "styles": {"body": {"fontName": "BodyFont", "fontSize": 10.5, "leading": 13.2}},
  "pages": [
    {"blocks": [
      {"type": "markdown", "path": "section.md", "style": "body", "gap_pt": 8},
      {"type": "figure_pdf", "path": "artwork.pdf", "page_index": 0, "width_mm": 180, "gap_pt": 5},
      {"type": "markdown", "path": "caption.md", "style": "caption", "gap_pt": 0}
    ]}
  ],
  "footer": {"text_path": "footer.txt", "page_numbers": true, "font_pt": 8, "y_mm": 9}
}
```

Paths are relative to the layout JSON unless absolute. `font_files`, `styles`, `page`, `metadata` and `footer` are optional. Registered font aliases may be used by styles; provide all Unicode-capable and bold/italic font mappings required by private content. The built-in default is Helvetica. Default style names are `title`, `heading`, `body`, `caption` and `small`. Markdown support is deliberately limited to paragraphs, headings, inline bold/italic/code, and links rendered as label plus URL. Markdown tables and fenced code are rejected so they cannot silently lose layout. The assembler records source basenames/hashes and block geometry, without duplicating prose into its audit.

## PDF checks and external visual-review registration

```text
python scripts/authoring/common/pdf_qa.py audit-render --pdf <private-output.pdf> --output-dir <private-qa-directory> --minimum-font 8 --expected-pages <page-count> --review-template <private-pending-review.json>
python scripts/authoring/common/pdf_qa.py validate-review --pdf <private-output.pdf> --automated-report <private-qa-directory>/AUTOMATED_PDF_QA.json --manual-review <private-completed-review.json> --output <private-review-registration.json>
```

The first command checks page count, blank pages, extracted text sizes and page boundaries, and renders every page at 96 and 144 dpi by default. `--dpis` can explicitly select other positive DPI values. Automated success is named `AUTOMATED_CHECKS_PASS_PENDING_MANUAL_REVIEW`. It does **not** establish legible glyphs, absence of text overlap, scientific accuracy or an actual visual review. A successful build or Word export likewise leaves review pending.

The generated review template contains `PENDING`, empty reviewer/date/attestation fields, and one row for every page/DPI image. After actually viewing every listed image, an external reviewer supplies their name, UTC review time, page-specific status/notes, immutable render hashes, and the attestation `I inspected these rendered page images.` The schema is:

```json
{
  "pdf_sha256": "<current PDF SHA256>",
  "reviewer": "<external reviewer>",
  "reviewed_utc": "<UTC timestamp>",
  "attestation": "I inspected these rendered page images.",
  "pages": [
    {"page": 1, "dpi": 96, "render_sha256": "<render SHA256>", "status": "<PENDING or reviewer decision>", "notes": "<review notes>"}
  ]
}
```

`validate-review` accepts only complete external PASS records covering every recorded page/DPI once, bound to the current PDF and existing render hashes. Its success is explicitly `EXTERNAL_MANUAL_REVIEW_REGISTERED_PASS`: it validates a supplied record and does not perform or manufacture the review. No command automatically fills a manual PASS.

## Scoped reference replacements

`replace_references.py` reads exact literal replacements from private JSON. It supports `.md`, `.txt` and `.docx` and writes a new file. Matching is simultaneous, longest old string first, to prevent cascading reference renumbering. Each rule has a unique identifier and an exact expected match count. A mismatch aborts before output writing. `--dry-run` validates scopes/counts without writing a revised document.

```text
python scripts/authoring/common/replace_references.py --input <private-input.docx> --output <private-revised.docx> --plan <private-replacements.json> --receipt <private-change-audit.json>
```

The following is a generic schema example, not an actual replacement map:

```json
{
  "input_sha256": "<optional exact source SHA256>",
  "reference_label_patterns": ["\\bAppendix [A-Z]\\d*\\b"],
  "rules": [
    {"id": "label_key", "old": "<old reference label>", "new": "<new reference label>", "expected": 1,
     "scope": {"parts": ["word/document.xml"], "paragraph_indices": [0], "paragraph_regex": "<optional paragraph pattern>"}}
  ]
}
```

Scopes use zero-based `paragraph_indices`/`paragraph_regex` for DOCX, or `line_indices`/`line_regex` for Markdown/text. DOCX rules default to `word/document.xml`; other supported parts are headers, footers, footnotes and endnotes, explicitly named in `parts`. A match crossing DOCX runs is edited in text nodes without clearing paragraphs or deleting drawings/run properties. Live Word fields in a matched paragraph cause rejection; flatten or edit those fields in a separately authorized workflow. DOCX image/object member hashes must remain unchanged. Replacement audits record rule IDs, counts and paragraph/line text hashes, rather than copying manuscript prose.

Each changed paragraph/line must retain its ordered numeric tokens. Optional `reference_label_patterns` exclude only the explicit reference-label syntax from that comparison. Keep those patterns narrow and inspect them as part of the private edit authorization: broad patterns would weaken the guard. The guard detects numeric changes; it is not a substitute for checking that an authorized textual replacement preserves scientific meaning. Neither reference replacement nor a numeric guard constitutes visual QA.

## Optional DOCX formatting and native Word export

```text
python scripts/authoring/common/format_docx.py format --input <private-input.docx> --output <private-formatted.docx> --styles <private-styles.json> --receipt <private-format-audit.json>
python scripts/authoring/common/format_docx.py word-pdf --input <private-formatted.docx> --output <private-word-export.pdf> --receipt <private-export-audit.json>
```

External style JSON may contain `section_mm` (`page_width`, `page_height`, `top_margin`, `bottom_margin`, `left_margin`, `right_margin`) and a `styles` map using Word style names. Each style supports `font_name`, `font_pt`, `bold`, `space_before_pt`, `space_after_pt` and `line_spacing`. Formatting normalizes title/subtitle/heading color to black, removes underlining and paragraph-border residues, and verifies all visible text and embedded media bytes are unchanged. It does not rewrite content. Example style keys are `Normal` and `Title`; style values are local author choices, not embedded manuscript content.

Word export creates a separate invisible Word instance, opens the source read-only, exports PDF, closes it without saving and verifies the DOCX hash is unchanged. This optional Windows path is not part of the portable synthetic test. Run the PDF audit/render and external review workflow after any layout-affecting edit/export.

## Synthetic smoke

```text
python scripts/authoring/common/smoke_synthetic.py
```

The test creates only generic text/shapes under ignored `.smoke/`. It exercises PDF assembly, render/structure checks, rejection of pending manual review, Markdown replacement, split-run DOCX replacement, unchanged media bytes, rejection of changed numeric tokens, formatting text preservation, and overflow rejection. Its success status is `SYNTHETIC_AUTOMATED_CHECKS_PASS`; manual visual review stays pending. The test data, JSON plans, PDFs, DOCX files and render images are never public release inputs.
