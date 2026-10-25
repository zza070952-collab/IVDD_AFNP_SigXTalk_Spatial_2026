# Read-only figure and document QA

`inspect_final_main_figures.py --project PRIVATE_PROJECT --output PRIVATE_QA`
reads the six figures in `PRIVATE_PROJECT/figures/`. It records actual PDF fonts,
text bounds, editable SVG text, TIFF DPI, source hashes and previews/contact sheet.
It never issues a manual visual PASS.

On Windows with licensed Microsoft Word, `export_word_pdf.ps1 -ProjectRoot
PRIVATE_PROJECT -Docx PRIVATE_DOCX -Pdf PRIVATE_PDF -Receipt PRIVATE_JSON` starts
its own hidden read-only Word instance, exports PDF and checks the original hash.
Only the instance it creates is closed. Inputs and outputs must be under the
explicit private project root. The script has no server/account paths.

`render_pdf_pages.py --pdf PRIVATE_PDF --output-root PRIVATE_QA
[--poppler-path POPPLER_BIN] [--include-text]` produces 96/144-dpi page images,
a contact sheet and hash manifest. Additional dependencies: pdf2image, pypdf,
Pillow and an external Poppler installation. Full page text is opt-in; all these
outputs remain private. Rendering does not establish actual visual review.
No Word or Poppler runtime is bundled.
