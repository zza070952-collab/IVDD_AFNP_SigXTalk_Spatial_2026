# New Figure 1 and main-figure renumbering

This optional own-code module reproduces the 25 October 2026 figure revision
from external frozen V4 inputs. It is separate from the scientific model
pipeline and does not train cNMF/SigXTalk/CellChat or redraw spatial depth samples.

| Final figure | Original content | Operation |
| --- | --- | --- |
| 1 | New study overview, modules A–E | Original editable vector drawing |
| 2 | Old Figure 1B–E | Replot as A–D; omit old scheme A |
| 3 | Old Figure 2 | Copy all four original formats byte-for-byte |
| 4 | Old Figure 3 | Copy all four original formats byte-for-byte |
| 5 | Old Figure 4 | Copy all four original formats byte-for-byte |
| 6 | Old Figure 5B–G | Replot as A–F; omit old scheme A |

The [Figure 1 guide](figure1/README.md) documents the schematic and its
non-data icons. The [Figure 2–6 guide](figures2_6/README.md) lists external
panel tables, read-only input checks and the independent legacy drawing oracle.
The [reference tools](references/README.md) apply simultaneous semantic figure
mapping; all actual manuscript prose, legends, paragraph anchors and source
documents are supplied in private files. The [QA tools](qa/README.md) provide
objective checks and render previews without declaring human visual approval.

Install the optional dependencies in a separate environment:

```bash
python -m pip install -r envs/authoring_requirements.txt
```

The optional Poppler page renderer additionally needs `pdf2image`, `pypdf` and
a system Poppler installation; native Word export needs Windows and licensed
Microsoft Word. No runtime or font binaries are redistributed.

Use explicit output locations outside this source checkout. No original
workbook/protocol, Ct/OD measurements, derived experimental tables, manuscripts,
full legends, page text, figures, review ZIPs or publisher reference images are
included. Runtime receipts may contain manuscript text or private input hashes;
keep them with the private project. The source publication gate reads actual
Git blobs and supplements manual review of embedded strings.

Expression support, predicted communication and anatomical spatial context stay
distinct. The experimental status remains
`INTERNAL_REVIEW_SOURCE_UNRESOLVED`; exporting or publishing code does not
authenticate laboratory sources. Automated layout PASS does not replace actual
rendered-page review. Portability and local component checks do not establish a
clean-machine reconstruction from all raw inputs.
