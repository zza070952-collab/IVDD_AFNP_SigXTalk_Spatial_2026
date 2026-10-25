# Private experimental-input analysis

This optional module publishes analysis code only. It contains no Ct/OD measurements, instrument files, sequences copied from a supplied protocol, figures, manuscript or result payload. Provide the workbook and protocol separately in a private workspace.

```bash
python scripts/authoring/v4/experiments/extract_inputs.py --workbook INPUT.xlsx --protocol PROTOCOL.docx --project ./analysis/experiment-review
python scripts/authoring/v4/experiments/analyze.py --project ./analysis/experiment-review
python scripts/authoring/v4/experiments/plot.py --project ./analysis/experiment-review
```

The original flat-table contract has worksheets `siRNA_qPCR_Raw`, `Program_qPCR_Raw` and `CCK8_Raw`. Source-sheet/row identifiers and undetermined Ct states are preserved. Input files are read-only, checked before/after and not copied into the repository. Oligo S/AS or F/R pairs must occur within the same original paragraph; the legacy adjacent-paragraph pairing error is not used. This fixes transcription, not assay validity.

The predeclared six-group, four-batch, technical-triplicate design is a method contract, not fabricated observations. The analyzer averages technical Ct, uses the two reference genes, estimates the two primary PLOD2 differences in differences with two-sided t intervals/tests and Holm correction, and retains the fixed secondary BH family. CCK-8 subtraction matches batch/time/coating blanks; it neither clips nonpositive corrected observations nor establishes equivalent viability. The figures derive statistics from runtime tables; no measured result is embedded as a fixture or shipped output.

Outputs keep `INTERNAL_REVIEW_SOURCE_UNRESOLVED`. Publishing the calculation code cannot certify laboratory provenance, instrument QC, protein-level effects or a causal mechanism. Original data and execution/QC records must be supplied and reconciled independently. Generated artifacts belong outside the source tree or under the ignored `analysis/` directory. Real-input validation logs are private and are not CI fixtures.
