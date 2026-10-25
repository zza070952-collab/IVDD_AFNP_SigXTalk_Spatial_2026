# Frozen all-program spatial cache bridge

`project_frozen_spatial_programs.py` exposes the existing public
`scripts/spatial_programs.py` scoring implementation for GEP1–GEP5 by default.
This is the same sparse, strict-one-to-one, fixed-member, zero-retaining route
used by the earlier two-program projection. No ROI geometry, specimen data or
expression values are embedded here.

```bash
python scripts/authoring/v4/cache/project_frozen_spatial_programs.py \
  --prepared48 INPUT/prepared48 \
  --roi-metadata INPUT/roi_validation/DISC_AF_NP_ROI_METADATA.tsv.gz \
  --coverage INPUT/HUMAN_MOUSE_ORTHOLOG_COVERAGE.tsv \
  --output-dir NEW_PROJECT/02_analysis/spatial_allprograms
```

An old-layout input root may instead be given with `--input-work INPUT_WORK`.
The default program table is the distributed
`resources/frozen_k5/FROZEN_K5_PROGRAMS.tsv`. An external `--programs-file` is
accepted only if the complete 125 ranked memberships agree with that frozen
definition. `--programs` can explicitly select a subset; it never fits or ranks
new genes. `--mapping-file` is an alias for the per-capture `--coverage` file.
The global mapping-only snapshot does not contain native-feature coverage and
cannot substitute for it.

The prepared directory needs each capture's sparse `counts48_symbols_qc.npz`,
`BIN48_METADATA_QC.tsv.gz` and `MOUSE_SYMBOL_FEATURES.tsv`. Counts are bins by
native symbols. ROI metadata must retain the same spot order/coordinates and
NP_CORE/AF_CORE labels; coverage needs the recorded native matched symbol,
`MGI_ONE_TO_ONE` status and detection counts for all requested program/candidate
queries. The existing public preparation, ROI assignment and ortholog-coverage
scripts produce these inputs. They are required external runtime data, not
provided by a source-only clone.

The score is the unweighted mean of log1p(counts per 10,000) across fixed
available strict-one-to-one members, including measured zeros. Unavailable
mapping remains unavailable. The implementation records region, candidate,
member and all-bin source tables, including `NP_PROGRAM_BY_REGION.tsv` consumed
by V4 and compact builders at `02_analysis/spatial_allprograms`. It preserves
the existing interpretation that structures/bins are nested in captures with
unresolved animal identity; no genotype P value is added.

Use a new or empty output directory separate from input directories. The
wrapper hashes inputs before/after projection and writes a receipt with logical
input roles, output hashes and NumPy/pandas/SciPy versions, without recording
absolute machine paths. It does not write or replace original inputs, train a
model, alter programs/ROIs, or run the 100-UMI draw. Generated data and receipts
are runtime outputs and should stay outside this code repository.

This addition was syntax/CLI checked without replaying the large sparse data.
It supplies a missing cache-generation entry point, not evidence that a new
clean-machine GEO-to-all-V4-final run has completed. Frozen prepared inputs and
complete derived-data publication remain separate access/reproduction issues.
