The spatial code exposes the complete preparation and scoring chain rather
than only a plotting facade. All configurable paths in the example YAML are
relative to a chosen project root or to the YAML itself; no personal machine or
server root is embedded in the public scripts. Runtime count caches, source
images, absolute-path hash logs and all-bin expression tables are excluded from
the public-code allowlist.

Public input and immutable geometry

- Public GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287499
- Inspected archive: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE287nnn/GSE287499/suppl/GSE287499_RAW.tar
  (41,461,760 bytes; SHA256
  `cc405fd3acd18abec4e04199ea043a35de06eec3142181045dd9c42b3cddc652`).
- Captures: GSM8746576_CTR and GSM8746577_cKO. The archive provides sparse
  counts, feature/barcode files, native 16-um positions, scale factors and real
  H&E images. No external raw sequencing reconstruction is required by this
  limited prepared-count analysis.
- Non-overlapping 48-um groups are defined by global floor(array_row/3) and
  floor(array_col/3). A retained group requires all nine native positions, all
  nine matrix barcodes and all nine eligible in-tissue flags. Counts are summed
  sparsely, repeated symbols are merged, centroids use the native coordinates,
  and QC requires at least 20 counts and 20 detected symbols.
- The cKO crop has 1,896 original in-tissue positions outside the true H&E image.
  Original coordinates and the original tissue mask are retained in
  ORIGINAL_POSITIONS; a separate eligibility mask excludes those positions from
  aggregation. No coordinate clipping, translation or axis swap is performed.
- Frozen JSON provides eight NP and sixteen bilateral AF polygons, source image
  hashes, capture-specific pixel scale and artifact exclusions. Only centers
  at least the 48-um bin half diagonal from ROI edges contribute. No CEP is
  assigned. Independent structural-marker checks follow fixed geometric
  assignment; unsupported regions remain UNRESOLVED instead of being redrawn.
- Limited targeted human/mouse mapping is retained with official MGI/Ensembl
  provenance and SHA256. MGI one-to-one protein-coding relationships define
  eligibility; Ensembl stable mouse IDs help identify released old symbols.
  Exact MGI symbols are a documented fallback. There is no title-case orthology
  inference. The frozen targeted table prevents mutable current sources from
  silently changing the completed analysis.

Raw mode (new private working cache)

```bash
python scripts/07_prepare_spatial.py --config configs/spatial.example.yaml \
  --project-root "$PROJECT" --from-raw --download-public
python scripts/08_spatial_depth_check.py --config configs/spatial.example.yaml \
  --project-root "$PROJECT" --from-raw
python scripts/plot_figure4.py \
  --spatial-input "$PROJECT/data_cache/spatial_stage" \
  --depth-input "$PROJECT/03_spatial_gse287499" \
  --roi-json configs/frozen_spatial_rois.json --out "$PROJECT/05_figures"
```

The first command uses a size/hash-verified archive already present in cache,
or downloads it from public NCBI when explicitly enabled. It executes sparse
preparation, frozen ROI assignment, ortholog coverage and fixed program/candidate
scoring, saving stage exit codes and logs. This raw entry point intentionally
does not claim to reconstruct the upstream human cNMF experiment; it consumes
the separately frozen K5 program table. Baseline spatial cache and new depth
output must be separate directories.

Cached mode (previously prepared, read-only baseline)

```bash
python scripts/07_prepare_spatial.py --config configs/spatial.example.yaml \
  --project-root "$PROJECT" --from-cached-inputs --cached-spatial "$BASELINE"
python scripts/08_spatial_depth_check.py --config configs/spatial.example.yaml \
  --project-root "$PROJECT" --cached-spatial "$BASELINE"
python scripts/plot_figure4.py --spatial-input "$BASELINE" \
  --depth-input "$PROJECT/03_spatial_gse287499" \
  --roi-json configs/frozen_spatial_rois.json --out "$PROJECT/05_figures"
```

Cached mode explicitly records that archive fetch, sparse preparation, ROI
assignment, mapping and scoring were skipped after validating baseline files.
The depth module nevertheless reconstructs original scores from the real sparse
counts, checks coordinates against frozen metadata and verifies protected input
hashes before and after the one fixed-count analysis. Figure4 validates source
image hashes, the fixed ROI hash, and original versus depth-table disc contrasts.

Required spatial Python packages are numpy, scipy, pandas, matplotlib, Pillow,
pyarrow (native parquet positions), and PyYAML. Use the recorded complete-run
environment for exact stochastic reproduction: NumPy sampling algorithms can
change across versions, even with the same seed. SciPy sparse formats avoid
dense whole-expression matrices. A per-gene vector is materialized only when
needed for an explicitly requested score/table. Reports record package versions.

Verified execution boundary

The actual V3 fixed-count analysis ran on the frozen real sparse matrices and
the amended Figure4 was regenerated from the resulting tables. Original
preparation/ROI/mapping/scoring ran in V2 and are being exposed here as full
parameterized scripts. The newly assembled all-in-one raw entry point has
static syntax checks and component technical tests; unless a subsequent run
log is supplied, it has not been rerun end-to-end from the public archive in
this release. The complete project must not be described as independently
100% from-raw reproduced on that basis. No technical test fixture is used as
research evidence.

Depth sensitivity is intentionally limited to one prespecified target, one
fixed seed and existing fixed cores. Eligible AF retention is low, and the
result substantially weakens the original NP contrast. It establishes neither
absence of the biological program nor an unconfounded genotype effect. Full
candidate detectability tables remain original and unadjusted.
