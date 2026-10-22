# Frozen K5 provenance and reconstruction

The private original archive retains the27,264-gene by five-program regression spectra, the five-program by2,000-HVG consensus loadings and original fitting artifacts. Public code distributes only125 ranked gene memberships, the original byte-identical lock and provenance; no fitted model coefficients or weights are included. No cNMF fitting was performed to prepare this code bundle.

The original lock was written at `2026-10-12T09:55:25.163048+00:00`, using seed `20261012`, and has SHA256 `12905c6475cd4f3dd07c660efeee15e60d59f592c8f74c6cc60202f549110da2`. The all-gene score export has SHA256 `6566967a1e714875a3a023724ff55e17de18a2849027b5abfe629dd8bc68e985`. See the machine-readable manifests in `resources/frozen_k5/`.

## Actual historical preprocessing and fit

The original code read public 10x disease NP integer matrices from eight GSE230808 donors. It converted CTGF/CYR61/WISP1/WISP2/WISP3 to CCN2/CCN1/CCN4/CCN5/CCN6 and summed duplicate submitted symbols. Basic QC retained 300-9,000 detected canonical-symbol features, mitochondrial counts no greater than 25%, and positive total counts. It explicitly eliminated stored sparse zeros before counting features.

For the original lineage audit, the code normalised to 10,000 counts/cell, applied log1p, selected 2,000 Seurat-flavour HVGs, computed 30 ARPACK PCs and used fixed KMeans with 20 clusters, n_init=10 and seed20261012. Marker-dominance rules removed coherent non-disc identities; uncertain cells stayed in the audit but outside program fitting. This was a working marker-based NP-lineage definition, not definitive lineage tracing or an AF-to-NP transition analysis.

The original NP-lineage population contained 35,293 cells. Water-filling donor balancing and seeded sampling selected exactly 12,000 fitting cells, with 1,500 from each of eight donors. Raw sparse counts were supplied to official cNMF1.7.1. The exact call was `prepare(counts_fn=..., components=[5,7,9], n_iter=20, seed=20261012, num_highvar_genes=2000, max_NMF_iter=1000)`, followed by four one-thread workers, `combine(skip_missing_files=False)`, the official K-selection plot and `consensus(k=K, density_threshold=.5, local_neighborhood_size=.3, build_ref=False)` for each candidate K. The parameter record identifies the Frobenius loss/default normalisation. No batch correction or external-outcome optimisation was used in that original fit record.

The preserved execution log reports completion, and the preserved cNMF exit-code file contains0. Stability values were0.9990038963,0.7477682575 and0.8874034008 for K5,K7,K9, respectively. The main K5 was selected in the historical lock for its stability/parsimony; it is not reselected by this application pipeline. The lock's `external_outcomes_seen=false` describes that historical lock, not a claim that the current authors have not seen later results.

## Rank and version consistency

Every frozen Top25 member matches the all-gene score export, the original Top25 table, the later ranking export and the independent-cohort projection input. The original frozen ordering is authoritative. GEP5 ranks10/11, AC013287.1 and AL160408.1, share the exported score0.0064955824037216. A new platform's default unstable sort can exchange that tie; this does not change Top25 membership. The pipeline preserves the lock order instead of silently redefining it.

The later discovery identity reconstruction used seed20261016: one Scrublet run per sample (expected doublet rate0.06,30 PCs), then normalisation/log1p, batch-aware2,000 HVGs, scaling capped at10,30 PCs,15 neighbours, Leiden resolution0.5/two iterations and fixed marker rules. Its retained population is35,769 broad stromal cells, of which11,084 have high-confidence NP-like labels. This is a different identity/QC scope from the original35,293 NP-lineage cells and the12,000 fitting cells. The code keeps these populations distinct.

Historical versions were Python3.10.20, cNMF1.7.1, numpy1.26.4, scipy1.14.1, pandas2.3.3, anndata0.10.8, scanpy1.11.5 and scikit-learn1.7.2. Discovery UCell used R4.3.1, UCell2.4.0, Matrix1.6-4 and hdf5r1.3.8. Version records and own-code hashes are preserved. A full original package-wheel/OS lock is unavailable; installing nominally identical versions is not proof of bit-identical numerical reconstruction.

## Executable replay without an implicit refit

The portable adapters in `scripts/human/discovery_qc.py`, `discovery_identity.py` and `discovery_cnmf.py` retain the original own-code scientific operations. Their only path adaptation uses configured new workspaces. The discovery wrapper validates Scrublet completion rather than treating a failed-doublet fallback as success.

```bash
python scripts/00_download_human_public_data.py --config config/human.example.yaml --dataset discovery
python scripts/01_discovery_np_cnmf.py --config config/human.example.yaml --step raw-qc
python scripts/01_discovery_np_cnmf.py --config config/human.example.yaml --step identity
python scripts/01_discovery_np_cnmf.py --config config/human.example.yaml --step ucell
# Optional historical fitting replay; not needed for the current frozen analysis:
python scripts/01_discovery_np_cnmf.py --config config/human.example.yaml --step cnmf --execute-historical-cnmf
```

The optional fit recreates the candidate outputs in a new workspace; it never overwrites the distributed lock or selects a new positive program. Raw reconstruction has not been rerun merely to create this code release. Any newly reconstructed spectra must be compared with the retained lock before claiming identity.

Method citation: Kotliar et al., eLife2019, doi:10.7554/eLife.43803. Official cNMF is an external dependency, not code developed here. No third-party package implementation is vendored in this bundle.


## Private fitted-artifact audit

Optional `private_model_artifacts` entries in a private human config identify external fitted artifacts by role, path and expected SHA256. `scripts/audit_frozen_models.py --config PRIVATE --project-root PROJECT --require` verifies actual bytes without publishing the files. The normal membership projection has no full-spectra dependency. An explicit original raw/cNMF replay can generate new private spectra; compare versions, membership and retained original hashes before claiming numerical identity. The existing full-artifact hashes describe private historical inputs and are not a licence or exception to the exclusion of model weights.
