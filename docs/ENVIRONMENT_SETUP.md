# Isolated environments

Use a new working directory and separate environments. Do not replace a shared server base environment. Configure executable paths in private module configurations; these paths are never committed.

1. The top dispatcher, human AUCell projection and spatial scripts use a scientific Python environment. Tu requires Python>=3.11 and the human module's pinned `envs/human_projection_requirements.txt`, including official decoupler2.2.0. NumPy/SciPy/pandas, matplotlib, Pillow, pyarrow and PyYAML are needed for the spatial route. Record the actual environment freeze. The completed spatial draw's recorded NumPy version matters for exact stochastic comparison.
2. Historical discovery/AF reconstruction uses the separately supplied `envs/historical_discovery_requirements.txt`. The original discovery Python was3.10.20. Historical discovery UCell uses R4.3.1, UCell2.4.0, Matrix1.6-4 and hdf5r1.3.8. An absent/incompatible UCell environment causes a failure, not an invented score substitute.
3. Official SigXTalk has its own R and older HGNN Python environment. Start from `envs/sigxtalk_r.bootstrap.yml` and read `envs/SIGXTALK_R_SETUP_NOTES.md`; the bootstrap file is a setup route, not a claim of a complete package lock. Use the exact final successful R session/installed-package receipt for an actual rerun. The current Seurat/SeuratObject pair must be compatible: an untested newer SeuratObject must not be silently swapped in. Configure `rscript`, `python`, `r_library_paths` and `official_source` in the B JSON. The older SigXTalk Python3.8 stack stays separate from the scientific Python environment.

```bash
# Example creation commands; install in isolated names chosen for this project.
python -m venv .venv-science
. .venv-science/bin/activate
python -m pip install -r envs/human_projection_requirements.txt
python -m pip install Pillow pyarrow PyYAML
# Use the historical discovery Python and compatible R separately.
# conda env create -f envs/sigxtalk_r.bootstrap.yml
```

The public downloader can fetch the exact pinned official SigXTalk archive and verify its critical Git blob identities. It stores that third-party source only under the private runtime root; the own-code repository does not vendor it. Set `downloads.include_official_demo=true` when the configured technical smoke consumes the official Zenodo RDS. No expected checksum is invented when only an actual downloaded hash is available.

The B wrapper pins the formal design to100epochs and seed20261017. Technical smoke length and its separate directory are declared in the B config before execution. An undertrained two-epoch smoke with no passing paths remains an unsuccessful closure attempt; the top dispatcher does not tune it, change its PA threshold or insert old predictions. A separately declared technical100-epoch run is still a subset check, not a complete published PBMC reproduction.

For a cached run, configure `human.tu_sparse_root` to actual Tu sparse caches and set the human/raw-selection module paths to genuine cell-level metadata and UCell/AUCell inputs. The candidate extractor requires actual expression matrices; summary scores cannot replace them. Configure `raw_selection.output_dir` to match B `input_dir`. The paired GSM and score-threshold TSVs are materialized from the verified public sample manifest and actual recomputed quartiles.

Every subprocess saves stdout/stderr, exit status and log SHA256. Do not call the workflow successful because packages imported or the technical tests passed. A real official model closure and the required measured outputs must exist.
