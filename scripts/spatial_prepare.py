#!/usr/bin/env python3
"""Prepare real GSE287499 16-um counts at 48 um without anatomical inference.

Scientific preparation only: source H&E is copied unchanged, aggregation uses
the original physical grid, and anatomical labels remain UNASSIGNED. This file
contains no program/candidate expression, ROI selection, or genotype tests.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import importlib.metadata as md
import io
import json
import os
from pathlib import Path
import platform
import re
import sys
import tarfile
import traceback

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np
import pandas as pd
from scipy import io as sio, sparse

CAPTURES = [("GSM8746576", "CTR"), ("GSM8746577", "cKO")]
FILES = {
    "matrix": "matrix.mtx.gz", "barcodes": "barcodes.tsv.gz",
    "features": "features.tsv.gz", "positions": "tissue_positions.parquet.gz",
    "scalefactors": "scalefactors_json.json.gz",
    "hires": "tissue_hires_image.png.gz", "lowres": "tissue_lowres_image.png.gz",
    "detected": "detected_tissue_image.jpg.gz", "fiducials": "aligned_fiducials.jpg.gz",
}


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda: f.read(4 * 1024**2), b""):
            h.update(data)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf8")


class Source:
    """Read named members without extractall, symlinks, or modifying the input."""
    def __init__(self, raw_tar=None, input_dir=None):
        self.tar = tarfile.open(raw_tar, "r:*") if raw_tar else None
        self.root = Path(input_dir).resolve() if input_dir else None
        self.index = {}
        entries = self.tar.getmembers() if self.tar else self.root.rglob("*")
        for entry in entries:
            if self.tar:
                if not entry.isfile():
                    continue
                name = Path(entry.name).name
            else:
                if not entry.is_file():
                    continue
                name = entry.name
            self.index.setdefault(name, []).append(entry)
        self.manifest = []

    def get(self, basename):
        entries = self.index.get(basename, [])
        if len(entries) != 1:
            raise ValueError(f"Expected exactly one source {basename}; found {len(entries)}")
        entry = entries[0]
        if self.tar:
            stream = self.tar.extractfile(entry)
            if stream is None:
                raise ValueError(f"Unreadable TAR member: {entry.name}")
            data = stream.read()
            source = entry.name
        else:
            data = entry.read_bytes()
            source = str(entry)
        self.manifest.append({"filename": basename, "source": source, "compressed_bytes": len(data)})
        return gzip.decompress(data) if basename.endswith(".gz") else data

    def close(self):
        if self.tar:
            self.tar.close()


def collapse_symbols(x, features):
    codes, symbols = pd.factorize(features.mouse_gene_symbol, sort=False)
    projection = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.int64), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(symbols)),
    )
    xs = (x @ projection).tocsr()
    assert np.array_equal(np.asarray(x.sum(axis=1)).ravel(), np.asarray(xs.sum(axis=1)).ravel())
    mapping = features.copy()
    mapping["symbol_aggregate_index"] = codes
    var = mapping.groupby("symbol_aggregate_index", sort=True).agg(
        mouse_gene_symbol=("mouse_gene_symbol", "first"),
        mouse_ensembl_ids=("mouse_ensembl_id", lambda z: ";".join(z)),
        n_original_features=("mouse_ensembl_id", "size"),
    ).reset_index()
    assert list(var.mouse_gene_symbol) == list(symbols)
    return xs, var, mapping


def aggregate48(x, barcodes, positions, gsm):
    """Sum non-overlapping global 3x3 groups; no missing-bin imputation."""
    pos = positions.copy().set_index("barcode", drop=False)
    assert pos.index.is_unique and barcodes.is_unique
    if not barcodes.isin(pos.index).all():
        raise ValueError("Matrix barcodes missing from physical positions")
    if pos.duplicated(["array_row", "array_col"]).any():
        raise ValueError("Duplicate original grid row/column")
    pos["agg_grid_row"] = pos.array_row.astype(np.int64) // 3
    pos["agg_grid_col"] = pos.array_col.astype(np.int64) // 3
    pos["grid_id"] = pos.agg_grid_row.astype(str) + "_" + pos.agg_grid_col.astype(str)
    pos["has_matrix_barcode"] = pos.index.isin(barcodes)
    groups = pos.groupby("grid_id", sort=False).agg(
        agg_grid_row=("agg_grid_row", "first"), agg_grid_col=("agg_grid_col", "first"),
        n_native_positions=("barcode", "size"), n_matrix_barcodes=("has_matrix_barcode", "sum"),
        n_in_tissue=("in_tissue", "sum"), x_fullres=("pxl_col_in_fullres", "mean"),
        y_fullres=("pxl_row_in_fullres", "mean"),
    )
    groups["complete_in_tissue"] = (groups.n_native_positions == 9) & (groups.n_matrix_barcodes == 9) & (groups.n_in_tissue == 9)
    chosen = groups.loc[groups.complete_in_tissue].copy()
    if len(chosen) == 0:
        raise ValueError("No complete nine-bin in-tissue groups")
    lookup = dict(zip(chosen.index, range(len(chosen))))
    source = pos.loc[barcodes]
    take = source.grid_id.isin(lookup)
    columns = np.flatnonzero(take.to_numpy())
    rows = source.loc[take, "grid_id"].map(lookup).to_numpy(dtype=np.int64)
    aggregation = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.int64), (rows, columns)), shape=(len(chosen), len(barcodes))
    )
    assert np.all(np.asarray(aggregation.sum(axis=1)).ravel() == 9)
    y = (aggregation @ x).tocsr()
    assert int(y.sum()) == int(x[columns].sum())
    chosen["spot_id"] = [f"{gsm}|agg48_r{rr:03d}_c{cc:03d}" for rr, cc in zip(chosen.agg_grid_row, chosen.agg_grid_col)]
    by_grid = dict(zip(chosen.index, chosen.spot_id))
    mapping = source[["barcode", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres", "in_tissue", "grid_id"]].copy()
    mapping["analysis_spot_id"] = mapping.grid_id.map(by_grid).fillna("")
    mapping["retained_complete_48um_group"] = mapping.analysis_spot_id != ""
    return y, chosen.reset_index(), groups.reset_index(), mapping.reset_index(drop=True)


def overview(folder, obs, scale, gsm, genotype):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "axes.titlesize": 9, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False})
    im = np.asarray(Image.open(folder / "source_images" / f"{gsm}_{genotype}_tissue_hires_image.png"))
    s = float(scale["tissue_hires_scalef"])
    extent = (0, im.shape[1] / s, im.shape[0] / s, 0)
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.4), layout="constrained")
    for ax in axes:
        ax.imshow(im, origin="upper", extent=extent)
        ax.set_aspect("equal")
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], 0)
        ax.set_xlabel("Original full-resolution x (pixels)")
        ax.set_ylabel("Original full-resolution y (pixels)")
    axes[0].set_title(f"{gsm} {genotype}: published H&E")
    scatter = axes[1].scatter(obs.x_fullres, obs.y_fullres, c=np.log1p(obs.n_counts),
                              s=3, marker="s", cmap="viridis", alpha=.75, linewidths=0,
                              rasterized=True)
    axes[1].set_title("Complete 48-µm groups: RNA coverage")
    fig.colorbar(scatter, ax=axes[1], shrink=.65, label="log(1 + UMI count)")
    fig.suptitle("One capture per genotype; anatomical compartments not yet assigned", fontsize=10)
    for ext in ["png", "pdf"]:
        fig.savefig(folder / f"{gsm}_{genotype}_HE_COVERAGE_OVERVIEW.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    obs[["spot_id", "x_fullres", "y_fullres", "n_counts", "n_genes", "pass_basic_qc"]].to_csv(
        folder / "HE_COVERAGE_OVERVIEW_SOURCE.tsv.gz", sep="\t", index=False)


def prepare_capture(source, out, gsm, genotype, min_counts, min_genes):
    folder = out / f"{gsm}_{genotype}"
    (folder / "source_images").mkdir(parents=True, exist_ok=True)
    prefix = f"{gsm}_{genotype}_"
    read = lambda key: source.get(prefix + FILES[key])
    print(f"{utc()} READ {gsm} {genotype}", flush=True)
    with io.BytesIO(read("matrix")) as f:
        original = sio.mmread(f)
    if not sparse.issparse(original):
        raise ValueError("Expected a sparse MatrixMarket coordinate matrix")
    x = original.T.tocsr()
    x.sum_duplicates()
    if not (np.isfinite(x.data).all() and (x.data >= 0).all() and (x.data == np.floor(x.data)).all()):
        raise ValueError("Original matrix is not finite nonnegative integer counts")
    x = x.astype(np.int64)
    barcodes = pd.Index(pd.read_csv(io.BytesIO(read("barcodes")), header=None, keep_default_na=False)[0].astype(str))
    features = pd.read_csv(io.BytesIO(read("features")), header=None, sep="\t", keep_default_na=False)
    if features.shape[1] != 3:
        raise ValueError(f"Expected three feature columns, got {features.shape[1]}")
    features.columns = ["mouse_ensembl_id", "mouse_gene_symbol", "feature_type"]
    assert len(barcodes) == x.shape[0] and len(features) == x.shape[1]
    assert features.mouse_gene_symbol.ne("").all() and features.mouse_ensembl_id.ne("").all()
    if not features.feature_type.eq("Gene Expression").all():
        raise ValueError("Non-gene-expression feature types require explicit review")
    pos = pd.read_parquet(io.BytesIO(read("positions")))
    expected = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    if not set(expected).issubset(pos.columns):
        raise ValueError(f"Missing position columns: {set(expected) - set(pos.columns)}")
    assert pos[expected].notna().all().all()
    assert pos.in_tissue.isin([0, 1]).all()
    scale = json.loads(read("scalefactors"))
    if float(scale.get("bin_size_um", -1)) != 16:
        raise ValueError(f"Expected original 16-um bins, got scale factors {scale}")
    mpp = float(scale["microns_per_pixel"])
    assert np.isfinite(mpp) and mpp > 0
    parts = pos.barcode.astype(str).str.extract(r"^s_016um_(\d+)_(\d+)-\d+$")
    if parts.isna().any().any():
        raise ValueError("Unexpected native barcode format")
    if not (np.array_equal(parts[0].astype(int), pos.array_row.astype(int)) and
            np.array_equal(parts[1].astype(int), pos.array_col.astype(int))):
        raise ValueError("Barcode row/column disagree with official positions; no automatic axis swap")
    coordinates = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    if not np.isfinite(coordinates).all():
        raise ValueError("Nonfinite physical coordinates")
    image_shapes = {}
    from PIL import Image
    for key in ["hires", "lowres", "detected", "fiducials"]:
        data = read(key)
        name = (prefix + FILES[key])[:-3]
        p = folder / "source_images" / name
        p.write_bytes(data)
        with Image.open(p) as im:
            im.verify()
        with Image.open(p) as im:
            image_shapes[key] = list(im.size)
    hires_scale = float(scale["tissue_hires_scalef"])
    assert hires_scale > 0
    x_max, y_max = np.asarray(image_shapes["hires"], dtype=float) / hires_scale
    native_in_tissue = pos.loc[pos.in_tissue == 1]
    # The released cKO grid extends above the cropped H&E (1,896 original
    # in-tissue positions have negative y). Keep the supplied coordinates and
    # original mask unchanged in ORIGINAL_POSITIONS; exclude those positions
    # only from aggregation eligibility. Never shift, swap or clip coordinates.
    in_image = pos.pxl_col_in_fullres.between(0, x_max) & pos.pxl_row_in_fullres.between(0, y_max)
    pos["inside_supplied_HE_bounds"] = in_image
    n_outside_in_tissue = int(((pos.in_tissue == 1) & ~in_image).sum())
    eligible_positions = pos.copy()
    eligible_positions.loc[~in_image, "in_tissue"] = 0
    xs, var, feature_map = collapse_symbols(x, features)
    y, obs, geometry, mapping = aggregate48(xs, barcodes, eligible_positions, gsm)
    obs["n_counts"] = np.asarray(y.sum(axis=1)).ravel()
    obs["n_genes"] = np.asarray((y > 0).sum(axis=1)).ravel()
    obs["pass_basic_qc"] = (obs.n_counts >= min_counts) & (obs.n_genes >= min_genes)
    obs["source_GSM"] = gsm
    obs["genotype"] = genotype
    obs["capture_id"] = f"{gsm}_{genotype}"
    obs["disc_id"] = "UNASSIGNED"
    obs["anatomical_region"] = "UNASSIGNED"
    obs["individual_animal_id"] = "NOT_IDENTIFIABLE_FROM_GEO"
    obs["native_bin_um"] = 16
    obs["analysis_bin_um"] = 48
    obs["x_um"] = obs.x_fullres * mpp
    obs["y_um"] = obs.y_fullres * mpp
    sparse.save_npz(folder / "counts48_symbols_all_complete.npz", y)
    sparse.save_npz(folder / "counts48_symbols_qc.npz", y[obs.pass_basic_qc.to_numpy()].tocsr())
    obs.to_csv(folder / "BIN48_METADATA_ALL_COMPLETE.tsv.gz", sep="\t", index=False)
    obs.loc[obs.pass_basic_qc].to_csv(folder / "BIN48_METADATA_QC.tsv.gz", sep="\t", index=False)
    var.to_csv(folder / "MOUSE_SYMBOL_FEATURES.tsv", sep="\t", index=False)
    feature_map.to_csv(folder / "ORIGINAL_FEATURE_TO_SYMBOL.tsv.gz", sep="\t", index=False)
    geometry.to_csv(folder / "ALL_GRID48_GEOMETRY.tsv.gz", sep="\t", index=False)
    mapping.to_csv(folder / "NATIVE_TO_48UM_MAPPING.tsv.gz", sep="\t", index=False)
    pos.to_csv(folder / "ORIGINAL_POSITIONS.tsv.gz", sep="\t", index=False)
    dump(folder / "ORIGINAL_SCALEFACTORS.json", scale)
    overview(folder, obs, scale, gsm, genotype)
    qc = {
        "source_GSM": gsm, "genotype": genotype, "capture_count": 1, "source_disc_count": 4,
        "individual_animal_identity": "NOT_IDENTIFIABLE_FROM_GEO", "native_bin_um": 16,
        "analysis_bin_um": 48, "native_matrix_n_barcodes": int(x.shape[0]),
        "original_feature_rows": int(x.shape[1]), "unique_mouse_symbols": int(xs.shape[1]),
        "native_nnz": int(x.nnz), "native_positions": len(pos),
        "native_in_tissue_positions": int(pos.in_tissue.sum()),
        "original_in_tissue_positions_outside_HE_excluded": n_outside_in_tissue,
        "coordinate_policy": "Original coordinates unchanged; full 3x3 groups require all nine original in-tissue positions inside supplied H&E bounds",
        "grid48_groups_total": len(geometry), "complete_nine_bin_in_tissue_groups": len(obs),
        "qc_min_UMI": min_counts, "qc_min_detected_symbols": min_genes,
        "qc_retained_48um_bins": int(obs.pass_basic_qc.sum()),
        "median_UMI_before_qc": float(obs.n_counts.median()),
        "median_genes_before_qc": float(obs.n_genes.median()),
        "zero_count_48um_bins": int(obs.n_counts.eq(0).sum()),
        "total_UMI_complete_groups": int(y.sum()), "sparse_48um_nnz": int(y.nnz),
        "microns_per_pixel": mpp, "hires_image_width_pixels": image_shapes["hires"][0],
        "hires_image_height_pixels": image_shapes["hires"][1],
        "anatomical_ROI_status": "UNASSIGNED_PENDING_HE_REVIEW", "program_scores_computed": False,
        "genotype_inference": "NOT_PERFORMED_CAPTURE_CONFOUNDED",
        "finished_utc": utc(),
    }
    dump(folder / "CAPTURE_PREPARATION_COMPLETE.json", qc)
    print(f"{utc()} PREPARED {gsm}: {len(obs)} geometry bins, {qc['qc_retained_48um_bins']} QC bins", flush=True)
    return qc


def self_test():
    """A tiny mathematical test, never a scientific input or plot."""
    rr, cc = np.indices((6, 6))
    p = pd.DataFrame({"barcode": [f"s_016um_{r:05d}_{c:05d}-1" for r, c in zip(rr.ravel(), cc.ravel())],
                      "array_row": rr.ravel(), "array_col": cc.ravel(), "in_tissue": 1,
                      "pxl_row_in_fullres": rr.ravel() * 4.0, "pxl_col_in_fullres": cc.ravel() * 4.0})
    a = sparse.csr_matrix(np.arange(72, dtype=np.int64).reshape(36, 2))
    y, obs, geometry, mapping = aggregate48(a, pd.Index(p.barcode), p, "TEST")
    assert y.shape == (4, 2) and int(y.sum()) == int(a.sum())
    assert np.allclose(obs.x_fullres, [4, 16, 4, 16]) and np.allclose(obs.y_fullres, [4, 4, 16, 16])
    p.loc[0, "in_tissue"] = 0
    y2, obs2, _, mapping2 = aggregate48(a, pd.Index(p.barcode), p, "TEST")
    assert len(obs2) == 3 and mapping2.retained_complete_48um_group.sum() == 27
    assert int(y2.sum()) == int(a.sum()) - int(y[0].sum())
    print("SELF_TEST_PASS: sparse count conservation, full nine-bin mask, global nonoverlap, physical centroids")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--raw-tar", type=Path)
    group.add_argument("--input-dir", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--min-counts", type=int, default=20)
    ap.add_argument("--min-genes", type=int, default=20)
    ap.add_argument("--capture", choices=[g for g, _ in CAPTURES], action="append")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.out is None or not (args.raw_tar or args.input_dir):
        ap.error("--out and either --raw-tar or --input-dir are required")
    out = args.out.resolve()
    if args.input_dir and (out == args.input_dir.resolve() or out.is_relative_to(args.input_dir.resolve())):
        raise ValueError("Output must not be inside the source input directory")
    out.mkdir(parents=True, exist_ok=True)
    source = Source(args.raw_tar, args.input_dir)
    status = {"status": "RUNNING", "started_utc": utc(), "command": sys.argv,
              "raw_tar": str(args.raw_tar.resolve()) if args.raw_tar else None,
              "raw_tar_sha256": sha256(args.raw_tar) if args.raw_tar else None,
              "raw_tar_bytes": args.raw_tar.stat().st_size if args.raw_tar else None,
              "script_sha256": sha256(__file__), "captures": [], "errors": [],
              "biological_association_results": False, "ROI_or_program_inputs_read": False}
    dump(out / "PREPARATION_STATUS.json", status)
    for gsm, genotype in CAPTURES:
        if args.capture and gsm not in args.capture:
            prior = out / f"{gsm}_{genotype}" / "CAPTURE_PREPARATION_COMPLETE.json"
            if prior.exists():
                old = json.loads(prior.read_text())
                old["reused_completed_preparation"] = True
                status["captures"].append(old)
            continue
        try:
            status["captures"].append(prepare_capture(source, out, gsm, genotype, args.min_counts, args.min_genes))
        except Exception as exc:
            status["errors"].append({"source_GSM": gsm, "error": str(exc), "traceback": traceback.format_exc()})
            print(traceback.format_exc(), flush=True)
        dump(out / "PREPARATION_STATUS.json", status)
    source.close()
    pd.DataFrame(source.manifest).to_csv(out / "READ_SOURCE_MEMBERS.tsv", sep="\t", index=False)
    pd.DataFrame(status["captures"]).to_csv(out / "SPATIAL_INPUT_QC.tsv", sep="\t", index=False)
    versions = {"python": platform.python_version()}
    for package in ["numpy", "scipy", "pandas", "pyarrow", "matplotlib", "Pillow"]:
        try:
            versions[package] = md.version(package)
        except md.PackageNotFoundError:
            versions[package] = "NOT_INSTALLED"
    dump(out / "PREPARATION_VERSIONS.json", versions)
    status["status"] = "DONE_INPUT_PREPARATION" if len(status["captures"]) == 2 and not status["errors"] else "PARTIAL_INPUT_PREPARATION"
    status["finished_utc"] = utc()
    dump(out / "PREPARATION_STATUS.json", status)
    print(json.dumps({k: status[k] for k in ["status", "errors", "finished_utc"]}, ensure_ascii=False), flush=True)
    return 0 if not status["errors"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
