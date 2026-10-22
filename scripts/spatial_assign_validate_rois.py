#!/usr/bin/env python3
"""Assign frozen H&E cores and describe independent structural-marker support.

The input polygon is never edited. Programs and candidate genes are not read.
Marker incompatibility leaves the final label UNRESOLVED, preserving the
original H&E label separately. This is provisional spatial localization, not
expert pathology, independent-capture inference or cell-type classification.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath
from matplotlib.patches import Circle, Polygon
from PIL import Image

NP_MARKERS = ["Krt19", "T", "Krt8", "Acan"]
AF_MARKERS = ["Chad", "Clec3a"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf8")


def checked_image(capture, path):
    if sha(path) != capture["source_image_sha256"]:
        raise ValueError(f"Original H&E hash mismatch for {capture['capture_id']}")
    with Image.open(path) as im:
        image = np.asarray(im).copy()
        if list(im.size) != capture["source_image_pixels"]:
            raise ValueError("Original H&E dimensions differ from frozen annotation")
    return image


def plot_frozen(capture, image, out):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(5.4, 8.5), layout="constrained")
    ax.imshow(image, origin="upper")
    palette = {"NP_CORE": "#286EAE", "AF_CORE": "#C67922"}
    for disc in capture["discs"]:
        for region in disc["regions"]:
            vertices = np.asarray(region["vertices"], dtype=float)
            ax.add_patch(Polygon(vertices, closed=True, fill=False, edgecolor=palette[region["region"]], linewidth=1.0))
        vertices = np.asarray(disc["regions"][0]["vertices"], dtype=float)
        ax.text(vertices[:, 0].mean(), vertices[:, 1].min() - 24, disc["disc_id"], ha="center", va="bottom",
                fontsize=8, bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none", "pad": 1})
    margin = np.sqrt(2) * 48 / (2 * capture["microns_per_pixel"])
    for circle in capture.get("exclusion_circles", []):
        ax.add_patch(Circle(circle["center"], circle["radius_pixels"] + margin, fill=False,
                            edgecolor="#B43C39", linestyle="--", linewidth=.8))
    x0, x1, y0, y1 = capture.get("overview_display_window_pixels", [0, image.shape[1], 0, image.shape[0]])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal")
    ax.set_xlabel("Original H&E x (pixels)")
    ax.set_ylabel("Original H&E y (pixels)")
    ax.set_title(f"{capture['capture_id']}: frozen provisional H&E cores\nBlue: NP core; ochre: AF core; dashed red: artifact exclusion", fontsize=9)
    for ext in ["png", "pdf"]:
        fig.savefig(out / f"{capture['capture_id']}_FROZEN_HE_CORE_OVERVIEW.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def edge_distance(points, vertices):
    minimum = np.full(len(points), np.inf)
    for first, last in zip(vertices, np.roll(vertices, -1, axis=0)):
        line = last - first
        denominator = float(line @ line)
        if denominator == 0:
            raise ValueError("Repeated adjacent polygon vertex")
        position = np.clip(((points - first) @ line) / denominator, 0, 1)
        nearest = first + position[:, None] * line
        minimum = np.minimum(minimum, np.sqrt(((points - nearest) ** 2).sum(axis=1)))
    return minimum


def assign_cores(obs, capture, margin):
    result = obs.copy()
    points = obs[["x_fullres", "y_fullres"]].to_numpy(dtype=float)
    n = len(obs)
    result["disc_id"] = "UNASSIGNED"
    result["histology_core_region"] = "UNASSIGNED"
    result["histology_core_side"] = "UNASSIGNED"
    result["assignment_status"] = "OUTSIDE_FROZEN_CORES"
    result["distance_to_core_edge_pixels"] = np.nan
    result["edge_margin_pixels"] = margin
    winner = np.full(n, -1, dtype=int)
    overlaps = np.zeros(n, dtype=bool)
    entries = []
    for disc in capture["discs"]:
        for region in disc["regions"]:
            if region["region"] not in ["NP_CORE", "AF_CORE"]:
                raise ValueError("Only frozen NP_CORE/AF_CORE are supported; no implicit CEP assignment")
            vertices = np.asarray(region["vertices"], dtype=float)
            if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 3 or not np.isfinite(vertices).all():
                raise ValueError("Invalid polygon vertices")
            inside = PolygonPath(np.vstack([vertices, vertices[0]]), closed=True).contains_points(points)
            distance = edge_distance(points, vertices)
            boundary = inside & (distance < margin)
            untouched_boundary = boundary & result.assignment_status.eq("OUTSIDE_FROZEN_CORES").to_numpy()
            result.loc[untouched_boundary, "disc_id"] = disc["disc_id"]
            result.loc[untouched_boundary, "assignment_status"] = "FROZEN_CORE_EDGE_BUFFER_UNASSIGNED"
            safe = inside & (distance >= margin)
            overlaps |= safe & (winner >= 0)
            index = len(entries)
            entries.append((disc["disc_id"], region["region"], region.get("side", "unspecified")))
            winner[safe] = index
            result.loc[safe, "distance_to_core_edge_pixels"] = distance[safe]
    for index, (disc_id, region, side) in enumerate(entries):
        take = winner == index
        result.loc[take, ["disc_id", "histology_core_region", "histology_core_side", "assignment_status"]] = [disc_id, region, side, "FROZEN_HE_CORE"]
    result.loc[overlaps, ["histology_core_region", "histology_core_side", "assignment_status"]] = ["UNRESOLVED", "UNRESOLVED", "OVERLAPPING_FROZEN_CORES_UNRESOLVED"]
    for circle in capture.get("exclusion_circles", []):
        excluded = np.sqrt(((points - np.asarray(circle["center"])) ** 2).sum(axis=1)) <= circle["radius_pixels"] + margin
        result.loc[excluded, ["disc_id", "histology_core_region", "histology_core_side", "assignment_status"]] = [circle["disc_id"], "UNRESOLVED", "UNRESOLVED", "ARTIFACT_EXCLUDED_WITH_BIN_MARGIN"]
    result["annotation_source"] = "new_frozen_provisional_H&E_core_polygons"
    result["expert_histopathology_verification"] = False
    return result


def marker_summary(x, obs, var, capture):
    lookup = dict(zip(var.mouse_gene_symbol, range(len(var))))
    qc = obs.pass_basic_qc.astype(str).str.lower().isin(["true", "1"]).to_numpy()
    library = np.asarray(x.sum(axis=1)).ravel().astype(float)
    rows = []
    for disc in capture["discs"]:
        for region in ["NP_CORE", "AF_CORE"]:
            morphology = obs.disc_id.eq(disc["disc_id"]).to_numpy() & obs.histology_core_region.eq(region).to_numpy()
            take = morphology & qc
            for marker in NP_MARKERS + AF_MARKERS:
                count = int(take.sum())
                present = marker in lookup
                row = dict(capture_id=capture["capture_id"], source_GSM=capture["source_GSM"], genotype=capture["genotype"],
                           disc_id=disc["disc_id"], histology_core_region=region, marker=marker,
                           marker_panel="NP" if marker in NP_MARKERS else "AF", feature_present=present,
                           n_geometry_bins=int(morphology.sum()), n_qc_bins=count, n_detected_bins=0,
                           detected_bin_fraction=0.0, total_raw_UMI=0, mean_log1p_per_10000=0.0)
                if present and count:
                    column = x[take, lookup[marker]].tocoo()
                    normalized = np.log1p(column.data.astype(float) / np.maximum(library[take][column.row], 1) * 10000)
                    row.update(n_detected_bins=int(np.count_nonzero(column.data)),
                               detected_bin_fraction=float(np.count_nonzero(column.data) / count),
                               total_raw_UMI=int(column.data.sum()), mean_log1p_per_10000=float(normalized.sum() / count))
                rows.append(row)
    return pd.DataFrame(rows)


def validate_markers(table, capture):
    output = []
    for disc in capture["discs"]:
        t = table.loc[table.disc_id == disc["disc_id"]]
        counts = t.groupby("histology_core_region").n_qc_bins.first()
        scores = t.groupby(["histology_core_region", "marker_panel"]).mean_log1p_per_10000.mean()
        row = dict(capture_id=capture["capture_id"], source_GSM=capture["source_GSM"], genotype=capture["genotype"],
                   disc_id=disc["disc_id"], NP_n_qc_bins=int(counts.get("NP_CORE", 0)), AF_n_qc_bins=int(counts.get("AF_CORE", 0)),
                   NP_panel_mean_NP_core=float(scores.get(("NP_CORE", "NP"), 0)),
                   NP_panel_mean_AF_core=float(scores.get(("AF_CORE", "NP"), 0)),
                   AF_panel_mean_NP_core=float(scores.get(("NP_CORE", "AF"), 0)),
                   AF_panel_mean_AF_core=float(scores.get(("AF_CORE", "AF"), 0)))
        for panel, own, other in [("NP", "NP_CORE", "AF_CORE"), ("AF", "AF_CORE", "NP_CORE")]:
            panel_table = t.loc[t.marker_panel == panel]
            own_table = panel_table.loc[panel_table.histology_core_region == own]
            row[f"{panel}_panel_features_present"] = int(own_table.feature_present.sum())
            row[f"{panel}_panel_genes_detected_in_own_core"] = int((own_table.n_detected_bins > 0).sum())
            if row["NP_n_qc_bins"] == 0 or row["AF_n_qc_bins"] == 0:
                status = "UNTESTABLE_NO_QC_BINS_IN_ONE_COMPARTMENT"
            elif int(panel_table.n_detected_bins.sum()) == 0:
                status = "UNTESTABLE_NO_INDEPENDENT_MARKER_SIGNAL"
            elif scores.get((own, panel), 0) > scores.get((other, panel), 0):
                status = "COMPATIBLE_DIRECTIONAL_MARKER_PATTERN"
            else:
                status = "UNRESOLVED_NONDISCRIMINATING_OR_OPPOSITE_MARKER_PATTERN"
            row[f"{panel}_marker_status"] = status
        row["interpretation"] = "Descriptive structural-marker check; no P value or independent-capture inference"
        output.append(row)
    return pd.DataFrame(output)


def process_capture(prepared, capture, out, config_sha):
    folder = prepared / capture["capture_id"]
    image = checked_image(capture, folder / "source_images" / capture["source_image_filename"])
    scale = json.loads((folder / "ORIGINAL_SCALEFACTORS.json").read_text())
    mpp = float(scale["microns_per_pixel"])
    if not np.isclose(mpp, capture["microns_per_pixel"], rtol=1e-9):
        raise ValueError("Frozen annotation pixel scale mismatch")
    margin = np.sqrt(2) * 48 / (2 * mpp)
    original = pd.read_csv(folder / "BIN48_METADATA_ALL_COMPLETE.tsv.gz", sep="\t", keep_default_na=False)
    var = pd.read_csv(folder / "MOUSE_SYMBOL_FEATURES.tsv", sep="\t", keep_default_na=False)
    x = sparse.load_npz(folder / "counts48_symbols_all_complete.npz").tocsr()
    assert x.shape == (len(original), len(var)) and original.spot_id.is_unique and var.mouse_gene_symbol.is_unique
    obs = assign_cores(original, capture, margin)
    table = marker_summary(x, obs, var, capture)
    validation = validate_markers(table, capture)
    status_lookup = {(r.disc_id, region): getattr(r, f"{panel}_marker_status") for r in validation.itertuples()
                     for panel, region in [("NP", "NP_CORE"), ("AF", "AF_CORE")]}
    obs["independent_marker_status"] = [status_lookup.get((d, r), "NOT_APPLICABLE") for d, r in zip(obs.disc_id, obs.histology_core_region)]
    is_core = obs.histology_core_region.isin(["NP_CORE", "AF_CORE"])
    compatible = obs.independent_marker_status.eq("COMPATIBLE_DIRECTIONAL_MARKER_PATTERN")
    obs["anatomical_region"] = "UNASSIGNED"
    obs.loc[obs.histology_core_region.eq("UNRESOLVED"), "anatomical_region"] = "UNRESOLVED"
    obs.loc[is_core & ~compatible, "anatomical_region"] = "UNRESOLVED"
    obs.loc[is_core & compatible, "anatomical_region"] = obs.loc[is_core & compatible, "histology_core_region"]
    obs["frozen_polygon_sha256"] = config_sha
    plot_frozen(capture, image, out)
    return obs, table, validation


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prepared48", type=Path)
    ap.add_argument("--polygons", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--preview-only", action="store_true")
    ap.add_argument("--image-root", type=Path)
    a = ap.parse_args()
    if not a.preview_only and not a.prepared48:
        ap.error("--prepared48 is required for real bin assignment/marker validation")
    config_sha = sha(a.polygons)
    config = json.loads(a.polygons.read_text(encoding="utf8"))
    assert config["candidate_and_program_expression_read"] is False
    assert config["validation_markers"] == {"NP": NP_MARKERS, "AF": AF_MARKERS}
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    metadata, markers, checks = [], [], []
    if a.preview_only:
        if not a.image_root:
            ap.error("--image-root is required for H&E-only preview")
        for capture in config["captures"]:
            images = list(a.image_root.rglob(capture["source_image_filename"]))
            if len(images) != 1:
                raise ValueError(f"Expected one exact original image, got {len(images)}")
            plot_frozen(capture, checked_image(capture, images[0]), out)
        print("HE_ONLY_PREVIEW_COMPLETE_NO_MATRIX_OR_MARKER_READ", flush=True)
        return
    for capture in config["captures"]:
        obs, table, validation = process_capture(a.prepared48, capture, out, config_sha)
        metadata.append(obs)
        markers.append(table)
        checks.append(validation)
    all_metadata = pd.concat(metadata, ignore_index=True)
    all_metadata.to_csv(out / "DISC_AF_NP_ROI_METADATA.tsv.gz", sep="\t", index=False)
    all_metadata.to_csv(out / "DISC_AF_NP_ROI_METADATA.tsv", sep="\t", index=False)
    pd.concat(markers, ignore_index=True).to_csv(out / "INDEPENDENT_MARKERS_BY_DISC_REGION.tsv", sep="\t", index=False)
    validation = pd.concat(checks, ignore_index=True)
    validation.to_csv(out / "INDEPENDENT_MARKER_VALIDATION_BY_DISC.tsv", sep="\t", index=False)
    dump(out / "SPATIAL_ROI_ASSIGNMENT_STATUS.json", {
        "status": "DONE_CONFIGURED_CAPTURES_ONLY", "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "command": sys.argv, "configured_capture_count": len(config["captures"]),
        "configured_disc_count": len(validation), "polygons_sha256": config_sha,
        "NP_compatible_discs": int(validation.NP_marker_status.eq("COMPATIBLE_DIRECTIONAL_MARKER_PATTERN").sum()),
        "AF_compatible_discs": int(validation.AF_marker_status.eq("COMPATIBLE_DIRECTIONAL_MARKER_PATTERN").sum()),
        "program_or_candidate_genes_read": False, "polygons_edited_after_marker_analysis": False,
        "expert_histopathology_verification": False, "genotype_significance_testing": False,
    })
    assert sha(a.polygons) == config_sha, "Frozen polygons changed during execution"
    print(validation.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
