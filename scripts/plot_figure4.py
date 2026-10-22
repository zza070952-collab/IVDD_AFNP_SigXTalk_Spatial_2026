#!/usr/bin/env python3
"""Figure 4: frozen spatial context and the single prespecified depth check.

All inputs are explicit. Maps and candidate summaries are unadjusted frozen
measurements; only panel e uses the completed fixed-count sensitivity analysis.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Circle, Polygon
from matplotlib.lines import Line2D
from PIL import Image

PROGRAM = "GEP2"
GENES = ["FN1", "ITGB1", "PRSS23"]
REGIONS = ["NP_CORE", "AF_CORE"]
CORE_PALETTE = {"NP_CORE": "#286EAE", "AF_CORE": "#C67922"}
CAPTURE_PALETTE = {"CTR": "#4477AA", "cKO": "#CC8844"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def draw_outlines(ax, capture, numbers=True, linewidth=.65, fontsize=5.5):
    margin = np.sqrt(2) * 48 / (2 * capture["microns_per_pixel"])
    for number, disc in enumerate(capture["discs"], 1):
        for region in disc["regions"]:
            ax.add_patch(Polygon(region["vertices"], closed=True, fill=False,
                                 edgecolor=CORE_PALETTE[region["region"]], linewidth=linewidth))
        if numbers:
            vertices = np.asarray(disc["regions"][0]["vertices"])
            ax.text(vertices[:, 0].mean(), vertices[:, 1].min() - 20, str(number),
                    fontsize=fontsize, ha="center", va="bottom", color="#222222",
                    bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none", "pad": .5})
    for circle in capture.get("exclusion_circles", []):
        ax.add_patch(Circle(circle["center"], circle["radius_pixels"] + margin, fill=False,
                            edgecolor="#B43C39", linewidth=linewidth, linestyle="--"))


def map_axes(ax, capture, image=None, scalebar=True):
    if image is not None:
        ax.imshow(image, origin="upper")
    x0, x1, y0, y1 = capture["overview_display_window_pixels"]
    ax.set_xlim(x0, x1)
    ax.set_ylim(y1, y0)
    ax.set_aspect("equal", adjustable="box", anchor="N")
    ax.axis("off")
    if scalebar:
        length = 500 / capture["microns_per_pixel"]
        start = x0 + .07 * (x1 - x0)
        y = y1 - .055 * (y1 - y0)
        ax.plot([start, start + length], [y, y], color="#222222", lw=1.0, solid_capstyle="butt")
        ax.text(start + length / 2, y - .012 * (y1 - y0), "500 µm", fontsize=5.5, ha="center", va="bottom")


def original_image(prepared, capture):
    choices = [prepared / capture["capture_id"] / "source_images" / capture["source_image_filename"],
               prepared / capture["source_image_filename"],
               prepared / capture["capture_id"] / capture["source_image_filename"]]
    path = next((p for p in choices if p.is_file()), None)
    if path is None:
        raise FileNotFoundError(f"Original H&E image missing: {capture['source_image_filename']}")
    if sha(path) != capture["source_image_sha256"]:
        raise ValueError(f"Frozen source-image hash mismatch: {capture['capture_id']}")
    with Image.open(path) as im:
        if list(im.size) != capture["source_image_pixels"]:
            raise ValueError("Source image dimensions differ from the frozen ROI record")
        return np.asarray(im).copy()


def independent_overlay(capture, image, out):
    fig, ax = plt.subplots(figsize=(5.0, 7.6), layout="constrained")
    map_axes(ax, capture, image)
    draw_outlines(ax, capture, linewidth=1, fontsize=8)
    ax.set_title(f"{capture['genotype']} capture: frozen H&E compartment cores\nBlue: NP; ochre: AF; numbers: disc positions", fontsize=9)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out / f"Figure4_{capture['genotype']}_HE_ROI_Overlay.{ext}", dpi=240, bbox_inches="tight")
    plt.close(fig)


def source_vertices(config):
    rows = []
    for capture in config["captures"]:
        for number, disc in enumerate(capture["discs"], 1):
            for region in disc["regions"]:
                for index, (x, y) in enumerate(region["vertices"]):
                    rows.append(dict(capture_id=capture["capture_id"], genotype=capture["genotype"],
                                     disc_number_in_capture=number, disc_id=disc["disc_id"],
                                     histology_region=region["region"], flank=region["side"],
                                     vertex_index=index, x_fullres=x, y_fullres=y,
                                     microns_per_pixel=capture["microns_per_pixel"],
                                     source_image_sha256=capture["source_image_sha256"]))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spatial-input", required=True, type=Path,
                    help="Read-only original program/map/candidate summaries")
    ap.add_argument("--depth-input", required=True, type=Path,
                    help="Completed single fixed-count sensitivity results")
    ap.add_argument("--roi-json", required=True, type=Path)
    ap.add_argument("--image-root", type=Path, help="Original image root; default spatial-input/prepared48")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    spatial = a.spatial_input.resolve()
    out = a.out.resolve()
    source_out = out / "source_data"
    source_out.mkdir(parents=True, exist_ok=True)
    polygons_path = a.roi_json.resolve()
    config = json.loads(polygons_path.read_text(encoding="utf8"))
    captures = sorted(config["captures"], key=lambda r: 0 if r["genotype"] == "CTR" else 1)
    if len(captures) != 2 or {r["genotype"] for r in captures} != {"CTR", "cKO"}:
        raise ValueError("Figure 4 requires the actual CTR and cKO captures")
    source_path = spatial / "SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz"
    program_path = spatial / "NP_PROGRAM_BY_REGION.tsv"
    gene_path = spatial / "CANDIDATE_GENES_BY_REGION.tsv"
    data = pd.read_csv(source_path, sep="\t", keep_default_na=False)
    programs = pd.read_csv(program_path, sep="\t")
    genes = pd.read_csv(gene_path, sep="\t")
    score_column = f"{PROGRAM}_mean_log1p_CP10k"
    required = {"spot_id", "capture_id", "genotype", "disc_id", "anatomical_region", "x_fullres", "y_fullres", score_column}
    if not required.issubset(data.columns) or not data.spot_id.is_unique:
        raise ValueError("Spatial source schema or spot identity check failed")
    score = pd.to_numeric(data[score_column], errors="coerce").to_numpy()
    if not np.isfinite(score).all() or (score < 0).any():
        raise ValueError("Program display requires finite actual nonnegative scores for all QC bins")
    program_norm = Normalize(vmin=0, vmax=max(float(score.max()), 1e-12))
    selected = programs.loc[programs.program.eq(PROGRAM) & programs.region.isin(REGIONS)].copy()
    if selected.duplicated(["capture_id", "disc_id", "region"]).any():
        raise ValueError("Repeated disc-region program summary")
    if selected.empty or not np.isfinite(selected.mean_program_score).all():
        raise ValueError("Actual GEP2 disc-region summary is unavailable")
    pair = selected.pivot(index=["capture_id", "genotype", "disc_id"], columns="region", values="mean_program_score").reset_index()
    if not set(REGIONS).issubset(pair.columns) or pair[REGIONS].isna().any().any():
        raise ValueError("A disc lacks an actual NP or AF program summary; do not create an artificial pair")
    depth_path = a.depth_input / "DEPTH_DISC_PAIR_COMPARISON.tsv"
    depth = pd.read_csv(depth_path, sep="\t")
    depth_status = json.loads((a.depth_input / "DEPTH_ANALYSIS_STATUS.json").read_text())
    if depth_status["ROI_sha256"] != sha(polygons_path):
        raise ValueError("Depth analysis used different frozen ROIs")
    depth["adequate_fixed_depth_pair"] = depth.adequate_fixed_depth_pair.astype(str).str.lower().eq("true")
    if depth.duplicated(["capture_id", "disc_id"]).any():
        raise ValueError("Duplicate disc in depth comparison")
    check = pair.merge(depth, on=["capture_id", "genotype", "disc_id"], validate="one_to_one")
    if len(check) != len(pair) or not np.allclose(check.NP_CORE-check.AF_CORE,
                                               check.original_all_NP_minus_AF, atol=1e-9, rtol=0):
        raise ValueError("Depth comparison does not reproduce original disc summaries")
    adequate = depth.loc[depth.adequate_fixed_depth_pair].copy()
    if len(adequate) != int(depth_status["adequate_paired_discs"]):
        raise ValueError("Depth adequacy counts disagree with completed analysis")
    if adequate.empty:
        raise ValueError("No adequate paired discs: do not draw a depth-confirmation panel")
    subset = genes.loc[genes.human_gene.isin(GENES) & genes.region.isin(REGIONS)].copy()
    subset["native_feature_assessable"] = subset.native_feature_assessable.astype(str).str.lower().isin(["true", "1"])
    if subset.duplicated(["capture_id", "disc_id", "region", "human_gene"]).any():
        raise ValueError("Repeated candidate disc-region summary")
    expected_keys = [(c["capture_id"], region, gene) for c in captures for region in REGIONS for gene in GENES]
    grouped = subset.groupby(["capture_id", "genotype", "region", "human_gene"], sort=False).agg(
        disc_equal_mean_log1p_CP10k=("mean_log1p_CP10k", "mean"),
        disc_equal_mean_fraction_detected=("fraction_detected_bins", "mean"),
        n_discs=("disc_id", "nunique"),
        n_assessable_discs=("native_feature_assessable", "sum"),
        mouse_symbols=("mouse_symbol", lambda z: ";".join(sorted(set(z.dropna().astype(str))))),
    ).reset_index()
    lookup = grouped.set_index(["capture_id", "region", "human_gene"])
    if not all(key in lookup.index for key in expected_keys):
        raise ValueError("Missing real candidate-region summaries")
    values = pd.to_numeric(grouped.disc_equal_mean_log1p_CP10k, errors="coerce").to_numpy()
    finite = values[np.isfinite(values)]
    if not len(finite):
        raise ValueError("None of the requested candidate orthologs is assessable")
    gene_norm = Normalize(vmin=0, vmax=max(float(finite.max()), 1e-12))

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": .7, "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False})
    fig = plt.figure(figsize=(8.4, 8.0))
    xs = [.045, .285, .525, .765]
    top_axes = [fig.add_axes([x, .555, .20, .385]) for x in xs]
    capture_counts, display_counts = {}, {}
    data["in_figure_view_window"] = False
    for index, capture in enumerate(captures):
        image = original_image(a.image_root or spatial / "prepared48", capture)
        independent_overlay(capture, image, out)
        cdata = data.loc[data.capture_id == capture["capture_id"]]
        if cdata.empty:
            raise ValueError(f"Missing actual spatial rows for {capture['capture_id']}")
        width, height = capture["source_image_pixels"]
        if not (cdata.x_fullres.between(0, width).all() and cdata.y_fullres.between(0, height).all()):
            raise ValueError("Outside-H&E QC bins must be flagged/excluded in preparation, not moved during plotting")
        x0, x1, y0, y1 = capture["overview_display_window_pixels"]
        displayed = cdata.x_fullres.between(x0, x1) & cdata.y_fullres.between(y0, y1)
        data.loc[cdata.index, "in_figure_view_window"] = displayed.to_numpy()
        he_ax, program_ax = top_axes[index], top_axes[index + 2]
        map_axes(he_ax, capture, image)
        draw_outlines(he_ax, capture)
        he_ax.set_title(f"{chr(97+index)}  {capture['genotype']} H&E + cores", loc="left", pad=6, fontsize=8)
        map_axes(program_ax, capture)
        # Display every real QC-bin center. This is neither smoothing nor an
        # interpolation of values into empty space or neighboring tissue.
        program_ax.scatter(cdata.x_fullres, cdata.y_fullres, c=cdata[score_column], cmap="viridis",
                           norm=program_norm, s=.75, marker="s", linewidths=0, rasterized=True)
        draw_outlines(program_ax, capture, linewidth=.5, fontsize=5)
        program_ax.set_title(f"{chr(99+index)}  {capture['genotype']} GEP2 (original)", loc="left", pad=6, fontsize=8)
        capture_counts[capture["capture_id"]] = len(cdata)
        display_counts[capture["capture_id"]] = int(displayed.sum())
    cax = fig.add_axes([.565, .515, .355, .012])
    cb = fig.colorbar(ScalarMappable(norm=program_norm, cmap="viridis"), cax=cax, orientation="horizontal")
    cb.set_label("GEP2 score: mean log(1 + CP10k)", fontsize=7, labelpad=2)
    cb.ax.tick_params(labelsize=6, length=2)
    fig.legend(handles=[Line2D([0], [0], color=CORE_PALETTE["NP_CORE"], lw=1, label="NP core"),
                        Line2D([0], [0], color=CORE_PALETTE["AF_CORE"], lw=1, label="AF core")],
               loc="upper left", bbox_to_anchor=(.055, .54), ncol=2, fontsize=7, handlelength=1.5, columnspacing=1)
    fig.text(.055, .467, "Frozen H&E cores; structural-marker compatibility. One capture per genotype.", fontsize=7)

    ax_e = fig.add_axes([.10, .135, .345, .265])
    jitter = np.linspace(-.055, .055, len(depth)) if len(depth) > 1 else [0]
    plotted_pairs = []
    stages = ["original_all_NP_minus_AF", "original_eligible_NP_minus_AF", "depth100_NP_minus_AF"]
    for offset, row in zip(jitter, depth.itertuples(index=False)):
        color = CAPTURE_PALETTE[row.genotype]
        if row.adequate_fixed_depth_pair:
            y = [getattr(row, c) for c in stages]
            ax_e.plot(np.arange(3)+offset, y, color=color, alpha=.72, linewidth=.85,
                      marker="o", markersize=3.8, markeredgecolor="white", markeredgewidth=.35)
        else:
            ax_e.scatter([offset], [row.original_all_NP_minus_AF], s=20,
                         facecolors="white", edgecolors="#777777", linewidths=.8, zorder=4)
        for stage_index, column in enumerate(stages):
            plotted_pairs.append(dict(capture_id=row.capture_id, genotype=row.genotype, disc_id=row.disc_id,
                                      stage=column, contrast_NP_minus_AF=getattr(row, column),
                                      adequate_fixed_depth_pair=row.adequate_fixed_depth_pair,
                                      displayed=stage_index == 0 or row.adequate_fixed_depth_pair,
                                      NP_retained_bins=row.NP_retained_bins, AF_retained_bins=row.AF_retained_bins,
                                      display_x_offset=offset))
    means = adequate[stages].mean().to_numpy()
    ax_e.plot(range(3), means, color="#222222", lw=1.2, marker="D", ms=4,
              markeredgecolor="white", markeredgewidth=.35, zorder=5)
    ax_e.axhline(0, color="#777777", lw=.6, ls="--", zorder=0)
    ax_e.set_xticks([0, 1, 2])
    ax_e.set_xticklabels([f"Original\n{len(depth)} discs", f"Eligible bins\n{len(adequate)} discs",
                         f"100 UMI\n{len(adequate)} discs"], fontsize=7)
    ax_e.set_xlim(-.18, 2.18)
    ax_e.set_ylabel("GEP2 NP − AF contrast\nmean log(1 + CP10k)")
    ax_e.set_title("e  NP contrast weakens at fixed depth", loc="left", pad=9)
    for capture in captures:
        n = int(adequate.capture_id.eq(capture["capture_id"]).sum())
        ax_e.plot([], [], "o-", color=CAPTURE_PALETTE[capture["genotype"]], ms=4, lw=.85,
                  label=f"{capture['genotype']}: {n} adequate pairs")
    ax_e.plot([], [], "D-", color="#222222", ms=4, lw=1, label=f"Mean, same {len(adequate)} discs")
    ax_e.legend(fontsize=6, loc="upper right", handlelength=1.4, borderaxespad=.15)
    n_positive = int((adequate.depth100_NP_minus_AF > 0).sum())
    fig.text(.10, .047, f"Same eligible bins: {means[1]:.3f} → {means[2]:.3f}; NP > AF: {n_positive}/{len(adequate)}.\n"
                         "Open point: original disc with insufficient AF retention.", fontsize=6.4, va="bottom")

    ax_f = fig.add_axes([.565, .145, .32, .255])
    columns = [(capture, region) for capture in captures for region in REGIONS]
    plotted_dots = []
    for col, (capture, region) in enumerate(columns):
        for row_index, gene in enumerate(GENES):
            row = lookup.loc[(capture["capture_id"], region, gene)]
            value = float(row.disc_equal_mean_log1p_CP10k)
            detection = float(row.disc_equal_mean_fraction_detected)
            if np.isfinite(detection) and not 0 <= detection <= 1:
                raise ValueError("Candidate detection fraction outside [0,1]")
            if not np.isfinite(value) or not np.isfinite(detection):
                ax_f.scatter(col, row_index, marker="x", s=24, c="#999999", linewidths=.8)
                display = "not_assessable"
            elif detection == 0:
                ax_f.scatter(col, row_index, s=15, facecolor="white", edgecolor="#999999", linewidths=.7)
                display = "zero_detected_hollow"
            else:
                ax_f.scatter(col, row_index, s=15 + 190 * detection, c=[value], cmap="magma",
                             norm=gene_norm, linewidths=.35, edgecolors="white")
                display = "actual_mean_and_detection"
            plotted_dots.append(dict(capture_id=capture["capture_id"], genotype=capture["genotype"], region=region,
                                     human_gene=gene, mouse_symbols=row.mouse_symbols,
                                     disc_equal_mean_log1p_CP10k=value, disc_equal_mean_fraction_detected=detection,
                                     n_discs=int(row.n_discs), n_assessable_discs=int(row.n_assessable_discs),
                                     display_status=display))
    ax_f.set_xticks(range(len(columns)))
    ax_f.set_xticklabels([c["genotype"] + "\n" + r.replace("_CORE", "") for c, r in columns])
    ax_f.set_yticks(range(len(GENES)))
    ax_f.set_yticklabels(GENES, fontstyle="italic")
    ax_f.set_xlim(-.5, len(columns) - .5)
    ax_f.set_ylim(len(GENES) - .5, -.5)
    ax_f.set_title("f  Candidate orthologs (unadjusted)", loc="left", pad=9)
    ax_f.tick_params(length=0, pad=6)
    for spine in ax_f.spines.values():
        spine.set_visible(False)
    gcax = fig.add_axes([.915, .18, .012, .175])
    gcb = fig.colorbar(ScalarMappable(norm=gene_norm, cmap="magma"), cax=gcax)
    gcb.set_label("Mean log(1 + CP10k)", fontsize=7, labelpad=4)
    gcb.ax.tick_params(labelsize=6, length=2)
    for detection in [.1, .5, 1.0]:
        ax_f.scatter([], [], s=15 + 190 * detection, color="#777777", label=f"{100*detection:.0f}%")
    ax_f.legend(title="Detected bins (disc mean)", title_fontsize=6.5, fontsize=6.5,
                bbox_to_anchor=(.5, -.22), loc="upper center", ncol=3, columnspacing=.9,
                handletextpad=.4, borderaxespad=0)
    fig.text(.565, .022, "Hollow: 0 detected; ×: not assessable", fontsize=6.5)
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(out / f"Figure4_Adult_Mouse_Spatial_Context.{ext}", dpi=300)
    plt.close(fig)

    source_vertices(config).to_csv(source_out / "Figure4ab_frozen_HE_vertices.tsv", sep="\t", index=False)
    spatial_columns = [c for c in ["spot_id", "capture_id", "genotype", "disc_id", "histology_core_region",
                                  "anatomical_region", "x_fullres", "y_fullres", "n_counts", "n_genes", score_column,
                                  "in_figure_view_window"] if c in data]
    data[spatial_columns].to_csv(source_out / "Figure4cd_actual_spatial_program_coordinates.tsv.gz", sep="\t", index=False)
    selected.to_csv(source_out / "Figure4e_original_disc_region_program_summary.tsv", sep="\t", index=False)
    depth.to_csv(source_out / "Figure4e_depth_disc_pair_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(plotted_pairs).to_csv(source_out / "Figure4e_connected_depth_points.tsv", sep="\t", index=False)
    pd.DataFrame({"stage": stages, "mean_NP_minus_AF_same_adequate_discs": means,
                  "n_adequate_discs": len(adequate)}).to_csv(source_out / "Figure4e_disc_equal_depth_means.tsv", sep="\t", index=False)
    subset.to_csv(source_out / "Figure4f_actual_candidate_disc_region_values.tsv", sep="\t", index=False)
    pd.DataFrame(plotted_dots).to_csv(source_out / "Figure4f_disc_equal_candidate_dots.tsv", sep="\t", index=False)
    disc_counts = "; ".join(f"{c['genotype']}: {pair.capture_id.eq(c['capture_id']).sum()} discs, {capture_counts[c['capture_id']]:,} QC bins" for c in captures)
    coverage = selected.groupby("capture_id").n_strict_1to1_native_members.unique().to_dict()
    coverage_text = "; ".join(f"{key}: {','.join(str(int(v)) for v in values)}" for key, values in coverage.items())
    legend = f"""# Figure 4 | Adult mouse disc spatial context of a fixed human NP program and candidate orthologs

**a,b**, Genuine original H&E images from GSE287499 with previously frozen NP-core (blue) and bilateral AF-core (ochre) polygons. Numbers identify disc positions within each capture and are not released individual animal identifiers. Only centers at least a 48-µm bin half diagonal from the frozen edges contribute to compartment summaries. Ambiguous interfaces, CEP/growth plate and artifact neighborhoods were left unassigned. Source pixels were not altered. Independent Krt19/T/Krt8/Acan and Chad/Clec3a panels were assessed after polygon freezing; Fn1 and Comp were excluded from anatomical validation. Scale bars, 500 µm, use each capture's own original pixel scale.

**c,d**, Actual QC 48-µm bins show the original, unadjusted GEP2 scores at their original coordinates with one common color scale. No interpolation or smoothing was used. These maps do not display the depth-controlled analysis in panel e. View windows contain all four disc structures and match the H&E panels; bins outside these windows remain in the underlying data, color-scale calculation and source TSV. Scores are an unweighted mean log(1+CP10k) over fixed available strict 1:1 ortholog members, with zeros retained. Member counts: {coverage_text}. Full source bin counts are technical coverage, not independent animals: {disc_counts}.

**e**, Disc-level NP-minus-AF GEP2 contrasts before and after the single predeclared depth sensitivity analysis. The first column includes all {len(depth)} original discs. Connected lines follow the same {len(adequate)} adequately covered discs through the original-all, original-eligible and exactly-100-UMI scores; black diamonds are equal-disc means over those same {len(adequate)} discs at every stage. Only original core bins with at least 100 total UMIs were eligible. One exact multivariate-hypergeometric draw per bin sampled 100 molecules without replacement across all genes, with fixed seed 20261022, before using the same frozen program members and normalization. Adequacy required at least five eligible bins per compartment and at least three paired discs per capture. The open original-only point has insufficient AF-bin retention and is retained visibly; its exploratory sampled estimate is present in the full table but excluded from adequate-pair summaries. The eligible-bin mean contrast fell from {means[1]:.4f} to {means[2]:.4f}; NP exceeded AF in {n_positive}/{len(adequate)} adequate disc pairs. This attenuation weakens the original compartment-localization claim. One draw at one depth is a sensitivity check, not an uncertainty interval. Discs are nested in capture; genotype and capture are inseparable. No bin-level, animal-level or genotype P value is shown.

**f**, Original unadjusted FN1, ITGB1 and PRSS23 measurements label their explicitly mapped mouse orthologs. Color is the equal-disc-weight mean of each region's mean normalized log expression; area indicates the equal-disc-weight fraction of detected QC bins. Each disc receives equal weight. Hollow dots denote zero detected expression and crosses denote unavailable feature/mapping data. The full 11-gene table remains available. This panel describes detectability and was not depth-corrected; sparse Prss23 detection does not establish spatial validation.

This is exploratory compartment localization in an adult conditional-knockout mouse model, not a human spatial validation cohort, an unconfounded genotype test, or proof of directional AF→NP communication. Sources: [GSE287499](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287499) and [Tsingas et al.](https://doi.org/10.1038/s41419-026-08942-4). PNG/PDF/SVG and per-panel source TSVs accompany the figure.
"""
    (out / "Figure4_LEGEND.md").write_text(legend, encoding="utf8")
    status = {"status": "DONE_FIGURE4_FROM_REAL_STAGE_C_OUTPUTS", "finished_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "actual_capture_bin_counts": capture_counts, "actual_original_discs": len(pair),
              "actual_adequate_depth_pairs": len(adequate), "depth_NP_greater_AF_pairs": n_positive,
              "depth_original_same_eligible_mean_contrast": float(means[1]),
              "depth100_mean_contrast": float(means[2]), "input_depth_comparison_sha256": sha(depth_path),
              "actual_bins_within_display_windows": display_counts,
              "actual_candidate_dot_rows": len(plotted_dots), "common_program_vmin": 0,
              "common_program_vmax": float(program_norm.vmax), "candidate_color_vmax": float(gene_norm.vmax),
              "polygons_sha256": sha(polygons_path), "input_program_summary_sha256": sha(program_path),
              "input_gene_summary_sha256": sha(gene_path), "ROI_or_program_recomputation": False,
              "smoothing_or_imputation": False, "genotype_significance_tests": False}
    (out / "FIGURE4_RENDER_STATUS.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf8")
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
