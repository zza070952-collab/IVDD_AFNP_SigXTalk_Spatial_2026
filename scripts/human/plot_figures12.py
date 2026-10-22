#!/usr/bin/env python3
"""Plot FCDB Figures 1–2 solely from actual Stage A tables.

Missing inputs fail clearly or are explicitly labelled. This script never
creates substitute observations and never changes frozen gene signatures.
"""
from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PROGRAMS = [f"GEP{i}" for i in range(1, 6)]
COLORS = {"NP": "#527E87", "AF": "#B9A88C", "Tu": "#597D9A", "Swahn": "#9A705F"}
IDENTITIES = ["NP_matrix_marker_supported", "disc_stromal_marker_supported",
              "immune_marker_supported", "endothelial_marker_supported",
              "mixed_lineage_marker_support", "unresolved_low_marker_support"]
IDENTITY_NAMES = {
    "NP_matrix_marker_supported": "NP matrix support",
    "disc_stromal_marker_supported": "Disc stromal support",
    "immune_marker_supported": "Immune support",
    "endothelial_marker_supported": "Endothelial support",
    "mixed_lineage_marker_support": "Mixed marker support",
    "unresolved_low_marker_support": "Unresolved support",
    "NP_chondrocyte_like": "NP chondrocyte-like",
    "other_disc_stromal": "Other disc stromal",
}
MARKER_GENES = ["ACAN", "COL2A1", "DCN", "COL1A1", "PTPRC", "LST1", "PECAM1", "VWF", "ACTA2", "TAGLN"]


def configure():
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 8.5,
        "xtick.labelsize": 8, "ytick.labelsize": 8.5, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .7,
        "pdf.fonttype": 42, "svg.fonttype": "none", "legend.frameon": False,
        "savefig.facecolor": "white"})


def read(path, required):
    if not path.exists():
        raise FileNotFoundError(f"Required actual source table missing: {path}")
    table = pd.read_csv(path, sep="\t")
    missing = set(required) - set(table.columns)
    if missing:
        raise ValueError(f"{path.name} missing columns: {sorted(missing)}")
    if table.empty:
        raise ValueError(f"Actual source table is empty: {path}")
    return table


def write(frame, directory, name):
    frame.to_csv(directory / name, sep="\t", index=False)


def letter(ax, value):
    ax.text(-.085, 1.04, value, transform=ax.transAxes, fontsize=12,
            weight="bold", va="bottom", ha="left")


def wrap(value, width=27):
    return "\n".join(textwrap.wrap(str(value), width=width))


def donor_sort(value):
    text = str(value)
    return int(text.rsplit("_", 1)[-1]) if text.startswith("Tu_") else text


def truth(series):
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def export(fig, directory, name):
    fig.savefig(directory / (name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(directory / (name + ".pdf"), bbox_inches="tight")
    fig.savefig(directory / (name + ".svg"), bbox_inches="tight")
    plt.close(fig)


def donor_gene_summary(source, gene_column):
    required = {"donor_id", "identity", gene_column, "fraction_detected", "mean_log1p_published_expression"}
    if not required.issubset(source.columns):
        raise ValueError("Incomplete actual gene-expression source")
    # A marker may occur in several input panels; that duplication must not
    # alter the donor weight. Collapse its repeated donor/identity row first.
    per_donor = source.groupby(["donor_id", "identity", gene_column], as_index=False).agg(
        fraction_detected=("fraction_detected", "mean"),
        mean_log1p_published_expression=("mean_log1p_published_expression", "mean"))
    grouped = per_donor.groupby(["identity", gene_column], as_index=False).agg(
        fraction_detected=("fraction_detected", "mean"),
        mean_log1p_published_expression=("mean_log1p_published_expression", "mean"),
        n_contributing_donors=("fraction_detected", "count"))
    return grouped


def dotplot(ax, table, genes, identities, gene_column, orientation="genes_x"):
    table = table.set_index(["identity", gene_column])
    vals = table.mean_log1p_published_expression.to_numpy()
    finite = vals[np.isfinite(vals)]
    vmax = float(finite.max()) if finite.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    norm = colors.Normalize(0, vmax)
    xs, ys, sizes, color_values = [], [], [], []
    for gi, gene in enumerate(genes):
        for ii, identity in enumerate(identities):
            x, y = (gi, ii) if orientation == "genes_x" else (ii, gi)
            key = (identity, gene)
            if key not in table.index:
                ax.plot(x, y, marker="x", color="#BBBBBB", markersize=3.5, mew=.7)
                continue
            row = table.loc[key]
            mean = float(row.mean_log1p_published_expression)
            detect = float(row.fraction_detected)
            if not np.isfinite(mean) or not np.isfinite(detect):
                ax.plot(x, y, marker="x", color="#BBBBBB", markersize=3.5, mew=.7)
                continue
            if detect <= 0:
                # A measured zero is a tiny open point; unmeasured is a cross.
                ax.scatter([x], [y], s=5, facecolors="none", edgecolors="#888888", linewidths=.4)
                continue
            xs.append(x); ys.append(y); sizes.append(115 * detect); color_values.append(mean)
    scatter = ax.scatter(xs, ys, s=sizes, c=color_values, cmap="Blues", norm=norm,
                         edgecolors="#4C5E65", linewidths=.25)
    if orientation == "genes_x":
        ax.set_xticks(range(len(genes)))
        ax.set_xticklabels(genes, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_yticks(range(len(identities)))
        ax.set_yticklabels([wrap(IDENTITY_NAMES.get(x, x), 22) for x in identities])
        ax.set_xlim(-.6, len(genes) - .4)
        ax.set_ylim(len(identities) - .5, -.6)
    else:
        ax.set_xticks(range(len(identities)))
        ax.set_xticklabels([wrap(IDENTITY_NAMES.get(x, x), 15) for x in identities], rotation=55,
                           ha="right", rotation_mode="anchor")
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.set_xlim(-.6, len(identities) - .4)
        ax.set_ylim(len(genes) - .5, -.6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return scatter


def size_legend(ax, loc="upper left", anchor=(0, -.22), figure_anchor=None):
    handles = [ax.scatter([], [], s=115*x, facecolors="#8EAFC5", edgecolors="#4C5E65", linewidths=.3)
               for x in [.25, .5, 1.0]]
    container = ax.figure if figure_anchor is not None else ax
    container.legend(handles, ["25%", "50%", "100%"], title="Mean donor detection", ncol=3,
              loc="lower left" if figure_anchor is not None else loc,
              bbox_to_anchor=figure_anchor if figure_anchor is not None else anchor,
              borderaxespad=0, handletextpad=.6, columnspacing=1.0, title_fontsize=8)


def figure1(args, results, output):
    counts = read(results / "FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv",
                  ["cohort", "donor_id", "gsm", "identity", "n_cells"])
    marker = read(results / "TU_MARKER_DOTPLOT_SOURCE.tsv",
                  ["donor_id", "identity", "gene", "fraction_detected", "mean_log1p_published_expression"])
    donor_counts = counts.groupby(["cohort", "donor_id"], as_index=False).n_cells.sum()
    donor_counts["NP_retained_cells"] = donor_counts.n_cells
    donor_counts["AF_retained_cells"] = np.nan
    donor_counts["AF_available"] = False
    af_note = "Donor overview shows analysed NP-source cells; AF sender input is reported with the formal communication analysis."
    if args.af_donor_counts:
        af = read(args.af_donor_counts, ["donor_id", "AF_retained_cells", "NP_retained_cells"])
        if af.donor_id.duplicated().any():
            raise ValueError("Duplicate AF donor count records")
        swahn = donor_counts.cohort.eq("Swahn_GSE230808")
        for i, row in donor_counts.loc[swahn].iterrows():
            match = af[af.donor_id == row.donor_id]
            if len(match) != 1:
                raise ValueError("Missing verified AF pairing for " + row.donor_id)
            if int(match.iloc[0].NP_retained_cells) != int(row.NP_retained_cells):
                raise ValueError("Legacy AF/NP count source disagrees with actual discovery NP count for " + row.donor_id)
            donor_counts.loc[i, "AF_retained_cells"] = int(match.iloc[0].AF_retained_cells)
            donor_counts.loc[i, "AF_available"] = True
        af_note = "Swahn AF counts are verified legacy retained stromal cells from matched donors; Tu sampled NP only."
    sw = donor_counts[donor_counts.cohort.eq("Swahn_GSE230808")].sort_values("donor_id")
    tu = donor_counts[donor_counts.cohort.eq("Tu_GSE165722")].copy()
    tu["sort"] = tu.donor_id.map(donor_sort)
    donor_counts = pd.concat([sw, tu.sort_values("sort")], ignore_index=True).drop(columns=["sort"], errors="ignore")
    identity_counts = counts.groupby(["cohort", "identity"], as_index=False).n_cells.sum()
    identity_counts["label"] = identity_counts.apply(lambda r: ("Tu: " if r.cohort.startswith("Tu") else "Swahn: ") + IDENTITY_NAMES.get(r.identity, r.identity), axis=1)
    identity_counts = identity_counts.sort_values(["cohort", "n_cells"], ascending=[False, False])
    dot = donor_gene_summary(marker[marker.gene.isin(MARKER_GENES)], "gene")
    identities = [x for x in IDENTITIES if x in dot.identity.unique()]
    genes = [x for x in MARKER_GENES if x in dot.gene.unique()]
    if not identities or not genes:
        raise ValueError("No actual Tu marker measurements for Figure 1")
    write(donor_counts, output, "Figure1a_donor_tissue_counts.tsv")
    write(identity_counts, output, "Figure1b_NP_source_identity_counts.tsv")
    write(dot, output, "Figure1c_marker_donor_weighted_dotplot.tsv")

    fig = plt.figure(figsize=(7.1, 9.1))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.15, 1.05, 1.25], hspace=.75,
                          left=.25, right=.92, top=.95, bottom=.11)
    ax = fig.add_subplot(gs[0]); letter(ax, "a")
    x = np.arange(len(donor_counts))
    np_values = donor_counts.NP_retained_cells.to_numpy()
    af_values = donor_counts.AF_retained_cells.fillna(0).to_numpy()
    ax.bar(x, np_values, width=.72, color=COLORS["NP"], label="NP-source cells")
    if donor_counts.AF_available.any():
        ax.bar(x, af_values, bottom=np_values, width=.72, color=COLORS["AF"], label="Matched AF stromal cells")
    labels = [r.donor_id.replace("Tu_", "Tu ") for r in donor_counts.itertuples()]
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=55, ha="right", rotation_mode="anchor")
    ax.set_ylabel("Retained cells")
    ax.set_title("Adult tissue coverage by donor", loc="left", pad=16)
    ax.grid(axis="y", color="#E9E9E9", linewidth=.5); ax.set_axisbelow(True)
    sw_n = len(sw)
    if sw_n and len(tu):
        ax.axvline(sw_n-.5, color="#888888", linestyle="--", linewidth=.7)
        ax.text((sw_n-1)/2, 1.02, f"Swahn: {sw_n} donors", ha="center", transform=ax.get_xaxis_transform())
        ax.text(sw_n+(len(tu)-1)/2, 1.02, f"Tu: {len(tu)} donors", ha="center", transform=ax.get_xaxis_transform())
    ax.legend(ncol=2, loc="upper right", bbox_to_anchor=(1, 1.35), fontsize=8)
    ax.set_ylim(0, float((np_values+af_values).max())*1.15)

    ax = fig.add_subplot(gs[1]); letter(ax, "b")
    y = np.arange(len(identity_counts))
    color = [COLORS["Tu"] if str(c).startswith("Tu") else COLORS["Swahn"] for c in identity_counts.cohort]
    ax.barh(y, identity_counts.n_cells, color=color, height=.65)
    ax.set_yticks(y); ax.set_yticklabels([wrap(x, 29) for x in identity_counts.label])
    ax.invert_yaxis(); ax.set_xlabel("Cells in NP-source analysis")
    for yy, val in zip(y, identity_counts.n_cells):
        ax.text(val+identity_counts.n_cells.max()*.012, yy, f"{int(val):,}", va="center", fontsize=8)
    ax.set_xlim(0, identity_counts.n_cells.max()*1.23)
    ax.set_title("Discovery identities and Tu marker-support identities", loc="left", pad=10)
    ax.grid(axis="x", color="#E9E9E9", linewidth=.5); ax.set_axisbelow(True)

    ax = fig.add_subplot(gs[2]); letter(ax, "c")
    scatter = dotplot(ax, dot, genes, identities, "gene")
    ax.set_title("Tu NP tissue: marker context", loc="left", pad=10)
    cbar = fig.colorbar(scatter, ax=ax, fraction=.035, pad=.035, aspect=20)
    cbar.set_label("Mean donor log1p expression", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    size_legend(ax, figure_anchor=(.25, .002))
    fig.text(.25, -.045, "Cohort-specific QC and identity rules; counts describe coverage, not independent patients.", fontsize=8)
    export(fig, output, "Figure1_Adult_Cohorts_and_Identity")
    total_discovery = int(counts.loc[counts.cohort.eq("Swahn_GSE230808"), "n_cells"].sum())
    np_discovery = int(counts.loc[counts.cohort.eq("Swahn_GSE230808") & counts.identity.eq("NP_chondrocyte_like"), "n_cells"].sum())
    legend = (f"Fig. 1 | Adult cohort and cell-identity context. (a) Actual retained cell counts by donor and tissue source. {af_note} "
              f"(b) NP-source analysis populations are shown using their actual cohort-specific identity definitions. The Swahn population contains {total_discovery:,} broad stromal cells, of which {np_discovery:,} carry the discovery high-confidence NP chondrocyte-like identity. "
              "Tu labels reflect transparent marker co-detection support and include mixed and unresolved cells; they are not a definitive integrated cell atlas. "
              "(c) Tu marker-expression context, averaged equally across contributing donors within each identity. Dot area denotes mean donor detection fraction; colour denotes mean donor log1p of the published processed expression. Tiny open dots denote measured zero detection and grey crosses denote unavailable measurements. "
              "Each donor is a biological sample; cells measure coverage. No between-cohort disease-effect P values are shown. QC, source platforms and annotation rules differ between cohorts.")
    (output / "Figure1_LEGEND.md").write_text(legend+"\n", encoding="utf-8")
    return {"donors_Swahn": len(sw), "donors_Tu": len(tu), "Swahn_NP_stromal_cells": total_discovery,
            "Swahn_high_confidence_NP_cells": np_discovery, "AF_counts_used": bool(args.af_donor_counts)}


def figure2(args, results, output):
    scores = read(results / "PROGRAM_SCORE_DONOR_SUMMARY.tsv", ["cohort", "donor_id", "gsm", "identity", "program", "score_median", "n_cells"])
    frozen = read(results / "FROZEN_K5_PROGRAMS.tsv", ["program", "rank", "gene"])
    clinical = read(results / "COHORT_AND_DONOR_TABLE.tsv", ["gsm", "donor_id", "geo_grade_label", "paper_pfirrmann_if_number_crosswalk_valid"])
    gene_source = read(results / "TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv", ["donor_id", "identity", "program", "original_gene", "gene_available", "fraction_detected", "mean_log1p_published_expression"])
    coverage = read(results / "PROGRAM_GENE_COVERAGE.tsv", ["program", "n_mapped_genes"])
    primary = scores[scores.cohort.eq("Tu_GSE165722") & scores.identity.eq(args.tu_primary_identity)].copy()
    if primary.empty:
        raise ValueError("No actual Tu primary-identity donor scores; cannot substitute another identity silently")
    if "score_method" not in scores:
        raise ValueError("Explicit score_method is required; AUCell must never be labelled UCell")
    tu_methods = primary.score_method.unique().tolist()
    if len(tu_methods) != 1:
        raise ValueError("Multiple Tu score methods in the actual source table")
    tu_unit = "AUCell" if "AUCell" in tu_methods[0] else "UCell"
    clinical["sort"] = clinical.donor_id.map(donor_sort)
    clinical = clinical.sort_values("sort")
    donors = clinical.donor_id.tolist()
    if primary.duplicated(["program", "donor_id"]).any():
        raise ValueError("Multiple scores for a primary donor/program")
    selected_genes = frozen[frozen.program.eq("GEP2")].sort_values("rank").gene.tolist()
    if len(selected_genes) != 25 or len(set(selected_genes)) != 25:
        raise ValueError("The exact frozen GEP2 Top25 list is required")
    selected_source = gene_source[gene_source.program.eq("GEP2")].copy()
    # Unavailable genes remain NaN, never converted into zero-expression evidence.
    missing_gene = ~truth(selected_source.gene_available)
    selected_source.loc[missing_gene, ["fraction_detected", "mean_log1p_published_expression"]] = np.nan
    gene_dot = donor_gene_summary(selected_source, "original_gene")
    identities = [x for x in IDENTITIES if x in gene_dot.identity.unique()]
    context = scores[scores.program.eq("GEP2")].copy()
    context["context_label"] = context.apply(lambda r: ("Tu: " if str(r.cohort).startswith("Tu") else "Swahn: ") + IDENTITY_NAMES.get(r.identity, r.identity), axis=1)
    write(primary, output, "Figure2a_Tu_fixed_program_donor_medians.tsv")
    write(clinical.drop(columns="sort"), output, "Figure2a_reported_and_putative_grade_labels.tsv")
    write(gene_dot, output, "Figure2b_GEP2_original25_donor_weighted_expression.tsv")
    write(context, output, "Figure2c_GEP2_cohort_identity_donor_medians.tsv")
    write(coverage, output, "Figure2_program_gene_coverage.tsv")

    fig = plt.figure(figsize=(7.1, 9.1))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 2.1], width_ratios=[1.2, 1],
                          hspace=.65, wspace=.9, left=.14, right=.95, top=.95, bottom=.16)
    ax = fig.add_subplot(gs[0, :]); letter(ax, "a")
    indexed = primary.set_index(["program", "donor_id"])
    xx, yy, zz = [], [], []
    for pi, program in enumerate(PROGRAMS):
        for di, donor in enumerate(donors):
            if (program, donor) not in indexed.index:
                ax.plot(di, pi, marker="x", color="#AAAAAA", markersize=5)
                continue
            value = float(indexed.loc[(program, donor), "score_median"])
            if not np.isfinite(value):
                ax.plot(di, pi, marker="x", color="#AAAAAA", markersize=5)
                continue
            xx.append(di); yy.append(pi); zz.append(value)
            ax.text(di, pi+.27, f"{value:.2f}", ha="center", va="top", fontsize=8, color="#46535B")
    dots = ax.scatter(xx, yy, c=zz, s=95, cmap="Blues", vmin=0, vmax=1, edgecolors="#496773", linewidths=.35)
    labels = []
    for row in clinical.itertuples():
        geo = str(row.geo_grade_label).replace("Grade ", "")
        paper = str(row.paper_pfirrmann_if_number_crosswalk_valid)
        labels.append(row.donor_id.replace("Tu_", "Tu ")+f"\n{geo} / {paper}*")
    ax.set_xticks(range(len(donors))); ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks(range(5)); ax.set_yticklabels(PROGRAMS)
    ax.set_xlim(-.6, len(donors)-.4); ax.set_ylim(4.65, -.55)
    ax.tick_params(length=0)
    ax.set_title("Fixed-program median scores in Tu NP-marker-supported cells", loc="left", pad=12)
    for spine in ax.spines.values(): spine.set_visible(False)
    cbar = fig.colorbar(dots, ax=ax, fraction=.025, pad=.03, aspect=22)
    cbar.set_label("Donor median " + tu_unit, fontsize=8)
    ax.set_xlabel("Donor; GEO grade / putative paper grade*", labelpad=7)
    ax.text(0, -.35, "* Sample-number crosswalk unresolved; no clinical grade-group test.",
            transform=ax.transAxes, va="top", fontsize=8)

    ax = fig.add_subplot(gs[1, 0]); letter(ax, "b")
    scatter = dotplot(ax, gene_dot, selected_genes, identities, "original_gene", orientation="genes_y")
    ax.set_title("Frozen GEP2: all 25 genes", loc="left", pad=10)
    cbar = fig.colorbar(scatter, ax=ax, fraction=.045, pad=.035, aspect=25)
    cbar.set_label("Mean donor log1p expression", fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    size_legend(ax, figure_anchor=(.14, .002))

    context_grid = gs[1, 1].subgridspec(2, 1, height_ratios=[.65, 1.35], hspace=.90)
    for cohort_index, cohort in enumerate(["Swahn_GSE230808", "Tu_GSE165722"]):
        ax = fig.add_subplot(context_grid[cohort_index, 0])
        if cohort_index == 0:
            letter(ax, "c")
        ordered = []
        names = ["NP_chondrocyte_like", "other_disc_stromal"] if cohort.startswith("Swahn") else IDENTITIES
        for identity in names:
            subset = context[context.cohort.eq(cohort) & context.identity.eq(identity)].sort_values("donor_id")
            subset = subset[np.isfinite(subset.score_median)]
            if not subset.empty:
                ordered.append(subset)
        labels = []
        for y, subset in enumerate(ordered):
            values = subset.score_median.to_numpy()
            n = len(values)
            offset = np.linspace(-.15, .15, n) if n > 1 else np.zeros(1)
            color = COLORS["Tu"] if cohort.startswith("Tu") else COLORS["Swahn"]
            ax.scatter(values, y+offset, c=color, s=21, edgecolors="white", linewidths=.35, zorder=3)
            ax.plot([values.mean(), values.mean()], [y-.25, y+.25], color="#273F48", linewidth=1.2, zorder=4)
            labels.append(wrap(IDENTITY_NAMES.get(subset.iloc[0].identity, subset.iloc[0].identity)+f" (n={n})", 21))
        ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(max(len(labels)-.5, .5), -.55)
        cohort_unit = tu_unit if cohort.startswith("Tu") else "UCell"
        # Independent axes explicitly prevent a between-method visual comparison.
        finite_scores = context.loc[context.cohort.eq(cohort), "score_median"].dropna()
        max_score = min(1.0, max(.1, float(finite_scores.max()) * 1.08)) if len(finite_scores) else 1.0
        ax.set_xlim(-.01, max_score)
        ax.set_xlabel("GEP2 donor median " + cohort_unit)
        ax.set_title(("Tu: " if cohort.startswith("Tu") else "Swahn: ") + cohort_unit,
                     loc="left", pad=9, fontsize=9)
        ax.grid(axis="x", color="#E5E8E9", linewidth=.6); ax.set_axisbelow(True)
    fig.text(.14, -.045, "Each context point is one donor; vertical mark = donor mean. AUCell and UCell scales are not compared.", fontsize=8)
    export(fig, output, "Figure2_Fixed_NP_Programs_CrossCohort")
    legend = ("Fig. 2 | Fixed NP programs assessed at donor level. "
              f"(a) Actual donor median scores for all five original 25-gene programs in Tu cells with {IDENTITY_NAMES.get(args.tu_primary_identity, args.tu_primary_identity)}. "
              "Each plotted circle represents a donor/program estimate; an unavailable donor/identity combination is marked by a grey cross. "
              "The first grade label is reported by GEO; the second is the putative paper Table 1 grade obtained by matching sample numbers only. The GSM-to-paper-patient crosswalk is unresolved, and no clinical grade association test is claimed. "
              "(b) The original frozen GEP2 Top25 list, unchanged in its original rank order. Point area denotes detection fraction averaged equally across contributing donors within each Tu identity; colour denotes equally donor-weighted mean log1p of published processed expression. Measured zero detection is an open dot; unavailable measurements are grey crosses. "
              "(c) GEP2 donor median scores within the actual Swahn discovery identities and Tu marker-support identities, on separate cohort-specific axes. Each point is one donor; the vertical mark is the mean of donor medians; n denotes contributing donors. "
              f"Swahn uses official discovery UCell scores and Tu uses {tu_unit} scores; the input table records whether each score set was reconstructed or reused. AUCell and UCell have different score definitions; their numerical scales are not compared across studies. The clinical question is assessed from within-cohort expression and identity context, with platform and processed-matrix differences retained. Marker support is not a definitive cell atlas classification. No cell-level disease-effect P values or positive-value filling are used.")
    (output / "Figure2_LEGEND.md").write_text(legend+"\n", encoding="utf-8")
    return {"Tu_primary_identity": args.tu_primary_identity, "Tu_primary_donors": primary.donor_id.nunique(),
            "frozen_GEP2_genes": len(selected_genes), "Tu_score_method": tu_methods[0],
            "context_axes": "separate cohort/method scales", "grade_inference_run": False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--af-donor-counts", type=Path)
    parser.add_argument("--tu-primary-identity", default="NP_matrix_marker_supported")
    args = parser.parse_args()
    results = args.results_dir or args.work / "01_human_programs"
    output = args.output_dir or args.work / "05_figures" / "Figure1_Figure2"
    output.mkdir(parents=True, exist_ok=True)
    configure()
    report = {"figure1": figure1(args, results, output), "figure2": figure2(args, results, output),
              "export": "PNG 300dpi and editable-text PDF/SVG; 180mm nominal width; font 8-10pt",
              "visual_QA": "PENDING_ACTUAL_RENDER_INSPECTION"}
    (output / "FIGURES12_RUN_STATUS.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
