#!/usr/bin/env python3
"""Donor-level summaries of actual Tu AUCell and official discovery UCell scores."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROGRAMS = ["GEP" + str(i) for i in range(1, 6)]


def truth(series):
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def summarize_cells(metadata, scores, cohort, mask, score_method, suffix, origin_override=None):
    if metadata.cell_id.duplicated().any() or scores.cell_id.duplicated().any():
        raise ValueError("Duplicate cell identifiers")
    frame = metadata.loc[mask].merge(scores, on="cell_id", validate="one_to_one")
    if len(frame) != int(mask.sum()):
        raise ValueError("Missing cell scores after one-to-one join")
    columns = [p + "_top25_" + suffix for p in PROGRAMS]
    if np.isinf(frame[columns].to_numpy()).any():
        raise ValueError("Infinite actual scores")
    frame["cohort"] = cohort
    for field in ["geo_grade_label", "paper_pfirrmann_if_number_crosswalk_valid", "crosswalk_status"]:
        if field not in frame:
            frame[field] = "not_applicable"
    groups = ["cohort", "donor_id", "gsm", "identity", "geo_grade_label",
              "paper_pfirrmann_if_number_crosswalk_valid", "crosswalk_status"]
    records = []
    for key, group in frame.groupby(groups, dropna=False, sort=True):
        base = dict(zip(groups, key))
        for program in PROGRAMS:
            values = group[program + "_top25_" + suffix].dropna().to_numpy()
            records.append({**base, "program": program, "n_cells": len(group), "n_valid_scores": len(values),
                            "score_mean": float(values.mean()) if values.size else np.nan,
                            "score_median": float(np.median(values)) if values.size else np.nan,
                            "score_q25": float(np.quantile(values, .25)) if values.size else np.nan,
                            "score_q75": float(np.quantile(values, .75)) if values.size else np.nan,
                            "fraction_score_gt_zero": float((values > 0).mean()) if values.size else np.nan,
                            "biological_replicate": "donor", "score_method": score_method,
                            "score_origin": origin_override or ("official_new_AUCell" if cohort.startswith("Tu") and suffix == "AUCell" else "official_new_UCell" if cohort.startswith("Tu") else "verified_legacy_UCell")})
    return pd.DataFrame(records), frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--discovery-metadata", type=Path, required=True)
    parser.add_argument("--discovery-scores", type=Path, required=True)
    parser.add_argument("--tu-score-method", choices=["AUCell", "UCell"], default="AUCell")
    parser.add_argument("--discovery-score-origin", choices=["verified_legacy_UCell", "reconstructed_official_UCell"], default="verified_legacy_UCell")
    args = parser.parse_args()
    out = args.work / "01_human_programs"
    tu_metadata = pd.read_csv(out / "TU_CELL_METADATA.tsv.gz", sep="\t")
    tu_score_file = out / ("TU_PROGRAM_AUCELL_CELL_SCORES.tsv.gz" if args.tu_score_method == "AUCell" else "TU_PROGRAM_UCELL_CELL_SCORES.tsv.gz")
    tu_scores = pd.read_csv(tu_score_file, sep="\t")
    tu_method = "AUCell_decoupler_2.2.0" if args.tu_score_method == "AUCell" else "UCell_R_official"
    tu_summary, tu_cells = summarize_cells(tu_metadata, tu_scores, "Tu_GSE165722", truth(tu_metadata.qc_pass), tu_method, args.tu_score_method)
    discovery_metadata = pd.read_csv(args.discovery_metadata, sep="\t")
    discovery_scores = pd.read_csv(args.discovery_scores, sep="\t")
    mask = truth(discovery_metadata.pass_qc)
    discovery_summary, discovery_cells = summarize_cells(discovery_metadata, discovery_scores, "Swahn_GSE230808", mask, "UCell_2.4.0", "UCell", args.discovery_score_origin)
    if discovery_metadata.donor_id.nunique() != 8 or set(discovery_metadata.geo_series) != {"GSE230808"}:
        raise ValueError("Unexpected discovery cohort/donor identity")
    combined = pd.concat([discovery_summary, tu_summary], ignore_index=True)
    combined.to_csv(out / "PROGRAM_SCORE_DONOR_SUMMARY.tsv", sep="\t", index=False)
    combined.to_csv(out / "FIGURE2_DONOR_PROGRAM_SOURCE.tsv", sep="\t", index=False)
    counts = pd.concat([
        tu_cells.groupby(["cohort", "donor_id", "gsm", "identity"]).size().rename("n_cells").reset_index(),
        discovery_cells.groupby(["cohort", "donor_id", "gsm", "identity"]).size().rename("n_cells").reset_index()
    ], ignore_index=True)
    counts.to_csv(out / "FIGURE1_COHORT_DONOR_IDENTITY_COUNTS.tsv", sep="\t", index=False)
    coverage = pd.read_csv(out / "PROGRAM_GENE_COVERAGE.tsv", sep="\t")
    notes = []
    for program in PROGRAMS:
        subset = tu_summary[tu_summary.program == program]
        cov = coverage[coverage.program == program]
        identity_means = subset.groupby("identity").score_median.mean().dropna()
        peak = identity_means.idxmax() if len(identity_means) else "no_cells"
        np_only = subset[subset.identity == "NP_matrix_marker_supported"]
        donor_medians = np_only.score_median.dropna().to_numpy()
        notes.append({"program": program, "Tu_donors_with_qc_cells": subset.donor_id.nunique(),
                      "Tu_donors_with_NP_matrix_supported_cells": np_only.donor_id.nunique(),
                      "Tu_donors_with_valid_NP_program_scores": np_only.loc[np_only.score_median.notna(), "donor_id"].nunique(),
                      "mapped_genes_min": int(cov.n_mapped_genes.min()),
                      "mapped_genes_max": int(cov.n_mapped_genes.max()),
                      "NP_donor_median_score_min": float(donor_medians.min()) if donor_medians.size else np.nan,
                      "NP_donor_median_score_max": float(donor_medians.max()) if donor_medians.size else np.nan,
                      "identity_with_highest_mean_donor_median": peak,
                      "grade_association_status": "UNTESTABLE_UNRESOLVED_GSM_PATIENT_CROSSWALK",
                      "expression_support_status": "LOW_TECHNICAL_COVERAGE" if cov.n_mapped_genes.min() < 10 else "AVAILABLE_FOR_BIOLOGICAL_REVIEW"})
    conclusion = pd.DataFrame(notes)
    conclusion.to_csv(out / "PROGRAM_EXPRESSION_CONTEXT_SUMMARY.tsv", sep="\t", index=False)
    lines = ["# Cross-cohort fixed-program expression assessment", "",
             f"The five original 25-gene programs were kept fixed. Tu uses actual {tu_method} scores; Swahn score provenance is {args.discovery_score_origin}. No cNMF refitting is part of this summary step.", "",
             f"Actual analysed Tu QC cells: {len(tu_cells):,} across {tu_cells.donor_id.nunique()} donors. Actual Swahn discovery stromal cells: {len(discovery_cells):,} across {discovery_cells.donor_id.nunique()} donors; score origin: {args.discovery_score_origin}.", "",
             "These separate cohorts/platforms and marker-support labels are not an integrated atlas. AUCell and UCell use different score definitions and their numerical scales are not compared. Interpret expression/identity context within each cohort; the legacy UCell version/normalisation must remain disclosed where available.", "",
             "## Grade boundary", "",
             "GEO reports grades I–IV whereas paper Table 1 reports Pfirrmann II–V. Matching SampleN and SN is not independently confirmed by age, sex, diagnosis or an author-provided accession crosswalk. Both labels remain in the donor table. No clinical grade-group P values, healthy-versus-disease inference, or significant replication claim was generated.", "",
             "## Actual expression context", "",
             "| Program | Tu donors with NP marker support | Mapped genes | NP donor-median score range | Highest donor-weighted identity |", "|---|---:|---|---|---|"]
    for row in notes:
        bounds = f"{row['NP_donor_median_score_min']:.4g}–{row['NP_donor_median_score_max']:.4g}" if np.isfinite(row['NP_donor_median_score_min']) else "not assessable"
        lines.append(f"| {row['program']} | {row['Tu_donors_with_NP_matrix_supported_cells']} | {row['mapped_genes_min']}–{row['mapped_genes_max']}/25 | {bounds} | {row['identity_with_highest_mean_donor_median']} |")
    lines += ["", "Expression coverage and cell identity inform CONSISTENT/MIXED/UNSUPPORTED interpretation. Final biological labels require review of the actual marker and program-gene dotplots; no arbitrary score threshold declares replication. Grade association remains UNTESTABLE. GEP2 is the pre-specified focus, not a required winner."]
    (out / "PROGRAM_CROSSCOHORT_CONCLUSION.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (out / "STAGE_A_COMPUTATION_STATUS.json").write_text(json.dumps({"status": "DONE_EXPRESSION_ONLY_GRADE_UNTESTABLE",
        "Tu_donors": int(tu_cells.donor_id.nunique()), "Tu_cells_qc_pass": len(tu_cells),
        "Swahn_donors": int(discovery_cells.donor_id.nunique()), "Swahn_cells_reused": len(discovery_cells),
        "Tu_score_method": tu_method, "Swahn_score_method": "UCell_2.4.0", "Swahn_score_origin": args.discovery_score_origin,
        "cross_method_numeric_scale_comparison": False,
        "clinical_grade_test_run": False, "biological_interpretation": "requires review of actual outputs"}, indent=2))
    print(json.dumps({"summaries_written": len(combined), "Tu_donors": int(tu_cells.donor_id.nunique()),
                      "Tu_qc_cells": len(tu_cells), "discovery_cells": len(discovery_cells)}), flush=True)


if __name__ == "__main__":
    main()
