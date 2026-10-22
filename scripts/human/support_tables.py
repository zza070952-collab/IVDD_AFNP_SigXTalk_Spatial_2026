#!/usr/bin/env python3
"""Build required Stage A donor/support tables and honest executed text."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--review-decisions", type=Path)
    args = parser.parse_args()
    out = args.work / "01_human_programs"
    score_path = out / "PROGRAM_SCORE_DONOR_SUMMARY.tsv"
    scores = pd.read_csv(score_path, sep="\t")
    coverage = pd.read_csv(out / "PROGRAM_GENE_COVERAGE.tsv", sep="\t")
    genes = pd.read_csv(out / "TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv", sep="\t")
    clinical = pd.read_csv(out / "COHORT_AND_DONOR_TABLE.tsv", sep="\t")
    sample_manifest = pd.read_csv(args.sample_manifest, sep="\t")
    review = {}
    review_status = "NO_REVIEW_DECISIONS_PROVIDED"
    if args.review_decisions:
        document = json.loads(args.review_decisions.read_text())
        if document.get("source_summary_sha256") == sha(score_path):
            review = document["program_decisions"]
            review_status = "DOCUMENTED_INTERPRETATION_BOUND_TO_EXACT_ACTUAL_SUMMARY"
        else:
            review_status = "REVIEW_REQUIRED_INPUT_HASH_CHANGED"
    support = []
    selected = scores[scores.program.eq("GEP2")]
    for cohort in sorted(selected.cohort.unique()):
        own = selected[selected.cohort.eq(cohort)]
        primary = "NP_matrix_marker_supported" if cohort.startswith("Tu") else "NP_chondrocyte_like"
        control = "immune_marker_supported" if cohort.startswith("Tu") else "other_disc_stromal"
        for donor in sorted(own.donor_id.unique()):
            donor_rows = own[own.donor_id.eq(donor)]
            np_row = donor_rows[donor_rows.identity.eq(primary)]
            control_row = donor_rows[donor_rows.identity.eq(control)]
            record = {"cohort": cohort, "donor_id": donor, "primary_identity": primary,
                      "control_identity": control, "grade_association": "UNTESTABLE_UNRESOLVED_TU_CROSSWALK" if cohort.startswith("Tu") else "NOT_A_CROSS_SERIES_DISEASE_TEST"}
            for label, frame in [("NP", np_row), ("context", control_row)]:
                record[label + "_cells"] = int(frame.n_cells.iloc[0]) if len(frame) else 0
                record[label + "_donor_median"] = float(frame.score_median.iloc[0]) if len(frame) else np.nan
            record["score_method"] = str(np_row.score_method.iloc[0]) if len(np_row) else str(donor_rows.score_method.iloc[0])
            record["NP_gt_context"] = record["NP_donor_median"] > record["context_donor_median"] if np.isfinite(record["NP_donor_median"]) and np.isfinite(record["context_donor_median"]) else None
            if cohort.startswith("Tu"):
                donor_genes = genes[genes.donor_id.eq(donor) & genes.identity.eq(primary) & genes.program.eq("GEP2")]
                for candidate in ["PRSS23", "PLOD2", "CCN2", "LOX", "LOXL2", "TAGLN", "MYL9"]:
                    row = donor_genes[donor_genes.original_gene.eq(candidate)]
                    record[candidate + "_detection_fraction"] = float(row.fraction_detected.iloc[0]) if len(row) else np.nan
                    record[candidate + "_released_symbol"] = str(row.mapped_gene.iloc[0]) if len(row) and pd.notna(row.mapped_gene.iloc[0]) else ""
            support.append(record)
    pd.DataFrame(support).to_csv(out / "GEP2_DONOR_SUPPORT.tsv", sep="\t", index=False)
    # Separate factual donor provenance from a putative grade crosswalk.
    manifests = []
    for row in clinical.to_dict("records"):
        manifests.append({"cohort": "Tu_GSE165722", "donor_id": row["donor_id"], "gsm": row["gsm"],
            "tissue": "NP", "platform": "BD_Rhapsody", "reported_grade": row["geo_grade_label"],
            "putative_paper_grade": row["paper_pfirrmann_if_number_crosswalk_valid"],
            "grade_mapping_status": row["crosswalk_status"], "submitted_cells": int(row["n_cells"]),
            "QC_cells": int(row["n_qc_pass"])})
    discovery_counts = scores[scores.cohort.eq("Swahn_GSE230808") & scores.program.eq("GEP2")].groupby(["donor_id", "gsm"]).n_cells.sum().reset_index()
    for row in discovery_counts.to_dict("records"):
        meta = sample_manifest[sample_manifest.gsm.eq(row["gsm"])]
        if len(meta) != 1 or meta.donor_id.iloc[0] != row["donor_id"]:
            raise ValueError("Discovery donor identity does not match official sample manifest")
        manifests.append({"cohort": "Swahn_GSE230808", "donor_id": row["donor_id"], "gsm": row["gsm"],
            "tissue": "NP", "platform": "10x_3prime_v3.1", "reported_grade": meta.thompson_grade.iloc[0],
            "putative_paper_grade": "not_applicable", "grade_mapping_status": "official_patient_id_and_disc",
            "submitted_cells": np.nan, "QC_cells": int(row["n_cells"])})
    pd.DataFrame(manifests).to_csv(out / "COHORT_AND_DONOR_MANIFEST.tsv", sep="\t", index=False)
    brief = []
    tu = scores[scores.cohort.eq("Tu_GSE165722")]
    for program in [f"GEP{i}" for i in range(1, 6)]:
        subset = tu[tu.program.eq(program)]
        own = subset[subset.identity.eq("NP_matrix_marker_supported")].score_median.dropna()
        cov = coverage[coverage.program.eq(program)]
        identity = subset.groupby("identity").score_median.mean().dropna()
        decision = review.get(program, {})
        brief.append({"program": program, "main_focus": program == "GEP2", "Tu_donors_with_NP_scores": len(own),
            "mapped_genes_min": int(cov.n_mapped_genes.min()), "mapped_genes_max": int(cov.n_mapped_genes.max()),
            "NP_donor_median_min": float(own.min()) if len(own) else np.nan,
            "NP_donor_median_max": float(own.max()) if len(own) else np.nan,
            "peak_mean_donor_score_identity": identity.idxmax() if len(identity) else "not_assessable",
            "expression_interpretation": decision.get("label", "REQUIRES_BIOLOGICAL_REVIEW"),
            "interpretation_scope": decision.get("scope", "Data-derived expression context; no arbitrary score threshold declares replication"),
            "grade_interpretation": "UNTESTABLE", "review_status": review_status})
    brief = pd.DataFrame(brief)
    brief.to_csv(out / "ALL_PROGRAMS_BRIEF.tsv", sep="\t", index=False)
    brief[~brief.main_focus].to_csv(out / "OTHER_PROGRAMS_BRIEF.tsv", sep="\t", index=False)
    algorithm = json.loads((out / "AUCELL_RUN_METADATA.json").read_text())
    lock = json.loads((args.frozen_dir / "PROGRAM_LOCK.json").read_text())
    scale = ["# Provenance and score scale", "",
        f"Frozen K={lock['K']}; original fitting seed={lock['seed']}; lock UTC={lock['locked_utc']}; lock SHA256={sha(args.frozen_dir / 'PROGRAM_LOCK.json')}.", "",
        "The original cNMF version was1.7.1; the preserved own-code raw QC/fitting and documented invocation retain the original scientific parameters. No cNMF recomputation is implicit in projection or table regeneration.", "",
        f"Tu implementation: {algorithm['implementation']}, version {algorithm['version']}; n_up={algorithm['n_up']}, tmin={algorithm['tmin']}. {algorithm['ties']}", "",
        algorithm["denominator"], "",
        "Discovery historical scoring used official R UCell2.4.0, maxRank1500, average ties, one core and200-cell chunks. Rebuilt identity preprocessing used seed20261016, whereas original cNMF fitting used seed20261012.", "",
        "AUCell and UCell scores differ in definition and feature handling. Do not directly compare their numerical scales or describe a score difference across cohorts as a disease effect. Plot separate cohort/method axes and use gene coverage, within-donor context and measured expression.", "",
        "The Tu GEO I-IV versus paper II-V mapping remains unconfirmed. No clinical grade-group test or disease-grade replication is performed. Gene-ID aliases require authoritative HGNC identity plus exactly one matching released feature; absent or ambiguous mappings remain unassessable.", "",
        "Cell-level metadata and scores are private runtime inputs to this repository; obtain public matrices from GEO and recreate them using the scripts. They are not redistributed as public code data."]
    (out / "PROVENANCE_AND_SCORE_SCALE.md").write_text("\n".join(scale)+"\n", encoding="utf-8")
    np_counts = scores[scores.program.eq("GEP2") & scores.identity.isin(["NP_matrix_marker_supported", "NP_chondrocyte_like"])].groupby("cohort").n_cells.sum().to_dict()
    methods = ["# Human programs: actually generated methods and results", "",
        "## Methods", "",
        "Frozen cNMF programs were projected without refitting. Tu processed gene-by-cell expression was parsed into sparse matrices; QC, transparent marker-support labels, authoritative unique gene-ID mappings and actual official AUCell output are preserved in the input/output manifests. Discovery scores and annotations carry their actual origin labels in PROGRAM_SCORE_DONOR_SUMMARY.tsv. Donor medians are the descriptive replicate summaries; cells are not independent patients.", "",
        "## Results", "",
        f"The regenerated donor manifest contains {clinical.donor_id.nunique()} Tu donors and {discovery_counts.donor_id.nunique()} discovery donors. Tu retained QC cells total {int(clinical.n_qc_pass.sum()):,}; NP-identity counts by cohort are {np_counts}.", "",
        "| Program | Mapped genes | Tu NP-donor median range | Expression review |", "|---|---|---|---|"]
    for row in brief.itertuples():
        methods.append(f"|{row.program}|{row.mapped_genes_min}-{row.mapped_genes_max}/25|{row.NP_donor_median_min:.4f}-{row.NP_donor_median_max:.4f}|{row.expression_interpretation}|")
    methods += ["", f"Review status: {review_status}. If the actual-summary hash differs from a prior review record, the script leaves biological review required rather than inventing a new positive conclusion.", "",
        "Grade association remains UNTESTABLE. The program-associated expression context does not independently establish disease effects, an AF-to-NP molecular causal chain, or activation of a specific receptor by a particular ligand. Mixed/weak programs and technical coverage limits remain reported."]
    (out / "STAGE_A_METHODS_RESULTS.md").write_text("\n".join(methods)+"\n", encoding="utf-8")
    print(json.dumps({"donor_rows": len(manifests), "programs": len(brief), "review_status": review_status}), flush=True)


if __name__ == "__main__":
    main()
