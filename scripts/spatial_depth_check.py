#!/usr/bin/env python3
"""One preregistered molecule-count depth sensitivity of frozen spatial cores.

Reads prior sparse counts/coordinates/scores without modifying them. Reproduces
the original GEP2 score first, then samples exactly 100 molecules without
replacement once per eligible bin. No threshold search or bin-level P value.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.metadata as md
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024**2), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf8")


def truth(series):
    return series.astype(str).str.lower().isin(["true", "1"])


def normalized_score(x, total, columns):
    z = x[:, columns].astype(float).multiply((10000 / np.maximum(total, 1))[:, None]).tocsr()
    z.data = np.log1p(z.data)
    return np.asarray(z.mean(axis=1)).ravel()


def stable_rng(seed, capture, spot_id):
    material = f"{seed}\0{capture}\0{spot_id}".encode("utf8")
    integer = int.from_bytes(hashlib.sha256(material).digest()[:16], "little")
    return np.random.default_rng(integer)


def sample_counts(counts, target, rng):
    """Exact multivariate hypergeometric, retaining no more than true counts."""
    counts = np.asarray(counts, dtype=np.int64)
    population = int(counts.sum())
    if population < target:
        raise ValueError("Cannot rarefy a library below the fixed target")
    if population == target:
        return counts.copy()
    if hasattr(rng, "multivariate_hypergeometric"):
        result = rng.multivariate_hypergeometric(counts, target, method="marginals")
    else:
        # Same exact distribution for older NumPy; no multinomial replacement.
        result = np.zeros_like(counts)
        remaining_draws = target
        for index, count in enumerate(counts[:-1]):
            if remaining_draws == 0:
                break
            result[index] = rng.hypergeometric(int(count), population - int(count), remaining_draws)
            remaining_draws -= int(result[index])
            population -= int(count)
        result[-1] = remaining_draws
    assert int(result.sum()) == target and (result >= 0).all() and (result <= counts).all()
    return result


def rho_row(capture, scope, score, total):
    valid = np.isfinite(score) & np.isfinite(total)
    s, t = np.asarray(score)[valid], np.asarray(total)[valid]
    constant = len(s) < 3 or np.unique(s).size < 2 or np.unique(t).size < 2
    rho = np.nan if constant else float(spearmanr(s, t).correlation)
    return {"capture_id": capture, "scope": scope, "n_bins_technical_only": int(valid.sum()),
            "spearman_rho_score_vs_UMI": rho, "P_value": "NOT_COMPUTED_BINS_NOT_BIOLOGICAL_N",
            "interpretation": "descriptive_depth_association" if not constant else "UNTESTABLE_CONSTANT_OR_TOO_FEW_VALUES"}


def expected_roi_rows(rois, roi_sha):
    rows = []
    for capture in rois["captures"]:
        for disc in capture["discs"]:
            for region in disc["regions"]:
                rows.append({"capture_id": capture["capture_id"], "source_GSM": capture["source_GSM"],
                             "genotype": capture["genotype"], "disc_id": disc["disc_id"],
                             "histology_region": region["region"], "flank": region["side"],
                             "vertices": json.dumps(region["vertices"], separators=(",", ":")),
                             "frozen_roi_sha256": roi_sha, "source_image_sha256": capture["source_image_sha256"],
                             "microns_per_pixel": capture["microns_per_pixel"],
                             "edge_margin_pixels": np.sqrt(2) * 48 / (2 * capture["microns_per_pixel"]),
                             "annotation_edited": False, "expert_histology_verification": False})
    return pd.DataFrame(rows)


def process_capture(root, capture, old_source, members, config, output):
    capture_id = capture["capture_id"]
    folder = root / "prepared48" / capture_id
    x = sparse.load_npz(folder / "counts48_symbols_qc.npz").tocsr()
    obs = pd.read_csv(folder / "BIN48_METADATA_QC.tsv.gz", sep="\t")
    var = pd.read_csv(folder / "MOUSE_SYMBOL_FEATURES.tsv", sep="\t")
    if x.shape != (len(obs), len(var)) or not obs.spot_id.is_unique or not var.mouse_gene_symbol.is_unique:
        raise ValueError("Sparse counts/metadata/feature identity mismatch")
    if not (np.isfinite(x.data).all() and (x.data >= 0).all() and (x.data == np.floor(x.data)).all()):
        raise ValueError("Depth sampling requires real nonnegative integer molecule counts")
    x = x.astype(np.int64)
    x.eliminate_zeros()
    old = old_source.loc[old_source.capture_id.eq(capture_id)].set_index("spot_id")
    if set(old.index) != set(obs.spot_id) or not old.index.is_unique:
        raise ValueError("Prior QC-score universe differs from frozen count metadata")
    old = old.loc[obs.spot_id].reset_index()
    if not np.allclose(old[["x_fullres", "y_fullres"]], obs[["x_fullres", "y_fullres"]], rtol=0, atol=1e-9):
        raise ValueError("Original score coordinates do not reproduce unchanged count geometry")
    total = np.asarray(x.sum(axis=1)).ravel().astype(np.int64)
    if not np.array_equal(total, old.n_counts.to_numpy(dtype=np.int64)):
        raise ValueError("Original UMI library sizes differ from the frozen score source")
    chosen = members.loc[members.capture_id.eq(capture_id) & members.program.eq(config["program"])].copy()
    chosen = chosen.loc[truth(chosen.strict_one_to_one_native_feature)]
    if chosen.human_gene.duplicated().any() or chosen.mouse_symbol.duplicated().any():
        raise ValueError("Fixed score membership contains duplicate human/mouse genes")
    index = dict(zip(var.mouse_gene_symbol, range(len(var))))
    columns = [index[symbol] for symbol in chosen.mouse_symbol]
    if not columns:
        raise ValueError("No frozen mapped program members available")
    score = normalized_score(x, total, columns)
    score_name = config["program"] + "_mean_log1p_CP10k"
    maximum_error = float(np.max(np.abs(score - old[score_name].to_numpy(dtype=float))))
    if maximum_error > config["score_reproduction_absolute_tolerance"]:
        raise ValueError(f"Original score reproduction failed: max error {maximum_error}")
    core = old.anatomical_region.isin(config["regions"]).to_numpy()
    relations = [rho_row(capture_id, "ALL_ORIGINAL_QC_BINS", score, total),
                 rho_row(capture_id, "FROZEN_NP_AF_CORES_COMBINED", score[core], total[core])]
    for region in config["regions"]:
        take = old.anatomical_region.eq(region).to_numpy()
        relations.append(rho_row(capture_id, region, score[take], total[take]))
    target = config["target_total_UMI"]
    eligible = core & (total >= config["minimum_original_total_UMI"])
    program_mask = np.zeros(x.shape[1], dtype=bool)
    program_mask[columns] = True
    values, thin_rows, thin_indices, indptr = [], [], [], [0]
    for i in np.flatnonzero(core):
        row = old.iloc[i]
        result = {"capture_id": capture_id, "source_GSM": capture["source_GSM"], "genotype": capture["genotype"],
                  "disc_id": row.disc_id, "region": row.anatomical_region, "spot_id": row.spot_id,
                  "x_fullres": float(row.x_fullres), "y_fullres": float(row.y_fullres),
                  "original_total_UMI": int(total[i]), "original_program_score": float(score[i]),
                  "eligible_at_fixed_depth": bool(eligible[i]), "target_total_UMI": target,
                  "sampled_total_UMI": np.nan, "sampled_program_score": np.nan,
                  "sampled_detected_genes": np.nan, "sampled_detected_program_members": np.nan,
                  "n_fixed_available_program_members": len(columns)}
        if eligible[i]:
            start, end = x.indptr[i:i+2]
            counts, genes = x.data[start:end], x.indices[start:end]
            thinned = sample_counts(counts, target, stable_rng(config["seed"], capture_id, row.spot_id))
            program_counts = thinned[program_mask[genes]]
            result.update(sampled_total_UMI=int(thinned.sum()),
                          sampled_program_score=float(np.log1p(program_counts * (10000 / target)).sum() / len(columns)),
                          sampled_detected_genes=int(np.count_nonzero(thinned)),
                          sampled_detected_program_members=int(np.count_nonzero(program_counts)))
            positive = thinned > 0
            thin_rows.extend(thinned[positive].tolist())
            thin_indices.extend(genes[positive].tolist())
            indptr.append(len(thin_rows))
        values.append(result)
    table = pd.DataFrame(values)
    sampled = sparse.csr_matrix((np.asarray(thin_rows, dtype=np.int64), np.asarray(thin_indices, dtype=np.int64),
                                 np.asarray(indptr, dtype=np.int64)), shape=(int(eligible.sum()), x.shape[1]))
    assert np.all(np.asarray(sampled.sum(axis=1)).ravel() == target)
    private = output / "private_depth_count_audit"
    private.mkdir(exist_ok=True)
    sparse.save_npz(private / f"{capture_id}_single_depth{target}_counts.npz", sampled)
    table.loc[table.eligible_at_fixed_depth, ["spot_id"]].to_csv(private / f"{capture_id}_sampled_matrix_row_ids.tsv", sep="\t", index=False)
    proof = {"capture_id": capture_id, "original_QC_bins": len(obs), "original_core_bins": int(core.sum()),
             "original_program_reproduction_max_abs_error": maximum_error,
             "n_frozen_available_members": len(columns), "eligible_core_bins": int(eligible.sum()),
             "sampled_total_UMI_per_bin": target, "sampled_matrix_shape": list(sampled.shape),
             "original_coordinates_verified_unchanged": True}
    return table, relations, proof


def summarize_disc(values, minimum):
    rows = []
    for keys, group in values.groupby(["capture_id", "source_GSM", "genotype", "disc_id", "region"], sort=True):
        eligible = group.loc[group.eligible_at_fixed_depth]
        n = len(eligible)
        row = dict(zip(["capture_id", "source_GSM", "genotype", "disc_id", "region"], keys))
        row.update(original_n_bins=len(group), retained_n_bins=n, retention_fraction=n / len(group),
                   original_UMI_min=int(group.original_total_UMI.min()), original_UMI_max=int(group.original_total_UMI.max()),
                   original_UMI_median=float(group.original_total_UMI.median()),
                   eligible_original_UMI_min=float(eligible.original_total_UMI.min()) if n else np.nan,
                   eligible_original_UMI_max=float(eligible.original_total_UMI.max()) if n else np.nan,
                   eligible_original_UMI_median=float(eligible.original_total_UMI.median()) if n else np.nan,
                   original_all_core_mean_score=float(group.original_program_score.mean()),
                   original_eligible_mean_score=float(eligible.original_program_score.mean()) if n else np.nan,
                   sampled_depth100_mean_score=float(eligible.sampled_program_score.mean()) if n else np.nan,
                   sampled_UMI_min=float(eligible.sampled_total_UMI.min()) if n else np.nan,
                   sampled_UMI_max=float(eligible.sampled_total_UMI.max()) if n else np.nan,
                   depth_coverage_status="ADEQUATE_DECLARED_DESCRIPTIVE_COVERAGE" if n >= minimum else "LOW_COVERAGE_RETAIN_IN_REPORT",
                   biological_unit="disc_nested_in_capture; bins_not_biological_replicates")
        rows.append(row)
    return pd.DataFrame(rows)


def make_pairs(disc, minimum):
    output = []
    for keys, group in disc.groupby(["capture_id", "genotype", "disc_id"], sort=True):
        row = dict(zip(["capture_id", "genotype", "disc_id"], keys))
        by_region = group.set_index("region")
        if not {"NP_CORE", "AF_CORE"}.issubset(by_region.index):
            raise ValueError("Frozen disc missing a compartment; do not create a fake pair")
        adequate = bool((by_region.retained_n_bins >= minimum).all())
        row["adequate_fixed_depth_pair"] = adequate
        for prefix, field in [("original_all", "original_all_core_mean_score"),
                              ("original_eligible", "original_eligible_mean_score"),
                              ("depth100", "sampled_depth100_mean_score")]:
            for region in ["NP_CORE", "AF_CORE"]:
                row[f"{prefix}_{region}"] = float(by_region.loc[region, field])
            row[f"{prefix}_NP_minus_AF"] = row[f"{prefix}_NP_CORE"] - row[f"{prefix}_AF_CORE"]
        row.update(NP_retained_bins=int(by_region.loc["NP_CORE", "retained_n_bins"]),
                   AF_retained_bins=int(by_region.loc["AF_CORE", "retained_n_bins"]))
        output.append(row)
    return pd.DataFrame(output)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spatial-input", required=True, type=Path, help="Read-only completed V2-compatible spatial stage")
    ap.add_argument("--roi-json", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    root, out = a.spatial_input.resolve(), a.out.resolve()
    if out == root or out.is_relative_to(root):
        raise ValueError("New analysis output must not be inside the old spatial project")
    out.mkdir(parents=True, exist_ok=True)
    config = json.loads(a.config.read_text(encoding="utf8"))
    if config["target_total_UMI"] != 100 or config["minimum_original_total_UMI"] != 100:
        raise ValueError("This predeclared analysis has one fixed depth of 100; no threshold scan")
    roi_sha = sha(a.roi_json)
    if roi_sha != config["expected_roi_sha256"]:
        raise ValueError("Frozen ROI hash changed; refuse new anatomical selection")
    rois = json.loads(a.roi_json.read_text(encoding="utf8"))
    source_path = root / "SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz"
    program_path = root / "NP_PROGRAM_BY_REGION.tsv"
    candidates_path = root / "CANDIDATE_GENES_BY_REGION.tsv"
    member_path = root / "PROGRAM_MEMBERS_USED_IN_SPATIAL_SCORE.tsv"
    roi_metadata_path = root / "roi_validation/DISC_AF_NP_ROI_METADATA.tsv.gz"
    protected = [a.roi_json, source_path, program_path, candidates_path, member_path, roi_metadata_path]
    for capture in rois["captures"]:
        folder = root / "prepared48" / capture["capture_id"]
        protected.extend(folder / name for name in ["counts48_symbols_qc.npz", "BIN48_METADATA_QC.tsv.gz",
                                                    "MOUSE_SYMBOL_FEATURES.tsv", "ALL_GRID48_GEOMETRY.tsv.gz",
                                                    "ORIGINAL_SCALEFACTORS.json"])
        image = folder / "source_images" / capture["source_image_filename"]
        if sha(image) != capture["source_image_sha256"]:
            raise ValueError("Original frozen H&E image changed")
        protected.append(image)
    before = {str(path.resolve()): sha(path) for path in protected}
    dump(out / "DEPTH_ANALYSIS_STATUS.json", {"status": "RUNNING", "started_utc": utc(),
                                             "config_sha256": sha(a.config), "ROI_sha256": roi_sha})
    old_source = pd.read_csv(source_path, sep="\t")
    if not old_source.spot_id.is_unique:
        raise ValueError("Duplicate frozen spatial source IDs")
    roi_metadata = pd.read_csv(roi_metadata_path, sep="\t").set_index("spot_id")
    aligned_roi = roi_metadata.loc[old_source.spot_id]
    if not np.array_equal(aligned_roi.disc_id.to_numpy(), old_source.disc_id.to_numpy()) or not np.array_equal(
            aligned_roi.anatomical_region.to_numpy(), old_source.anatomical_region.to_numpy()):
        raise ValueError("Original score labels differ from frozen post-validation ROI metadata")
    if not np.allclose(aligned_roi[["x_fullres", "y_fullres"]], old_source[["x_fullres", "y_fullres"]], rtol=0, atol=1e-9):
        raise ValueError("Original post-validation ROI coordinates differ from score geometry")
    members = pd.read_csv(member_path, sep="\t", keep_default_na=False)
    values, relations, proofs = [], [], []
    for capture in rois["captures"]:
        v, r, p = process_capture(root, capture, old_source, members, config, out)
        values.append(v)
        relations.extend(r)
        proofs.append(p)
    values = pd.concat(values, ignore_index=True)
    disc = summarize_disc(values, config["minimum_retained_bins_per_disc_region"])
    old_programs = pd.read_csv(program_path, sep="\t")
    old_programs = old_programs.loc[old_programs.program.eq(config["program"])]
    check = disc.merge(old_programs, left_on=["capture_id", "disc_id", "region"],
                       right_on=["capture_id", "disc_id", "region"], validate="one_to_one", suffixes=("", "_old"))
    if len(check) != len(disc) or not np.allclose(check.original_all_core_mean_score, check.mean_program_score,
                                              rtol=0, atol=config["score_reproduction_absolute_tolerance"]):
        raise ValueError("Frozen per-disc regional scores failed reproduction")
    pairs = make_pairs(disc, config["minimum_retained_bins_per_disc_region"])
    retained = values.loc[values.eligible_at_fixed_depth]
    retention = values.groupby(["capture_id", "source_GSM", "genotype", "region"], sort=True).agg(
        original_n_bins=("spot_id", "size"), retained_n_bins=("eligible_at_fixed_depth", "sum"),
        original_UMI_min=("original_total_UMI", "min"), original_UMI_max=("original_total_UMI", "max"),
        original_UMI_median=("original_total_UMI", "median")).reset_index()
    erange = retained.groupby(["capture_id", "region"]).agg(
        eligible_original_UMI_min=("original_total_UMI", "min"), eligible_original_UMI_max=("original_total_UMI", "max"),
        sampled_UMI_min=("sampled_total_UMI", "min"), sampled_UMI_max=("sampled_total_UMI", "max")).reset_index()
    retention = retention.merge(erange, how="left", on=["capture_id", "region"], validate="one_to_one")
    retention["retention_fraction"] = retention.retained_n_bins / retention.original_n_bins
    adequacy = pairs.groupby("capture_id").adequate_fixed_depth_pair.sum()
    adequate_pairs = pairs.loc[pairs.adequate_fixed_depth_pair]
    sufficient = bool((adequacy >= config["minimum_paired_discs_per_capture"]).all())
    positive = int((adequate_pairs.depth100_NP_minus_AF > 0).sum())
    n_adequate = len(adequate_pairs)
    before_delta = float(adequate_pairs.original_eligible_NP_minus_AF.mean()) if n_adequate else None
    after_delta = float(adequate_pairs.depth100_NP_minus_AF.mean()) if n_adequate else None
    if not sufficient:
        outcome = "DEPTH_SENSITIVITY_INCONCLUSIVE"
    elif positive == n_adequate:
        outcome = "DESCRIPTIVE_DIRECTION_PERSISTS_SINGLE_DEPTH_CHECK"
    else:
        outcome = "DEPTH_CONTROLLED_DIRECTION_MIXED_OR_REVERSED"
    attenuation = bool(before_delta is not None and after_delta is not None and after_delta < before_delta)
    for name, table in [
        ("DEPTH_SENSITIVITY_VALUES.tsv", values), ("DISC_LEVEL_PROGRAM_AND_DEPTH.tsv", disc),
        ("DEPTH_DISC_PAIR_COMPARISON.tsv", pairs), ("DEPTH_RETENTION_BY_CAPTURE_REGION.tsv", retention),
        ("GEP2_SCORE_UMI_ASSOCIATION.tsv", pd.DataFrame(relations)),
        ("ROI_FIXED_PROVENANCE.tsv", expected_roi_rows(rois, roi_sha)),
    ]:
        table.to_csv(out / name, sep="\t", index=False, na_rep="NA")
    candidates = pd.read_csv(candidates_path, sep="\t")
    candidates.to_csv(out / "CANDIDATE_GENES_PER_DISC_REGION.tsv", sep="\t", index=False)
    after = {str(path.resolve()): sha(path) for path in protected}
    if after != before:
        raise ValueError("Protected old input changed during analysis")
    pd.DataFrame([{"path": path, "sha256_before": digest, "sha256_after": after[path], "unchanged": True}
                  for path, digest in before.items()]).to_csv(out / "PRIVATE_FROZEN_INPUT_HASH_AUDIT.tsv", sep="\t", index=False)
    versions = {"python": platform.python_version(), **{p: md.version(p) for p in ["numpy", "scipy", "pandas"]}}
    status = {"status": "DONE_SINGLE_PREDECLARED_DEPTH_CHECK", "outcome": outcome, "finished_utc": utc(),
              "config_sha256": sha(a.config), "ROI_sha256": roi_sha, "target_total_UMI": 100,
              "random_seed": config["seed"], "single_draw_only": True, "versions": versions,
              "capture_proofs": proofs, "adequate_paired_discs_by_capture": {k: int(v) for k, v in adequacy.items()},
              "adequate_paired_discs": n_adequate, "NP_greater_AF_after_depth100": positive,
              "mean_paired_delta_original_same_eligible_bins": before_delta, "mean_paired_delta_depth100": after_delta,
              "observed_mean_difference_attenuated": attenuation, "protected_old_inputs_unchanged": True,
              "bin_level_P_values": False, "genotype_inference": False, "threshold_tuning": False}
    dump(out / "DEPTH_ANALYSIS_STATUS.json", status)
    retention_text = retention.to_csv(sep="\t", index=False, na_rep="NA")
    relation_text = pd.DataFrame(relations).to_csv(sep="\t", index=False, na_rep="NA")
    report = f"""# Fixed 100-UMI depth sensitivity

Status: **{outcome}**. One declared target of 100 molecules per original bin with at least 100 total UMIs was used, seed {config['seed']}. Sampling was exact multivariate hypergeometric without replacement across all detected genes. This is one single stochastic technical sensitivity check, not a threshold search, an animal-level test or an independent validation cohort.

The frozen ROI hash is `{roi_sha}`. Original per-bin GEP2 scores, total counts, coordinates and regional means reproduced within the declared tolerance before sampling. Protected original inputs were verified unchanged afterward. The same frozen available ortholog members and zero-retaining mean log1p CP10k definition were used after sampling; only molecule depth changed for retained bins.

## Retention and comparable depth

```tsv
{retention_text}```

All retained libraries contain exactly 100 molecules. The original eligible UMI ranges are reported rather than silently trimmed to an extra matching range. Eligibility changes the evaluated subset, so original-all-core, original-same-eligible and sampled-100 summaries are separate. Every excluded or low-coverage disc remains in DISC_LEVEL_PROGRAM_AND_DEPTH.tsv. A declared descriptive safeguard requires at least {config['minimum_retained_bins_per_disc_region']} retained bins in each compartment and {config['minimum_paired_discs_per_capture']} paired discs per capture. Actual adequate pairs: {dict(status['adequate_paired_discs_by_capture'])}.

## Original score versus library size, quantified before sampling

```tsv
{relation_text}```

These Spearman coefficients describe a technical bin-level relationship; no bin-level P value is computed and bins are not counted as biological replicates.

## Disc-context result

Among {n_adequate} adequately retained disc pairs, {positive} had NP>AF after the one fixed-depth draw. Mean NP-minus-AF contrast in the same eligible original bins was {before_delta}; after sampling it was {after_delta}. The observed mean difference was lower after sampling: {attenuation}. These are descriptive values, with no requirement that the direction remain positive and no attempt to optimize it. If coverage is insufficient, the outcome remains INCONCLUSIVE regardless of the numeric directions.

The finite single-draw variation, small cores and selective retention limit precision. If the contrast weakens or reverses, use anatomical expression context rather than corrected NP specificity. One capture per genotype makes capture and genotype inseparable; eight disc structures do not supply replicated capture arrays. This check does not prove cell adjacency, signaling direction or a molecular causal chain.
"""
    (out / "DEPTH_SENSITIVITY_REPORT.md").write_text(report, encoding="utf8")
    interpretation = f"""# Revised spatial interpretation

The frozen adult mouse anatomy and original program/candidate values were retained. The new predeclared 100-UMI check returned **{outcome}**, with {positive}/{n_adequate} adequately retained pairs showing NP>AF and an original-same-subset versus sampled mean contrast of {before_delta} versus {after_delta}. Coverage qualification is listed per disc and capture; this is not a significance test.

The original eight-disc NP-favored GEP2 pattern is descriptive anatomical expression context and was not independently depth-controlled at V2. Interpret its robustness using the actual retained subset and observed change above, without claiming that depth correction proves NP specificity. GEP5's contrasting AF context, FN1's expression in NP as well as AF, PRSS23's sparse core detection and CCN2's AF-weighted pattern remain in the unmodified regional tables. Low mouse expression does not establish absence in human disease.

NP and AF labels are provisional H&E cores corroborated by structural-marker direction, without expert pathology validation. No CEP was assigned. Genotype is fully confounded with capture, no cellular adjacency was tested, and neither a predicted receptor-complex path nor spatial expression establishes an FN1→ITGB1→PLOD2 causal chain.
"""
    (out / "SPATIAL_REVISED_INTERPRETATION.md").write_text(interpretation, encoding="utf8")
    print(json.dumps(status, indent=2), flush=True)


if __name__ == "__main__":
    main()
