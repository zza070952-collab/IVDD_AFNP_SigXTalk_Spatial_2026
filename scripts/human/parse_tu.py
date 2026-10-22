#!/usr/bin/env python3
"""Stream the eight published Tu gene-by-cell TSVs into sparse matrices.

The published values are treated as processed/normalized expression, even when
integer-valued. Nothing is rounded or imputed. Original files remain read only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from array import array
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


MARKERS = {
    "NP_matrix": ["ACAN", "COL2A1", "SOX9", "KRT19", "CA12", "COMP", "FMOD", "CRTAC1"],
    "disc_ECM": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "DCN", "LUM", "COL12A1"],
    "immune": ["PTPRC", "LST1", "TYROBP", "FCER1G", "CD3D", "CD3E", "CD79A", "MS4A1"],
    "myeloid": ["LST1", "TYROBP", "FCER1G", "AIF1", "S100A8", "S100A9", "LYZ"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR", "CLDN5", "RAMP2"],
    "contractile": ["ACTA2", "TAGLN", "MYL9", "RGS5", "MCAM", "PDGFRB"],
}


def write_tsv(frame, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def load_sparse_tsv(path, gsm):
    values, indices = array("d"), array("i")
    pointers, genes = [0], []
    integer_values, max_value = True, 0.0
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        if header[0].strip().lower() != "gene":
            raise ValueError(f"{path}: expected first header field gene")
        cell_labels = header[1:]
        if not cell_labels or len(cell_labels) != len(set(cell_labels)):
            raise ValueError(f"{path}: missing/duplicate cell labels")
        for row_index, line in enumerate(handle):
            gene, separator, numeric = line.rstrip("\r\n").partition("\t")
            if not separator or not gene:
                raise ValueError(f"{path}:{row_index + 2}: malformed row")
            row = np.fromstring(numeric, sep="\t", dtype=np.float64)
            if row.size != len(cell_labels):
                raise ValueError(f"{path}:{row_index + 2}: wrong number of fields")
            if not np.isfinite(row).all() or (row < 0).any():
                raise ValueError(f"{path}:{row_index + 2}: nonfinite/negative expression")
            integer_values = integer_values and bool(np.equal(row, np.floor(row)).all())
            max_value = max(max_value, float(row.max(initial=0)))
            nonzero = np.flatnonzero(row)
            values.extend(row[nonzero])
            indices.extend(nonzero.astype(np.int32))
            pointers.append(len(values))
            genes.append(gene)
            if (row_index + 1) % 5000 == 0:
                print(json.dumps({"gsm": gsm, "parsed_genes": row_index + 1,
                                  "nnz": len(values)}), flush=True)
    matrix = sparse.csr_matrix((np.asarray(values), np.asarray(indices),
                                np.asarray(pointers, dtype=np.int64)),
                               shape=(len(genes), len(cell_labels)))
    # Duplicate submitted gene symbols are combined explicitly by summation.
    duplicates = len(genes) - len(set(genes))
    if duplicates:
        unique = list(dict.fromkeys(genes))
        row_lookup = {gene: i for i, gene in enumerate(unique)}
        transform = sparse.csr_matrix((np.ones(len(genes)),
                                       ([row_lookup[g] for g in genes], range(len(genes)))),
                                      shape=(len(unique), len(genes)))
        matrix = (transform @ matrix).tocsr().astype(np.float64)
        genes = unique
    return matrix, genes, cell_labels, {
        "all_submitted_values_integer": integer_values,
        "max_submitted_value": max_value,
        "duplicate_gene_rows_summed": duplicates,
        "published_matrix_nature": "normalized counts per official GSM metadata",
        "storage_dtype": "float64; no rounding/imputation",
    }


def match_genes(requested, genes, alias_rows):
    """Exact first, then unique case match; optional explicit aliases last."""
    lookup = {}
    for gene in genes:
        lookup.setdefault(gene.upper(), []).append(gene)
    aliases = {}
    for record in alias_rows:
        aliases.setdefault(record["original_gene"], []).append(record["mapped_gene"])
    output = []
    for original in requested:
        candidates = ([original] if original in genes else lookup.get(original.upper(), []))
        origin = "exact" if original in genes else "unique_case_match"
        if len(candidates) != 1:
            candidates = [x for x in aliases.get(original, []) if x in genes]
            origin = "explicit_alias_table"
        mapped = candidates[0] if len(candidates) == 1 else ""
        output.append((original, mapped, origin if mapped else "not_mapped"))
    return output


def annotate(matrix, genes):
    """Conservative marker-support labels, not a fitted cell atlas classifier."""
    gene_index = {g.upper(): i for i, g in enumerate(genes)}
    result = {}
    for panel, requested in MARKERS.items():
        rows = sorted({gene_index[g] for g in requested if g in gene_index})
        result[panel + "_detected_markers"] = np.asarray((matrix[rows] > 0).sum(axis=0)).ravel() if rows else np.zeros(matrix.shape[1], dtype=int)
        result[panel + "_panel_genes_available"] = len(rows)
        result[panel + "_mean_log1p_expression"] = np.asarray(matrix[rows].log1p().mean(axis=0)).ravel() if rows else np.zeros(matrix.shape[1])
    table = pd.DataFrame(result)
    immune = table.immune_detected_markers >= 2
    endothelial = table.endothelial_detected_markers >= 2
    np_support = table.NP_matrix_detected_markers >= 2
    stromal = (table.disc_ECM_detected_markers >= 2) | np_support
    labels = np.full(len(table), "unresolved_low_marker_support", dtype=object)
    labels[stromal] = "disc_stromal_marker_supported"
    labels[np_support & ~immune & ~endothelial] = "NP_matrix_marker_supported"
    labels[immune & ~stromal & ~endothelial] = "immune_marker_supported"
    labels[endothelial & ~immune & ~stromal] = "endothelial_marker_supported"
    labels[(immune & stromal) | (endothelial & stromal) | (immune & endothelial)] = "mixed_lineage_marker_support"
    table["identity"] = labels
    table["NP_matrix_support"] = np_support
    table["disc_stromal_support"] = stromal
    table["immune_support"] = immune
    table["endothelial_support"] = endothelial
    table["annotation_method"] = "transparent marker co-detection support; no forced cluster identity"
    return table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--frozen-signatures", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--alias-map", type=Path)
    parser.add_argument("--min-features", type=int, default=200)
    parser.add_argument("--max-mito-fraction", type=float, default=0.20)
    args = parser.parse_args()
    results = args.work / "01_human_programs"
    cache = args.work / "data_cache" / "Tu_sparse"
    results.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    signatures = pd.read_csv(args.frozen_signatures, sep="\t")
    if sorted(signatures.program.unique()) != ["GEP1", "GEP2", "GEP3", "GEP4", "GEP5"]:
        raise ValueError("Exactly the five frozen programs are required")
    if not signatures.groupby("program").size().eq(25).all():
        raise ValueError("Exactly 25 frozen genes/program are required")
    write_tsv(signatures, results / "FROZEN_K5_PROGRAMS.tsv")
    crosswalk = pd.read_csv(args.crosswalk, sep="\t").set_index("gsm")
    alias_rows = pd.read_csv(args.alias_map, sep="\t").to_dict("records") if args.alias_map else []
    all_metadata, sample_rows, mappings, coverage, marker_expression = [], [], [], [], []
    program_gene_expression = []
    for sample in range(1, 9):
        gsm = "GSM" + str(5048707 + sample)
        input_file = args.input_dir / f"{gsm}_Sample{sample}.counts.tsv.gz"
        if not input_file.exists():
            found = list(args.input_dir.rglob(input_file.name))
            if len(found) != 1:
                raise FileNotFoundError(input_file)
            input_file = found[0]
        matrix, genes, cells, details = load_sparse_tsv(input_file, gsm)
        sample_dir = cache / gsm
        sample_dir.mkdir(parents=True, exist_ok=True)
        cellname_file = input_file.with_name(f"{gsm}_Sample{sample}.cellname.txt.gz")
        if not cellname_file.exists():
            found = list(args.input_dir.rglob(cellname_file.name))
            if len(found) != 1:
                raise FileNotFoundError(cellname_file)
            cellname_file = found[0]
        barcode_table = pd.read_csv(cellname_file, sep="\t", dtype=str)
        if not {"CellName", "CellIndex"}.issubset(barcode_table.columns):
            raise ValueError(f"Unexpected companion barcode columns: {cellname_file}")
        if barcode_table.CellIndex.duplicated().any() or barcode_table.CellName.duplicated().any():
            raise ValueError(f"Duplicate companion barcode mapping: {cellname_file}")
        barcode_lookup = barcode_table.set_index("CellIndex").CellName
        if set(barcode_lookup.index) != set(cells):
            raise ValueError(f"Matrix/companion barcode identity mismatch: {gsm}")
        sparse.save_npz(sample_dir / "expression_gene_by_cell.npz", matrix)
        # MatrixMarket is an interoperable sparse input for the actual R UCell.
        mmwrite(str(sample_dir / "expression_gene_by_cell.mtx"), matrix)
        write_tsv(pd.DataFrame({"gene": genes}), sample_dir / "genes.tsv")
        cell_ids = [gsm + ":" + c for c in cells]
        submitted_barcodes = barcode_lookup.loc[cells].tolist()
        write_tsv(pd.DataFrame({"cell_id": cell_ids, "submitted_cell_index": cells,
                               "submitted_cell_barcode": submitted_barcodes}), sample_dir / "cells.tsv")
        observed = np.asarray((matrix > 0).sum(axis=0)).ravel()
        totals = np.asarray(matrix.sum(axis=0)).ravel()
        mt_rows = [i for i, g in enumerate(genes) if g.upper().startswith("MT-")]
        mito = np.asarray(matrix[mt_rows].sum(axis=0)).ravel() if mt_rows else np.zeros(len(cells))
        fraction = np.divide(mito, totals, out=np.zeros(len(cells), dtype=float), where=totals > 0)
        meta = annotate(matrix, genes)
        meta.insert(0, "cell_id", cell_ids)
        meta.insert(1, "gsm", gsm)
        meta.insert(2, "donor_id", "Tu_" + str(sample))
        meta["submitted_cell_index"] = cells
        meta["submitted_cell_barcode"] = submitted_barcodes
        meta["cohort"] = "Tu_GSE165722"
        meta["tissue"] = "NP surgical tissue"
        meta["expression_sum_published"] = totals
        meta["n_features_detected"] = observed
        meta["mitochondrial_expression_fraction"] = fraction
        meta["qc_pass"] = (observed >= args.min_features) & (fraction <= args.max_mito_fraction) & (totals > 0)
        for field, value in crosswalk.loc[gsm].items():
            if field not in meta:
                meta[field] = value
        write_tsv(meta, sample_dir / "metadata.tsv.gz")
        all_metadata.append(meta)
        details.update({"gsm": gsm, "input_path": str(input_file), "input_size_bytes": input_file.stat().st_size,
                        "n_genes": len(genes), "n_cells": len(cells), "nnz": matrix.nnz,
                        "n_qc_pass": int(meta.qc_pass.sum()), "n_mito_genes_available": len(mt_rows),
                        "qc_min_features": args.min_features, "qc_max_mito_fraction": args.max_mito_fraction})
        sample_rows.append(details)
        gene_index = {g: i for i, g in enumerate(genes)}
        for program, group in signatures.groupby("program", sort=True):
            mapped = match_genes(group.sort_values("rank").gene.tolist(), genes, alias_rows)
            for original, mapped_gene, basis in mapped:
                row = {"gsm": gsm, "program": program, "original_gene": original, "mapped_gene": mapped_gene, "mapping_basis": basis}
                if mapped_gene:
                    gene_row = matrix.getrow(gene_index[mapped_gene])
                    row["n_detected_all_cells"] = gene_row.nnz
                    row["fraction_detected_all_cells"] = gene_row.nnz / len(cells)
                else:
                    row["n_detected_all_cells"] = 0
                    row["fraction_detected_all_cells"] = 0.0
                mappings.append(row)
                for identity in sorted(meta.identity.unique()):
                    mask = (meta.identity == identity).to_numpy() & meta.qc_pass.to_numpy()
                    if not mask.any():
                        continue
                    vec = matrix.getrow(gene_index[mapped_gene]).toarray().ravel() if mapped_gene else None
                    program_gene_expression.append({"gsm": gsm, "donor_id": "Tu_" + str(sample),
                        "program": program, "original_gene": original, "mapped_gene": mapped_gene,
                        "identity": identity, "n_cells": int(mask.sum()), "gene_available": bool(mapped_gene),
                        "mean_log1p_published_expression": float(np.log1p(vec[mask]).mean()) if mapped_gene else np.nan,
                        "fraction_detected": float((vec[mask] > 0).mean()) if mapped_gene else np.nan})
            present = [g for _, g, _ in mapped if g]
            detectable = sum(matrix.getrow(gene_index[g]).nnz > 0 for g in present)
            coverage.append({"gsm": gsm, "donor_id": "Tu_" + str(sample), "program": program,
                             "n_frozen_genes": 25, "n_mapped_genes": len(present),
                             "fraction_mapped": len(present) / 25,
                             "n_detected_genes_any_cell": detectable,
                             "technical_coverage_status": "assessable" if len(present) >= 10 else "low_coverage"})
        for panel, markers in MARKERS.items():
            for gene in markers:
                rows = [i for i, g in enumerate(genes) if g.upper() == gene]
                if not rows:
                    continue
                vec = matrix.getrow(rows[0]).toarray().ravel()
                for identity in sorted(meta.identity.unique()):
                    mask = (meta.identity == identity).to_numpy() & meta.qc_pass.to_numpy()
                    if mask.any():
                        marker_expression.append({"gsm": gsm, "donor_id": "Tu_" + str(sample), "identity": identity,
                                                  "marker_panel": panel, "gene": gene, "n_cells": int(mask.sum()),
                                                  "mean_log1p_published_expression": float(np.log1p(vec[mask]).mean()),
                                                  "fraction_detected": float((vec[mask] > 0).mean())})
        print(json.dumps({"sample_complete": gsm, "cells": len(cells), "qc_pass": int(meta.qc_pass.sum()),
                          "matrix_nature": details["published_matrix_nature"]}), flush=True)
    metadata = pd.concat(all_metadata, ignore_index=True)
    write_tsv(metadata, results / "TU_CELL_METADATA.tsv.gz")
    write_tsv(pd.DataFrame(sample_rows), results / "TU_MATRIX_AND_QC_MANIFEST.tsv")
    write_tsv(pd.DataFrame(mappings), results / "PROGRAM_GENE_MAPPING.tsv")
    write_tsv(pd.DataFrame(coverage), results / "PROGRAM_GENE_COVERAGE.tsv")
    write_tsv(pd.DataFrame(marker_expression), results / "TU_MARKER_DOTPLOT_SOURCE.tsv")
    write_tsv(pd.DataFrame(program_gene_expression), results / "TU_PROGRAM_GENE_DOTPLOT_SOURCE.tsv")
    donor_table = crosswalk.reset_index().copy()
    donor_table["donor_id"] = ["Tu_" + str(i) for i in range(1, 9)]
    donor_table = donor_table.merge(pd.DataFrame(sample_rows)[["gsm", "n_cells", "n_qc_pass", "n_genes"]], on="gsm", validate="one_to_one")
    write_tsv(donor_table, results / "COHORT_AND_DONOR_TABLE.tsv")
    write_tsv(pd.DataFrame([{"panel": panel, "gene": gene} for panel, genes in MARKERS.items() for gene in genes]), results / "ANNOTATION_MARKER_PANEL.tsv")
    (results / "TU_PARSE_STATUS.json").write_text(json.dumps({"status": "DONE", "samples": 8,
        "cells_total": len(metadata), "cells_qc_pass": int(metadata.qc_pass.sum()),
        "clinical_grade_crosswalk": "unresolved; no grade-group inferential test",
        "expression_values": "published normalized counts, no raw-count assumption",
        "annotation_boundary": "marker-support labels; mixed and unresolved cells retained"}, indent=2))


if __name__ == "__main__":
    main()
