#!/usr/bin/env python3
"""Project frozen programs onto sparse 48-um counts and fixed independent ROIs.

Scores include zeros and use the fixed available strict-one-to-one member set.
Only a requested gene column or score vector is materialized as a dense vector.
No ROI fitting, smoothing, genotype test, or bin-level significance test occurs.
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy import sparse

CANDIDATES = "FN1 ITGB1 PRSS23 CCN2 PLOD2 THBS1 LOX TGFBI GJA1 CD44 IGFBP3".split()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prepared48', required=True, type=Path)
    ap.add_argument('--roi-metadata', required=True, type=Path)
    ap.add_argument('--coverage', required=True, type=Path)
    ap.add_argument('--programs-file', required=True, type=Path)
    ap.add_argument('--programs', nargs='+', default=['GEP2', 'GEP5'])
    ap.add_argument('--candidates', nargs='+', default=CANDIDATES)
    ap.add_argument('--out', required=True, type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rois = pd.read_csv(a.roi_metadata, sep='\t', keep_default_na=False)
    coverage = pd.read_csv(a.coverage, sep='\t', keep_default_na=False)
    frozen = pd.read_csv(a.programs_file, sep='\t', keep_default_na=False)
    if rois.spot_id.duplicated().any() or frozen.duplicated(['program', 'gene']).any():
        raise ValueError('Ambiguous frozen spot/program member identity')
    program_rows, gene_rows, sources, mapping_rows = [], [], [], []
    folders = sorted(p.parent for p in a.prepared48.glob('*/counts48_symbols_qc.npz'))
    if not folders:
        raise ValueError('No actual prepared spatial sparse matrices')
    for folder in folders:
        x = sparse.load_npz(folder/'counts48_symbols_qc.npz').tocsr()
        obs = pd.read_csv(folder/'BIN48_METADATA_QC.tsv.gz', sep='\t')
        var = pd.read_csv(folder/'MOUSE_SYMBOL_FEATURES.tsv', sep='\t')
        if x.shape != (len(obs), len(var)) or not obs.spot_id.is_unique or not var.mouse_gene_symbol.is_unique:
            raise ValueError('Sparse count/metadata alignment or uniqueness failure')
        idx = dict(zip(var.mouse_gene_symbol, range(len(var))))
        anno = rois.set_index('spot_id').loc[obs.spot_id].reset_index()
        if anno.spot_id.tolist() != obs.spot_id.tolist():
            raise ValueError('ROI metadata order mismatch')
        for col in ['x_fullres', 'y_fullres']:
            if not np.allclose(anno[col].to_numpy(), obs[col].to_numpy(), atol=1e-9, rtol=0):
                raise ValueError('Original spatial coordinates changed')
        lib = np.asarray(x.sum(axis=1)).ravel().astype(float)
        if np.any(lib <= 0):
            raise ValueError('QC count matrix contains empty bins')
        z = x.astype(float).multiply((10000/lib)[:, None]).tocsr()
        z.data = np.log1p(z.data)
        local = coverage.loc[coverage.capture_id.eq(folder.name)].set_index('human_gene')
        if not local.index.is_unique:
            raise ValueError('Ambiguous human query in coverage')
        columns = ['spot_id', 'capture_id', 'genotype', 'disc_id', 'histology_core_region',
                   'anatomical_region', 'x_fullres', 'y_fullres', 'n_counts', 'n_genes']
        score_source = anno[columns].copy()
        programs = {}
        for program in a.programs:
            genes = frozen.loc[frozen.program.eq(program)].sort_values('rank').gene.tolist()
            if not genes:
                raise ValueError(f'Program missing from frozen table: {program}')
            matched = []
            for gene in genes:
                r = local.loc[gene]
                symbol = str(r.matched_native_mouse_symbol)
                exists = symbol in idx and r.mapping_status == 'MGI_ONE_TO_ONE'
                if exists:
                    matched.append(idx[symbol])
                mapping_rows.append(dict(capture_id=folder.name, program=program, human_gene=gene,
                                         mouse_symbol=symbol, strict_one_to_one_native_feature=exists,
                                         n_detected_qc_bins=int(r.n_detected_bins)))
            matched = list(dict.fromkeys(matched))
            scores = np.asarray(z[:, matched].mean(axis=1)).ravel() if matched else np.full(len(obs), np.nan)
            programs[program] = (scores, matched, len(genes))
            score_source[program+'_mean_log1p_CP10k'] = scores
        gene_values = {}
        for gene in a.candidates:
            r = local.loc[gene]
            symbol = str(r.matched_native_mouse_symbol)
            present = symbol in idx and r.mapping_status == 'MGI_ONE_TO_ONE'
            values = z[:, idx[symbol]].toarray().ravel() if present else np.full(len(obs), np.nan)
            raw = x[:, idx[symbol]].toarray().ravel() if present else np.full(len(obs), np.nan)
            gene_values[gene] = (values, raw, symbol, present)
            score_source[gene+'_log1p_CP10k'] = values
        for (disc, region), group in anno.groupby(['disc_id', 'anatomical_region'], sort=True):
            if region not in ['NP_CORE', 'AF_CORE']:
                continue
            take = group.index.to_numpy()
            base = dict(capture_id=folder.name, genotype=obs.genotype.iloc[0], disc_id=disc,
                        region=region, n_bins=len(take), median_UMI=float(np.median(lib[take])),
                        summary_unit='disc_nested_in_capture; animal_ID_unavailable',
                        between_genotype_inference='not_performed_capture_confounded')
            for program, (values, matched, frozen_n) in programs.items():
                detected = int(np.count_nonzero(np.asarray(x[take][:, matched].sum(axis=0)))) if matched else 0
                program_rows.append({**base, 'program': program,
                    'mean_program_score': float(np.mean(values[take])),
                    'median_program_score': float(np.median(values[take])),
                    'n_frozen_human_members': frozen_n, 'n_strict_1to1_native_members': len(matched),
                    'n_members_detected_in_region': detected,
                    'score_method': 'unweighted mean log1p CP10k over fixed available strict1to1 members; zeros retained'})
            for gene, (values, raw, symbol, present) in gene_values.items():
                detected = int(np.count_nonzero(raw[take])) if present else 0
                gene_rows.append({**base, 'human_gene': gene, 'mouse_symbol': symbol,
                    'native_feature_assessable': present, 'mean_log1p_CP10k': float(np.mean(values[take])),
                    'fraction_detected_bins': float(np.mean(raw[take]>0)) if present else np.nan,
                    'total_raw_UMI': float(np.sum(raw[take])) if present else np.nan,
                    'n_detected_bins': detected,
                    'interpretation': 'detectable' if detected else 'not_detectable_in_region' if present else 'not_assessable_mapping_or_feature'})
        sources.append(score_source)
    pd.DataFrame(program_rows).to_csv(a.out/'NP_PROGRAM_BY_REGION.tsv', sep='\t', index=False)
    pd.DataFrame(gene_rows).to_csv(a.out/'CANDIDATE_GENES_BY_REGION.tsv', sep='\t', index=False)
    pd.DataFrame(mapping_rows).to_csv(a.out/'PROGRAM_MEMBERS_USED_IN_SPATIAL_SCORE.tsv', sep='\t', index=False)
    pd.concat(sources, ignore_index=True).to_csv(a.out/'SPATIAL_PROGRAM_AND_GENE_SOURCE.tsv.gz', sep='\t', index=False)
    rois.to_csv(a.out/'DISC_AF_NP_ROI_METADATA.tsv', sep='\t', index=False)
    if a.prepared48.joinpath('SPATIAL_INPUT_QC.tsv').exists():
        target = a.out/'SPATIAL_INPUT_QC.tsv'
        if target.resolve() != a.prepared48.joinpath('SPATIAL_INPUT_QC.tsv').resolve():
            shutil.copy2(a.prepared48/'SPATIAL_INPUT_QC.tsv', target)
    status = dict(status='DONE_EXPLORATORY_REGION_QUANTIFICATION',
                  finished_utc=dt.datetime.now(dt.timezone.utc).isoformat(), captures=len(sources),
                  programs=a.programs, disc_program_rows=len(program_rows), candidate_region_rows=len(gene_rows),
                  histology_polygons_changed=False, genotype_p_values=False, imputation_or_smoothing=False,
                  frozen_programs_sha256=sha(a.programs_file), coverage_sha256=sha(a.coverage),
                  roi_metadata_sha256=sha(a.roi_metadata))
    (a.out/'SPATIAL_QUANTIFICATION_STATUS.json').write_text(json.dumps(status, indent=2)+'\n', encoding='utf8')
    print(json.dumps(status, indent=2))

if __name__ == '__main__':
    main()
