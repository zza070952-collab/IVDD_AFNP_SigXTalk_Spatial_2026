from pathlib import Path
import csv, gzip, hashlib, json, re
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('--project-root',type=Path,required=True,help='Read-only V4 source cache')
ap.add_argument('--output-root',type=Path,required=True,help='Independent private artifact output')
ap.add_argument('--v3-review-root',type=Path,required=True)
ap.add_argument('--work-dir',type=Path)
a=ap.parse_args()
ROOT=a.project_root.resolve()
V3=a.v3_review_root.resolve()
OUT=a.output_root.resolve()
if OUT==ROOT or OUT in ROOT.parents or ROOT in OUT.parents:raise ValueError('Output must be independent of the read-only V4 input cache')
(OUT/'06_supplementary').mkdir(parents=True,exist_ok=True)
WORK=(a.work_dir or OUT/'02_analysis'/'workbook_build').resolve()
WORK.mkdir(parents=True,exist_ok=True)
HUM=ROOT/'03_tables'/'human'; CS=ROOT/'03_tables'/'communication_spatial'; EXP=ROOT/'03_tables'/'experimental'/'internal'/'source_unresolved'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        rd=csv.reader(f,delimiter='\t');rows=list(rd)
    assert rows and rows[0],p
    assert all(len(r)==len(rows[0]) for r in rows),p
    return rows[0],rows[1:]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def public_text(v):
    # Export a portable logical reference rather than a user or host path.
    if re.search(r"(?:[A-Za-z]:[/\\]|/(?:home|Users)/)",str(v)):
        return 'source_file/'+str(v).replace('\\','/').rsplit('/',1)[-1]
    return v

missing={'','NA','NaN','nan','NAN','None','null'}
numeric=re.compile(r'^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$')
ident=re.compile(r'(^id$|_id$|_key$|^barcode$|^gsm$|^gene$|^symbol$|^donor$|^patient$|^batch$|^group$|^capture$|^disc$|^sample$|^source$|^target$|^receptor$|^ssc$|^ligand$|^program$|^sequence$)',re.I)
def kind(col,vals):
    keep=[v for v in vals if v not in missing]
    if not keep or ident.search(col):return 'text'
    if all(re.match(r'^\d{4}-\d{2}-\d{2}(?:T[0-9:.+-]+Z?)?$',v) for v in keep):return 'datetime'
    if all(v.lower() in ['true','false'] for v in keep):return 'boolean'
    if all(numeric.match(v) and not re.match(r'^0\d+',v) for v in keep):
        return 'integer' if all(float(v).is_integer() and abs(float(v))<2**53 for v in keep) else 'number'
    return 'text'
def unit(col,k):
    c=col.lower()
    if k=='datetime':return 'UTC timestamp or source calendar date'
    if 'sha' in c:return 'SHA256 identifier'
    if 'ct' in c and any(x in c for x in ['deltact','mean_ct','ct_','_ct','cq']):return 'PCR cycles'
    if 'od450' in c or 'od_' in c or 'blank_corrected' in c:return 'OD450 absorbance'
    if re.search(r'(^p$|^p_|_p$|p_adjusted|p_two|q_value|fdr|pvalue|p_value)',c):return 'probability'
    if 'fold' in c or 'ratio' in c:return 'dimensionless ratio'
    if 'fraction' in c or 'proportion' in c or 'retention' in c:return 'dimensionless fraction unless source header specifies percent'
    if 'hours' in c or c.endswith('_h'):return 'hours'
    if c.endswith('_nm'):return 'nM'
    if 'score' in c or 'ucell' in c or 'aucell' in c:return 'method-specific program score; UCell/AUCell not comparable across methods'
    if 'prs' in c:return 'official pathway regulation strength'
    if 'fidelity' in c or 'specificity' in c:return 'official SigXTalk metric'
    if any(x in c for x in ['pixel','x_full','y_full','image_x','image_y']):return 'image pixels; see coordinate source table'
    if any(x in c for x in ['count','n_cells','umi','n_bins','n_genes','n_paths','n_batches']):return 'count; nested records are not independent donors/animals'
    if k=='text':return 'identifier, category, status or source-defined text'
    return 'source-defined quantity; see column name and sheet source'
def metadata(path):
    if path.is_relative_to(ROOT):return path.relative_to(ROOT).as_posix()
    if path.is_relative_to(OUT):return 'export/'+path.relative_to(OUT).as_posix()
    return 'V3/'+path.relative_to(V3).as_posix()
humanmap={r['panel_id']:r for r in csv.DictReader((HUM/'HUMAN_DATA_TO_FIGURE_MAP.tsv').open(encoding='utf-8-sig'),delimiter='\t')}
commmap={r['figure']+r['panel_id']:r for r in csv.DictReader((ROOT/'08_review'/'communication_spatial'/'PANEL_SOURCE_DATA_MAP.tsv').open(encoding='utf-8-sig'),delimiter='\t')}
def add_sheet(name,path,internal=False,panel=None):
    heads,rows=read(path)
    types=[kind(h,[r[i] for r in rows]) for i,h in enumerate(heads)]
    out=[]
    for raw in rows:
        row=[]
        for h,v,t in zip(heads,raw,types):
            if v=='':row.append(None)
            elif t in ['number','integer'] and v not in missing:row.append(int(float(v)) if t=='integer' else float(v))
            elif t=='boolean' and v not in missing:row.append(v.lower()=='true')
            else:row.append(public_text(v) if not internal else v)
        out.append(row)
    if internal and 'SOURCE_STATUS' not in heads:
        heads.append('SOURCE_STATUS');types.append('text')
        for row in out:row.append('INTERNAL_REVIEW_SOURCE_UNRESOLVED')
    m=humanmap.get(panel,{}) if panel else {}
    if not m:m=commmap.get(panel,{}) if panel else {}
    return {'name':name,'source':metadata(path),'source_sha256':sha(path),'headers':heads,'rows':out,'types':types,'units':[unit(h,t) for h,t in zip(heads,types)],'panel':panel or 'supplementary table','analysis_unit':m.get('unit',m.get('statistical_unit','See donor/batch/disc/capture fields; rows are not automatically independent')),'independent_n':m.get('independent_n',m.get('independent_n','See source metadata and manuscript')),'transformation':m.get('transform',m.get('transformation','Static typed export of script-generated source table')),'test':m.get('statistical_test','See source statistical-result columns; no additional Excel test'),'status':'INTERNAL_REVIEW_SOURCE_UNRESOLVED' if internal else 'COMPUTATIONAL_SOURCE_VERIFIED','missing':'Blank retained as unavailable source entry. NA/NaN and explicit state fields retained. Zero is not a substitute for missing data.','typed_matrix_sha256':hashlib.sha256(json.dumps([heads]+out,ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()).hexdigest()}

supp=[]
for file,name in [
('S1_PUBLIC_COHORTS.tsv','S1_Cohorts'),('S2_DONOR_IDENTITY_COUNTS.tsv','S2_DonorCounts'),('S2_TU_QC.tsv','S2_TuQC'),('S2_COMPLETE_MARKER_EXPRESSION.tsv','S2_Markers'),('S2_COMMUNICATION_DONOR_ELIGIBILITY.tsv','S2_CommunicationDonors'),('S3_ALL125_FROZEN_PROGRAM_MEMBERS.tsv','S3_All125Programs'),('S3_CNMf_FROZEN_K_SELECTION_PARAMETERS.tsv','S3_FrozenKParameters'),('S3_PROGRAM_COVERAGE.tsv','S3_ProgramCoverage'),('S3_GENE_MAPPING.tsv','S3_GeneMapping'),('S3_AUCELL_EFFECTIVE_DENOMINATORS.tsv','S3_AUCellCoverage'),('S4_DONOR_PROGRAM_SCORES.tsv','S4_DonorScores'),('S4_GEP2_DONOR_SUPPORT.tsv','S4_GEP2Support'),('S4_PROGRAM_INTERPRETATION.tsv','S4_ProgramInterpretation'),('S4_ALL_PROGRAM_GENE_EXPRESSION.tsv','S4_AllProgramGeneExpression'),('S4_CANDIDATE_EXPRESSION.tsv','S4_CandidateExpression')]:supp.append(add_sheet(name,HUM/file))
for file,name in [
('S5_S6_COMMUNICATION_COUNTS.tsv','S5S6_RunCounts'),('S5_ALL_AF_TO_NP_LR.tsv','S5_All103AF_NP_LR'),('S5_ALL_SOURCE_TARGET_LR.tsv','S5_AllSourceTargetLR'),('S5_FULL_TARGET_SELECTION_UNIVERSE.tsv','S5_FullTargetUniverse'),('S5_DONOR_GROUP_COUNTS.tsv','S5_DonorCellInput'),('S5_LR_COMPLEX_COMPONENTS.tsv','S5_ComplexComponents'),('S6_ALL_25_GEP2_TESTABILITY_AND_PRS_STATES.tsv','S6_All25TargetStates'),('S6_CANDIDATE_TESTABILITY.tsv','S6_CandidateTestability'),('S6_REPRESENTATIVE_OFFICIAL_METRICS.tsv','S6_OfficialMetrics'),('S6_RECEPTOR_NODE_AXES_AND_LR_COMPONENTS.tsv','S6_ReceptorNodeAudit'),('S7_ALL_PROGRAM_ORTHOLOG_COVERAGE.tsv','S7_OrthologCoverage'),('S7_ALL_PROGRAM_COVERAGE_SUMMARY.tsv','S7_CoverageSummary'),('S7_ALL_PROGRAM_MEMBERS_USED.tsv','S7_AllProgramMembers'),('S7_CAPTURE_INPUT_QC.tsv','S7_CaptureQC'),('S7_FROZEN_ROI_PROVENANCE.tsv','S7_FrozenROI'),('S7_HUMAN_MOUSE_STRICT_MAPPING.tsv','S7_HumanMouseMap'),('S7_ALL_INDEPENDENT_MARKERS.tsv','S7_IndependentMarkers'),('S7_IMAGE_INTEGRITY_PROVENANCE.tsv','S7_ImageProvenance'),('S7_MARKER_COMPATIBILITY.tsv','S7_MarkerCompatibility'),('S8_ALL_FIVE_PROGRAMS_DISC_REGION.tsv','S8_AllProgramsDiscRegion'),('S8_CANDIDATE_DISC_REGION_VALUES.tsv','S8_CandidatesDiscRegion'),('S8_DEPTH_DISC_PAIRS.tsv','S8_DepthDiscPairs'),('S8_DEPTH_RETENTION_EACH_DISC_REGION.tsv','S8_DepthRetention'),('S8_SINGLE_FIXED_DRAW_BIN_VALUES.tsv','S8_SingleDrawBins')]:supp.append(add_sheet(name,CS/file))
s9source=V3/'04_integrated'/'MAIN_CANDIDATE_EVIDENCE_MATRIX.tsv'
s9full=OUT/'06_supplementary'/'Supplementary_Data_3_Complete_Candidate_Evidence.tsv.gz'
with s9full.open('wb') as out:
    with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as gz:gz.write(s9source.read_bytes())
assert gzip.decompress(s9full.read_bytes())==s9source.read_bytes()
s9=add_sheet('S9_CandidateEvidence',s9source)
excluded=[h for h in s9['headers'] if h.endswith('_json') or h in ['role_specific_database_SSC_symbols','role_specific_receiver_targets','role_specific_path_keys']]
take=[i for i,h in enumerate(s9['headers']) if h not in excluded]
for field in ['headers','types','units']:s9[field]=[s9[field][i] for i in take]
s9['rows']=[[r[i] for i in take] for r in s9['rows']]
s9['source']=s9full.relative_to(OUT).as_posix();s9['source_sha256']=sha(s9full)
s9['transformation']='All 10 candidates and scalar evidence retained. Nested JSON and long member/path lists indexed in S9_ExternalFields and supplied without changes in Supplementary Data 3.'
s9['typed_matrix_sha256']=hashlib.sha256(json.dumps([s9['headers']]+s9['rows'],ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()).hexdigest()
supp.append(s9)
external_head=['excluded_field','reason','complete_source_file','compressed_SHA256','uncompressed_TSV_SHA256']
external_rows=[[h,'Nested/long text supplied externally to preserve full content and Excel cell limits',s9['source'],sha(s9full),sha(s9source)] for h in excluded]
supp.append({'name':'S9_ExternalFields','source':s9['source'],'source_sha256':sha(s9full),'headers':external_head,'rows':external_rows,'types':['text']*5,'units':['source field name','exclusion rationale','relative package reference','SHA256 identifier','SHA256 identifier'],'panel':'Supplementary Table S9','analysis_unit':'10 candidate genes','independent_n':'Role-specific evidence; no candidate score','transformation':'Index of externally supplied fields; full gzip decompresses exactly to authoritative TSV','test':'none','status':'COMPUTATIONAL_SOURCE_VERIFIED','missing':'No fields deleted from complete gzip source','typed_matrix_sha256':hashlib.sha256(json.dumps([external_head]+external_rows,ensure_ascii=False,separators=(',',':')).encode()).hexdigest()})

main=[]
for fig,panels in [(1,'ABCDE'),(2,'ABCDEF')]:
    for panel in panels:main.append(add_sheet(f'Figure{fig}{panel}',HUM/f'Figure{fig}{panel}_source.tsv',panel=f'Figure{fig}{panel}'))
for key,r in commmap.items():
    if r['figure'] not in ['Figure3','Figure4']:continue
    file=ROOT/Path(r['source_data_file'].replace('\\','/'))
    main.append(add_sheet(key,file,panel=key))

internal_files=[
('QPCR_CONTRASTS_AND_CI.tsv','qPCR_Contrasts_CI'),('QPCR_BATCH_CONTRASTS.tsv','qPCR_BatchContrasts'),('QPCR_BATCH_DELTACt.tsv','qPCR_BatchDeltaCt'),('QPCR_TECHNICAL_MEANS.tsv','qPCR_TechnicalMeans'),('QPCR_NTC_NO_RT_RECORDED.tsv','qPCR_NTC_NoRT'),('CCK8_BATCH_SUMMARY.tsv','CCK8_BatchSummary'),('CCK8_BLANK_CORRECTED.tsv','CCK8_BlankCorrected'),('CCK8_MATCHED_BLANKS.tsv','CCK8_MatchedBlanks'),('EXPERIMENT_DESIGN.tsv','ExperimentalDesign'),('EXPERIMENT_QC_AVAILABILITY.tsv','QC_Availability'),('SUPPLIED_SEQUENCES_EXACT.tsv','ReagentSequences'),('siRNA_qPCR_Raw.tsv','Raw_siRNA_qPCR'),('Program_qPCR_Raw.tsv','Raw_Program_qPCR'),('CCK8_Raw.tsv','Raw_CCK8')]
internal_files += [('SUPPLIED_REAGENTS_AND_CONDITIONS.tsv','SuppliedReagents'),('SUPPLIED_PROTOCOL_CONDITIONS.tsv','SuppliedProtocol')]
cells=[add_sheet(name,EXP/file,True) for file,name in internal_files]
s10=[add_sheet('S10_Design',EXP/'EXPERIMENT_DESIGN.tsv',True),add_sheet('S10_Sequences',EXP/'SUPPLIED_SEQUENCES_EXACT.tsv',True),add_sheet('S11_qPCR_DeltaCt',EXP/'QPCR_BATCH_DELTACt.tsv',True),add_sheet('S11_qPCR_Contrasts',EXP/'QPCR_CONTRASTS_AND_CI.tsv',True),add_sheet('S11_qPCR_BatchContrasts',EXP/'QPCR_BATCH_CONTRASTS.tsv',True),add_sheet('S11_AssayQC',EXP/'EXPERIMENT_QC_AVAILABILITY.tsv',True),add_sheet('S11_NTC_NoRT',EXP/'QPCR_NTC_NO_RT_RECORDED.tsv',True),add_sheet('S12_CCK8_MatchedBlanks',EXP/'CCK8_MATCHED_BLANKS.tsv',True),add_sheet('S12_CCK8_BlankCorrected',EXP/'CCK8_BLANK_CORRECTED.tsv',True),add_sheet('S12_CCK8_BatchSummary',EXP/'CCK8_BATCH_SUMMARY.tsv',True)]
s10.insert(2,add_sheet('S10_Reagents',EXP/'SUPPLIED_REAGENTS_AND_CONDITIONS.tsv',True))
s10.insert(3,add_sheet('S10_Protocol',EXP/'SUPPLIED_PROTOCOL_CONDITIONS.tsv',True))
specs=[{'file':'06_supplementary/Supplementary_Tables.xlsx','title':'Supplementary Tables S1-S9','internal':False,'sheets':supp},{'file':'06_supplementary/Source_Data_Main_Figures.xlsx','title':'Source data for Figures 1-4','internal':False,'sheets':main},{'file':'06_supplementary/internal/source_unresolved/Source_Data_Cell_Experiments.xlsx','title':'Cell-experiment source values for internal review','internal':True,'sheets':cells},{'file':'06_supplementary/internal/source_unresolved/Supplementary_Tables_S10-S12.xlsx','title':'Supplementary Tables S10-S12 for internal review','internal':True,'sheets':s10}]
(WORK/'workbook_payload.json').write_text(json.dumps({'root':str(OUT),'books':specs},ensure_ascii=False,allow_nan=False,separators=(',',':')),encoding='utf-8')
summary=[{'file':b['file'],'source_sheets':len(b['sheets']),'rows':sum(len(s['rows']) for s in b['sheets']),'cells':sum((len(s['rows'])+1)*len(s['headers']) for s in b['sheets'])} for b in specs]
(WORK/'INPUT_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
