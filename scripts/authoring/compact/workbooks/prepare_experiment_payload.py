"""Lossless raw observations and frozen summaries; no statistical recomputation."""
import argparse,csv,hashlib,json,re
from pathlib import Path
from docx import Document
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--project',type=Path,required=True)
ap.add_argument('--revision',type=Path,required=True)
ap.add_argument('--protocol',type=Path,required=True)
a=ap.parse_args();R=a.project.resolve();C=a.revision.resolve();E=R/'03_tables/experimental/internal/source_unresolved'
(C/'source_tables').mkdir(parents=True,exist_ok=True);(C/'audit').mkdir(exist_ok=True)
STATUS='INTERNAL_REVIEW_SOURCE_UNRESOLVED'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(n):
 p=E/n
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def number(x):return None if x in ['',None,'NA','nan'] else float(x)
def sheet(name,headers,rows,types,description,source):
 return dict(name=name,headers=headers,rows=rows,types=types,description=description,source=source,source_sha256=sha(E/source) if (E/source).exists() else '',status=STATUS)
raw=json.loads((E/'RAW_EXTRACT.json').read_text(encoding='utf-8'))
qp=raw['siRNA_qPCR_Raw']+raw['Program_qPCR_Raw'];assert len(qp)==632
h=list(qp[0]);qr=[[r[k] for k in h]+[STATUS] for r in qp]
types=['integer' if k in ['FN_hours','Tech_rep','Source_excel_row'] else 'number' if k=='Raw_Ct' else 'text' for k in h]+['text']
q=sheet('Raw_qPCR_Ct',h+['Source_status'],qr,types,'All 632 supplied wells: 576 sample, 24 NTC and 32 No-RT records. Blank Ct plus Cq_status=Undetermined is retained; never Ct=40. B01–B04 are batches of one commercial immortalized cell source, not donors; technical wells are not independent replicates. Original input checksums belong to the private extraction receipt.','RAW_EXTRACT.json')
blank={r['Blank_key']:r for r in read('CCK8_MATCHED_BLANKS.tsv')}
cc=raw['CCK8_Raw'];h2=list(cc[0]);ccr=[]
for r in cc:
 b=blank[r['Blank_key']]
 ccr.append([r[k] for k in h2]+[number(b['matched_blank_mean_OD450']),int(r['FN_hours'])+2,STATUS])
ct=['integer' if k in ['FN_hours','Tech_rep','Source_excel_row'] else 'number' if k=='Raw_OD450' else 'text' for k in h2]+['number','integer','text']
c=sheet('Raw_CCK8_OD450',h2+['Matched_blank_mean_OD450','Actual_read_hours','Source_status'],ccr,ct,'All 288 supplied cell and matched cell-free wells. FN_hours=CCK8 addition at 4/24/48 h after reseeding; actual reads at 6/26/50 h. Matched blanks use the same batch, time and coating. No deletion or clipping of observations.','RAW_EXTRACT.json')
headers=['Record_type','Batch_ID','Group_ID','siRNA','FN_condition','Gene','Contrast','Measure','Value','CI95_low','CI95_high','P_two_sided','P_adjusted','Adjustment','Independent_n','Technical_n','Units','Definition','Source_file','Source_status']
rows=[]
def add(kind,r,measure,value,units,definition,src,lo=None,hi=None,p=None,padj=None,adj='NA',n=4,tn=None):
 rows.append([kind,r.get('Batch_ID','NA'),r.get('Group_ID','NA'),r.get('siRNA','NA'),r.get('FN_condition','NA'),r.get('Gene',r.get('gene','NA')),r.get('contrast','NA'),measure,number(value),number(lo),number(hi),number(p),number(padj),adj,n,tn,units,definition,src,STATUS])
for r in read('QPCR_TECHNICAL_MEANS.tsv'):
 for measure,unit,definition in [('mean_Ct','cycles','Arithmetic mean of the recorded technical Ct wells'),('technical_sd','cycles','Sample SD across the recorded technical wells; not instrument assay validation')]:
  add('qPCR_technical_summary',r,measure,r[measure],unit,definition,'QPCR_TECHNICAL_MEANS.tsv',tn=int(r['technical_n']))
for r in read('QPCR_BATCH_DELTACt.tsv'):
 for measure,unit,definition in [('mean_reference_Ct','cycles','Mean of the same-batch/group RPLP0 and TBP technical means'),('deltaCt','cycles','Target mean Ct minus mean reference Ct'),('deltaDeltaCt','cycles','DeltaCt minus same-batch G1 deltaCt'),('fold_vs_same_batch_G1','fold','2^(-deltaDeltaCt), display only')]:
  add('qPCR_batch_result',r,measure,r[measure],unit,definition,'QPCR_BATCH_DELTACt.tsv',tn=int(r['technical_n']))
for r in read('QPCR_BATCH_CONTRASTS.tsv'):
 # Preserve every cached batch-level contrast using its original field label.
 for k in ['deltaCt_contrast','fold_ratio']:
  if k in r:add('qPCR_batch_contrast',r,k,r[k],'cycles' if k=='deltaCt_contrast' else 'fold','Cached prespecified batch contrast; D1=(G4-G3)-(G2-G1), D2=(G6-G5)-(G2-G1)','QPCR_BATCH_CONTRASTS.tsv')
for r in read('QPCR_CONTRASTS_AND_CI.tsv'):
 add('qPCR_contrast_statistics',r,'fold_ratio',r['fold_ratio'],'fold','Exponentiated mean DeltaCt contrast; CI bounds reversed on exponentiation. Primary PLOD2 D1/D2 Holm family=2; secondary BH family=20; FN_siNC supportive PLOD2 contrast unadjusted.','QPCR_CONTRASTS_AND_CI.tsv',r['fold_CI95_low'],r['fold_CI95_high'],r['p_two_sided'],r['p_adjusted'],r['adjustment'],int(r['n_batches']))
 add('qPCR_contrast_statistics',r,'mean_deltaCt_contrast',r['mean_deltaCt_contrast'],'cycles','Two-sided one-sample t test of four batch contrasts; df=3. This numeric analysis does not authenticate source.','QPCR_CONTRASTS_AND_CI.tsv',r['deltaCt_CI95_low'],r['deltaCt_CI95_high'],r['p_two_sided'],r['p_adjusted'],r['adjustment'],int(r['n_batches']))
for r in read('CCK8_BATCH_SUMMARY.tsv'):
 r['contrast']='addition_'+r['FN_hours']+'h_read_'+r['actual_read_hours_after_reseeding']+'h'
 for measure,unit,definition in [('mean_corrected_OD450','OD450','Mean cell OD450 minus matched blank mean'),('mean_raw_OD450','OD450','Arithmetic mean of the three cell wells'),('matched_blank_mean_OD450','OD450','Mean of matched batch/time/coating blank wells'),('ratio_vs_same_batch_time_G1','ratio','Corrected group mean divided by corrected same-batch/time G1'),('ratio_vs_own_4h','ratio','Corrected group mean divided by its own same-batch/group 4h-addition value')]:
  add('CCK8_batch_result',r,measure,r[measure],unit,definition,'CCK8_BATCH_SUMMARY.tsv',tn=int(r['tech_n']))
stats=sheet('Batch_Results_Statistics',headers,rows,['number' if k in ['Value','CI95_low','CI95_high','P_two_sided','P_adjusted'] else 'integer' if k in ['Independent_n','Technical_n'] else 'text' for k in headers],'All cached qPCR batch results and prespecified contrasts plus all CCK8 batch/time/group summaries; static verified values. NA indicates inapplicable fields, not zero. n=4 batch identifiers; 3 technical wells are averaged. Missing instrument QC and source linkage remain unresolved.','QPCR_CONTRASTS_AND_CI.tsv')

# Extract sequence pairs within the SAME original DOCX paragraph, not adjacent paragraphs.
doc=a.protocol.resolve()
paras=[p.text for p in Document(doc).paragraphs]
s6head=['Record_type','Record_ID','Gene_or_group','siRNA','FN_condition','Forward_or_sense_5to3','Reverse_or_antisense_5to3','PrimerBank_ID','Amplicon_bp','Parameter','Value','Units','Source_paragraph_1based','Source_SHA256','Source_status','Notes']
s6=[]
for i,p in enumerate(paras):
 m=re.match(r'^(siITGB1_[12]|siNC)\s+S\s+([ACGU]+dTdT)\s+AS\s+([ACGU]+dTdT)\s*$',p)
 if m:s6.append(['siRNA',m[1],'ITGB1' if m[1]!='siNC' else 'Negative_control',m[1],'NA',m[2],m[3],'NA',None,'Stock',20,'uM',i+1,sha(doc),STATUS,'Both strands 5 to 3; exactly copied from the same supplied paragraph; assay authenticity pending'])
 m=re.match(r'^([^（]+)（([0-9]+c[0-9]+)，([0-9]+) bp）\s+F\s+([ACGT]+)\s+R\s+([ACGT]+)\s*$',p)
 if m:
  gene='ITGB1_c3_region_check' if '区域' in m[1] else m[1]
  s6.append(['qPCR_primer',gene,gene,'NA','NA',m[4],m[5],m[2],int(m[3]),'Working',10,'uM',i+1,sha(doc),STATUS,'ITGB1 c2 overlaps si2 target; independent c3 validation records unavailable' if gene.startswith('ITGB1') else 'Exact supplied sequence; no independent assay efficiency or specificity record supplied'])
assert len([r for r in s6 if r[0]=='siRNA'])==3
assert len([r for r in s6 if r[0]=='qPCR_primer'])==9
seen=set()
for r in read('EXPERIMENT_DESIGN.tsv'):
 if r['Group_ID'] in seen:continue
 seen.add(r['Group_ID'])
 s6.append(['Experimental_group',r['Group_ID'],r['Group_ID'],r['siRNA'],r['FN_condition'],'NA','NA','NA',None,'siRNA',10,'nM',None,sha(doc),STATUS,'Four batch identifiers B01-B04 from one SV40T immortalized commercial source; RNA at 72h post-transfection / 24h post-reseeding'])
# Optional protocol-condition fields are provided as a private external table;
# they are not filled from unverified hard-coded study prose.
conditions=E/'SUPPLIED_PROTOCOL_CONDITIONS.tsv'
if conditions.exists():
 for row in read(conditions.name):
  key=row.get('step',row.get('condition',row.get('item',row.get('parameter','recorded_condition'))))
  value=row.get('documented_condition',row.get('value',row.get('recorded_value',row.get('supplied_text',''))))
  if value:s6.append(['Core_condition',key,'NA','NA','NA','NA','NA','NA',None,key,value,row.get('units','source-defined'),None,sha(doc),STATUS,'Copied from private protocol-condition input; execution is not authenticated'])
six=dict(name='S6_Oligos_and_Design',headers=s6head,rows=s6,types=['integer' if x in ['Amplicon_bp','Source_paragraph_1based'] else 'text' for x in s6head],description='Exact within-paragraph sequence pairs from the supplied protocol; six-group design and essential conditions. Source unresolved. No purchasing list. Primer sequences appear only here.',source='Supplied experimental protocol',source_sha256=sha(doc),status=STATUS)

out={'experiment_sheets':[q,c,stats],'S6_sheet':six,'raw_counts':{k:len(v) for k,v in raw.items()},'sequence_pairs':12}
(C/'source_tables/EXPERIMENT_PAYLOAD.json').write_text(json.dumps(out,ensure_ascii=False,allow_nan=False,indent=2),encoding='utf-8')

for s in out['experiment_sheets']+[six]:
 with (C/'source_tables'/(s['name']+'.tsv')).open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(s['headers']);w.writerows(s['rows'])
print(json.dumps({'sheets':[{k:s[k] for k in ['name']}|{'rows':len(s['rows']),'columns':len(s['headers'])} for s in out['experiment_sheets']+[six]],'sequence_pairs':12},indent=2))
