"""Read private Ct/OD and protocol inputs; write private, provenance-gated extracts."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import openpyxl
from sequence_parser import docx_paragraphs,parse_sequence_pairs
STATUS='INTERNAL_REVIEW_SOURCE_UNRESOLVED'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def extract(workbook:Path,protocol:Path,project:Path):
    if not workbook.is_file() or not protocol.is_file():raise FileNotFoundError('Private workbook and protocol are required')
    initial=[sha(workbook),sha(protocol)]
    out=project/'03_tables/experimental/internal/source_unresolved';out.mkdir(parents=True,exist_ok=True)
    audit=project/'00_inventory';audit.mkdir(parents=True,exist_ok=True)
    wb=openpyxl.load_workbook(workbook,read_only=True,data_only=False)
    required={'siRNA_qPCR_Raw','Program_qPCR_Raw','CCK8_Raw'}
    if set(wb.sheetnames)!=required:raise ValueError('Unexpected raw observation worksheet contract')
    raw={};sheets=[]
    for ws in wb:
        values=list(ws.values);headers=list(values[0])
        if len(headers)!=len(set(headers)) or any(h is None for h in headers):raise ValueError('Invalid raw headers')
        if any(isinstance(v,str) and v.startswith('=') for row in values for v in row):raise ValueError('Raw observations must contain source values, not formulas')
        rows=[]
        for i,values_row in enumerate(values[1:],2):
            row=dict(zip(headers,values_row));row.update(Source_sheet=ws.title,Source_excel_row=i);rows.append(row)
        raw[ws.title]=rows;sheets.append({'sheet':ws.title,'rows':len(rows),'columns':len(headers)})
        with (out/(ws.title+'.tsv')).open('w',encoding='utf-8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
    wb.close()
    (out/'RAW_EXTRACT.json').write_text(json.dumps(raw,ensure_ascii=False,indent=2),encoding='utf-8')
    paragraphs=docx_paragraphs(protocol);pairs=parse_sequence_pairs(paragraphs)
    if not pairs:raise ValueError('No complete same-paragraph oligo pairs found')
    with (out/'SUPPLIED_SEQUENCES_EXACT.tsv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(pairs[0]),delimiter='\t');writer.writeheader();writer.writerows(pairs)
    # The supplied protocol remains a private input. No transcription is committed.
    assert initial==[sha(workbook),sha(protocol)],'An original input changed'
    receipt={'source_status':STATUS,'analytical_status':'READ_ONLY_EXTRACTION_DONE','sheets':sheets,'sequence_pairs':len(pairs),'input_sha256':{'workbook':initial[0],'protocol':initial[1]},'experimental_authenticity_verified':False,'instrument_records_verified':False,'submission_eligible':False}
    (audit/'EXPERIMENT_INPUT_STATUS.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    return receipt

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workbook',type=Path,required=True);parser.add_argument('--protocol',type=Path,required=True);parser.add_argument('--project',type=Path,required=True)
    a=parser.parse_args();print(json.dumps(extract(a.workbook.resolve(),a.protocol.resolve(),a.project.resolve()),indent=2))
if __name__=='__main__':main()
