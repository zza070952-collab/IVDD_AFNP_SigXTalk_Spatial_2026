from pathlib import Path
import csv, datetime, hashlib, json, math, re, zipfile
from xml.etree import ElementTree as ET

import argparse
ap=argparse.ArgumentParser()
ap.add_argument('--payload',type=Path,required=True)
a=ap.parse_args()
payload=json.loads(a.payload.read_text(encoding='utf-8'))
ROOT=Path(payload['root']); QA=ROOT/'08_review'/'workbook_qa'
N={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def colnum(a):
    n=0
    for c in a:n=n*26+ord(c)-64
    return n-1
def equal(a,b):
    if isinstance(a,bool):return isinstance(b,bool) and a==b
    if isinstance(a,(float,int)):
        return isinstance(b,(float,int)) and not isinstance(b,bool) and math.isclose(a,b,rel_tol=2e-13,abs_tol=2e-13)
    return a==b
reports=[]
for b in payload['books']:
    p=ROOT/b['file']
    if not p.exists():raise RuntimeError('Workbook not exported: '+str(p))
    report={'file':b['file'],'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size,'sheets':[],'errors':[],'formula_count':0,'merged_regions':0,'hidden_sheets':0}
    text_cells=[]
    with zipfile.ZipFile(p) as z:
        assert z.testzip() is None
        shared=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            shared=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',N)]
        wb=ET.fromstring(z.read('xl/workbook.xml'))
        rel={x.attrib['Id']:x.attrib['Target'] for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        expected={s['name']:s for s in b['sheets']}
        for sh in wb.findall('s:sheets/s:sheet',N):
            name=sh.attrib['name']; rid=sh.attrib['{'+N['r']+'}id']; target=rel[rid]
            target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
            xml=ET.fromstring(z.read(target)); matrix={};native_types={}
            for c in xml.findall('s:sheetData/s:row/s:c',N):
                address=c.attrib['r']; t=c.attrib.get('t','n'); v=c.find('s:v',N)
                if t=='s':value=shared[int(v.text)] if v is not None else None
                elif t=='inlineStr':value=''.join(c.find('s:is',N).itertext())
                elif t=='str':value=v.text if v is not None else None
                elif t=='b':value=v is not None and v.text=='1'
                elif t=='d':
                    dt=datetime.datetime.fromisoformat(v.text.replace('Z','+00:00'))
                    if dt.tzinfo is not None:dt=dt.astimezone(datetime.timezone.utc).replace(tzinfo=None)
                    value=(dt-datetime.datetime(1899,12,30)).total_seconds()/86400
                elif t=='e':value=v.text if v is not None else None; report['errors'].append(name+'!'+address+' unexpectedExcelError')
                elif v is not None and v.text is not None:
                    try:value=float(v.text);value=int(value) if value.is_integer() else value
                    except ValueError:value=v.text
                else:value=None
                matrix[address]=value;native_types[address]=t
                if isinstance(value,str):text_cells.append(value)
                if c.find('s:f',N) is not None:report['formula_count']+=1
            merges=xml.find('s:mergeCells',N)
            if merges is not None:report['merged_regions']+=len(merges)
            if sh.attrib.get('state','visible')!='visible':report['hidden_sheets']+=1
            pane=xml.find('s:sheetViews/s:sheetView/s:pane',N)
            freeze={'xSplit':pane.attrib.get('xSplit') if pane is not None else None,'ySplit':pane.attrib.get('ySplit') if pane is not None else None}
            entry={'sheet':name,'static_cells':len(matrix),'freeze':freeze,'comparison':'METADATA_SHEET'}
            if name in expected:
                src=expected[name]; rows=[src['headers']]+src['rows'];fail=[]
                for i,row in enumerate(rows,1):
                    for j,value in enumerate(row):
                        if i>1 and value is not None and src['types'][j]=='datetime':
                            dt=datetime.datetime.fromisoformat(value.replace('Z','+00:00'))
                            if dt.tzinfo is not None:dt=dt.astimezone(datetime.timezone.utc).replace(tzinfo=None)
                            value=(dt-datetime.datetime(1899,12,30)).total_seconds()/86400
                        n=j+1;a=''
                        while n:a=chr(65+(n-1)%26)+a;n=(n-1)//26
                        address=a+str(i)
                        if not equal(value,matrix.get(address)):fail.append({'cell':address,'expected':value,'actual':matrix.get(address)})
                extras=[a for a,v in matrix.items() if v is not None and (int(re.search(r'\d+',a).group())>len(rows) or colnum(re.match(r'[A-Z]+',a).group())>=len(src['headers']))]
                entry.update(expected_data_rows=len(src['rows']),columns=len(src['headers']),mismatches=len(fail),extra_cells=len(extras),source_sha256=src['source_sha256'],comparison='PASS' if not fail and not extras else 'FAIL')
                if fail:entry['first_mismatches']=fail[:5];report['errors'].append(name+' source-matrix mismatch')
                if extras:report['errors'].append(name+' unexpectedextra cells')
                if freeze!={'xSplit':'1','ySplit':'1'}:report['errors'].append(name+' missingfreeze rows/columns')
            report['sheets'].append(entry)
        for fn in z.namelist():
            if fn.startswith('xl/tables/') and fn.endswith('.xml'):
                table=ET.fromstring(z.read(fn))
                if table.find('s:autoFilter',N) is None:report['errors'].append(fn+' missingautofilter')
        report['native_table_count']=len([fn for fn in z.namelist() if fn.startswith('xl/tables/') and fn.endswith('.xml')])
        report['all_sheets_have_native_tables']=report['native_table_count']==len(report['sheets'])
        if not report['all_sheets_have_native_tables']:report['errors'].append('Missing native data tables')
        if not b['internal']:
            private=[s for s in text_cells if re.search(r'[A-Za-z]:[/\\]Users[/\\]|/(?:home|Users)/[^/\s]+/|/mnt/data/',s,re.I)]
            report['private_path_cells']=len(private)
            if private:report['errors'].append('privatepaths in publicworkbook')
    if report['formula_count']:report['errors'].append('Unexpectedformula in staticdata export')
    if report['merged_regions']:report['errors'].append('Mergedresearch data cells')
    if report['hidden_sheets']:report['errors'].append('Hiddensheets')
    report['status']='PASS' if not report['errors'] else 'FAIL'
    reports.append(report)
(QA/'WORKBOOK_XLSX_NATIVE_VALIDATION.json').write_text(json.dumps(reports,indent=2,ensure_ascii=False),encoding='utf-8')
with (QA/'WORKBOOK_VALIDATION_SUMMARY.tsv').open('w',encoding='utf-8',newline='') as f:
    cols=['file','status','bytes','sha256','sheet_count','data_rows','native_table_count','formula_count','merged_regions','hidden_sheets']
    w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
    for r in reports:w.writerow({**{k:r[k] for k in ['file','status','bytes','sha256','native_table_count','formula_count','merged_regions','hidden_sheets']},'sheet_count':len(r['sheets']),'data_rows':sum(s.get('expected_data_rows',0) for s in r['sheets'])})
print(json.dumps([{'file':r['file'],'status':r['status'],'errors':r['errors']} for r in reports],indent=2))
if any(r['status']!='PASS' for r in reports):raise SystemExit(1)
