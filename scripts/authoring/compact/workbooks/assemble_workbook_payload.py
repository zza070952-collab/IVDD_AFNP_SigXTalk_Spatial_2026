import json,sys
from pathlib import Path
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--revision',type=Path,required=True);ap.add_argument('--work-dir',type=Path,required=True);ap.add_argument('--experimental',action='store_true');a=ap.parse_args()
C=a.revision.resolve();W=a.work_dir.resolve();W.mkdir(parents=True,exist_ok=True)
e=json.loads((C/'source_tables/EXPERIMENT_PAYLOAD.json').read_text(encoding='utf-8'))
sheets=[]
for pattern in ([] if '--experimental' in sys.argv else ['human/*PAYLOAD*.json','communication_spatial/*PAYLOAD*.json']):
 for p in C.glob(pattern):
  d=json.loads(p.read_text(encoding='utf-8'));ss=d.get('sheets',d.get('supplementary_sheets',[]))
  sheets.extend(ss)
names={s['name'] for s in sheets}
if sheets:
 assert len(names)==len(sheets),names
 assert len(sheets)==5,names
sheets.append(e['S6_sheet'])
out={'output':str((C/'INTERNAL_REVIEW').resolve()),'qa':str((C/'audit/workbook_render').resolve()),'books':[{'filename':'Additional_file_2_Supplementary_Tables.xlsx','sheets':sheets},{'filename':'Additional_file_3_Experimental_Source_Data.xlsx','sheets':e['experiment_sheets']}]}
(W/'COMPACT_WORKBOOK_PAYLOAD.json').write_text(json.dumps(out,ensure_ascii=False,allow_nan=False),encoding='utf-8')
print([(b['filename'],len(b['sheets'])) for b in out['books']])
