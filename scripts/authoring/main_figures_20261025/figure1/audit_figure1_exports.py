"""Read-only structural checks and PDF raster QA for the original Figure 1."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import fitz
from PIL import Image

ap=argparse.ArgumentParser();ap.add_argument('--figure-root',type=Path,required=True)
ap.add_argument('--audit-root',type=Path,required=True);args=ap.parse_args()
args.audit_root.mkdir(parents=True,exist_ok=True)
stem=args.figure_root/'Figure1_Study_Design_Workflow'
doc=fitz.open(str(stem)+'.pdf');page=doc[0]
spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans'] if s['text'].strip()]
pdftext=page.get_text()
page.get_pixmap(matrix=fitz.Matrix(2.5,2.5),alpha=False).save(args.audit_root/'Figure1_PDF_actual_render.png')
svg=ET.parse(str(stem)+'.svg').getroot();ns={'s':'http://www.w3.org/2000/svg'}
svgtexts=svg.findall('.//s:text',ns)
measurements={}
for suffix in ['png','tiff']:
    with Image.open(str(stem)+'.'+suffix) as im:
        measurements[suffix]={'pixels':list(im.size),'dpi':[float(x) for x in im.info.get('dpi',[])], 'mode':im.mode}
required=['GSE230808','GSE165722','GSE287499','CP-H097Y','6 eligible AF','UCell','AUCell',
          'CellChat','SigXTalk','ITGB1','HIF1A','PLOD2','schematic','siITGB1-1','siITGB1-2','CCK-8']
checks={
 'one_pdf_page':len(doc)==1,
 'width_180mm':abs(page.rect.width*25.4/72-180)<.02,
 'height_220mm':abs(page.rect.height*25.4/72-220)<.02,
 'minimum_pdf_font_8pt':min(s['size'] for s in spans)>=7.99,
 'pdf_text_selectable':len(pdftext)>800,
 'svg_editable_text':len(svgtexts)>50,
 'svg_contains_no_embedded_raster':len(svg.findall('.//s:image',ns))==0,
 'required_labels':all(x in pdftext for x in required),
 'no_result_p_values':not bool(re.search(r'P\s*[<=>]\s*0\.',pdftext)),
 'all_rasters_at_least_300dpi':all(len(m['dpi'])==2 and min(m['dpi'])>=299 for m in measurements.values())
}
receipt={'checks':checks,'automated_status':'PASS' if all(checks.values()) else 'FAIL',
         'pdf_size_mm':[page.rect.width*25.4/72,page.rect.height*25.4/72],
         'minimum_pdf_font_pt':min(s['size'] for s in spans),'pdf_fonts':sorted(set(s['font'] for s in spans)),
         'svg_text_count':len(svgtexts),'raster_metadata':measurements,
         'files':{f.name:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size}
                  for f in sorted(args.figure_root.glob('Figure1_Study_Design_Workflow.*'))},
         'manual_visual_review':'PENDING; actual PDF render and full figure must be inspected separately'}
(args.audit_root/'Figure1_EXPORT_STRUCTURAL_QA.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps({k:receipt[k] for k in ['checks','automated_status','minimum_pdf_font_pt']},indent=2))
