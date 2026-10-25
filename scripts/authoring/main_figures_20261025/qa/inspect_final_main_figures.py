#!/usr/bin/env python3
"""Read-only structural checks and Python PDF/TIFF previews for six final figures.

This command generates only QA previews and objective measurements. Manual
review is recorded separately, after the generated pixels have actually been
opened and assessed.
"""
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime,timezone
import xml.etree.ElementTree as ET

import fitz
from PIL import Image, ImageDraw, ImageFont


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--project',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    source=a.project/'figures';out=a.output;out.mkdir(parents=True,exist_ok=True);records=[];pictures=[]
    for n in range(1,7):
        stem='Figure1_Study_Design_Workflow' if n==1 else 'Figure'+str(n)
        files={ext:source/(stem+'.'+ext) for ext in ['pdf','svg','png','tiff']}
        record={'figure':'Figure'+str(n),'files':{ext:{'path':p.relative_to(a.project).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for ext,p in files.items()},'manual_visual_status':'PENDING_ACTUAL_REVIEW'}
        with fitz.open(files['pdf']) as doc:
            if len(doc)!=1:raise ValueError(stem+' must have exactly one PDF page.')
            page=doc[0];spans=[s for b in page.get_text('dict')['blocks'] if 'lines'in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
            outside=[{'text':s['text'],'bbox_pt':s['bbox']} for s in spans if s['bbox'][0]<-.2 or s['bbox'][1]<-.2 or s['bbox'][2]>page.rect.width+.2 or s['bbox'][3]>page.rect.height+.2]
            record['pdf']={'pages':len(doc),'width_mm':page.rect.width*25.4/72,'height_mm':page.rect.height*25.4/72,'text_spans':len(spans),'text_characters':len(page.get_text()),'font_names':sorted(set(s['font'] for s in spans)),'minimum_text_pt':min(s['size'] for s in spans),'text_below_8pt':[{'text':s['text'],'size':s['size']} for s in spans if s['size']<7.99],'outside_page_text':outside,'fonts':[{ 'xref':f[0],'extension':f[1],'type':f[2],'basefont':f[3],'encoding':f[5]} for f in page.get_fonts()],'embedded_raster_images':len(page.get_images())}
            for dpi,label in [(100,'100pct'),(200,'200pct')]:
                dest=out/(f'Figure{n}_{label}_PDF.png');page.get_pixmap(matrix=fitz.Matrix(dpi/72,dpi/72),alpha=False).save(dest)
                record[label+'_preview']={'path':dest.name,'sha256':sha(dest),'dpi':dpi}
            pictures.append(out/(f'Figure{n}_100pct_PDF.png'))
        with Image.open(files['tiff']) as image:
            tiff_dpi=image.info.get('dpi');tiff_dpi=[float(v) for v in tiff_dpi] if tiff_dpi else None
            record['tiff']={'dimensions_px':list(image.size),'mode':image.mode,'compression':image.info.get('compression'),'dpi':tiff_dpi,'actual_width_mm':image.width/tiff_dpi[0]*25.4 if tiff_dpi else None,'actual_height_mm':image.height/tiff_dpi[1]*25.4 if tiff_dpi else None}
            copy=image.copy().convert('RGB');copy.thumbnail((709,1100));dest=out/(f'Figure{n}_TIFF_check.png');copy.save(dest);record['tiff_check_preview']={'path':dest.name,'sha256':sha(dest)}
        with Image.open(files['png']) as image:
            record['png']={'dimensions_px':list(image.size),'mode':image.mode,'dpi':list(image.info['dpi']) if image.info.get('dpi') else None}
        svg=ET.parse(files['svg']).getroot();record['svg']={'text_elements':len(svg.findall('.//{http://www.w3.org/2000/svg}text')),'path_elements':len(svg.findall('.//{http://www.w3.org/2000/svg}path')),'embedded_image_elements':len(svg.findall('.//{http://www.w3.org/2000/svg}image')),'width':svg.get('width'),'height':svg.get('height')}
        issues=[]
        if abs(record['pdf']['width_mm']-180)>.02:issues.append('PDF width differs from 180 mm')
        if record['pdf']['outside_page_text']:issues.append('text outside PDF page')
        if record['pdf']['text_below_8pt']:issues.append('some PDF text below 8 pt')
        if not record['svg']['text_elements']:issues.append('SVG text is not editable')
        if not record['tiff']['dpi'] or min(record['tiff']['dpi'])<299.5:issues.append('TIFF metadata DPI is missing or below 300')
        for ext,p in files.items():
            if sha(p)!=record['files'][ext]['sha256']:raise RuntimeError('A final figure changed during read-only QA: '+str(p))
        record['automated_status']='PASS' if not issues else 'REVIEW';record['automated_issues']=issues;records.append(record)
    # Contact sheet is an assessment aid, not a modification of any main figure.
    thumb_width=580;header=46;pad=24;loaded=[]
    for p in pictures:
        with Image.open(p) as image:
            im=image.convert('RGB');im=im.resize((thumb_width,round(im.height*thumb_width/im.width)),Image.Resampling.LANCZOS);loaded.append(im)
    row_heights=[max(im.height for im in loaded[:3]),max(im.height for im in loaded[3:])]
    sheet=Image.new('RGB',(3*(thumb_width+pad)+pad,sum(row_heights)+2*header+3*pad),'#E6E8EA');draw=ImageDraw.Draw(sheet)
    try:font=ImageFont.truetype('DejaVuSans.ttf',22)
    except OSError:font=ImageFont.load_default()
    y=pad
    for row in range(2):
        for col in range(3):
            i=row*3+col;x=pad+col*(thumb_width+pad);draw.text((x,y+6),'Figure '+str(i+1),fill='#202020',font=font);sheet.paste(loaded[i],(x,y+header))
        y+=header+row_heights[row]+pad
    contact=out/'MAIN_FIGURES_CONTACT_SHEET.png';sheet.save(contact)
    receipt={'generated_utc':datetime.now(timezone.utc).isoformat(),'workflow':'Python PyMuPDF PDF rendering + Pillow TIFF metadata/QA thumbnails/contact sheet','selected_plotting_backend':'Python','final_artwork_modified':False,'figures':records,'contact_sheet':{'path':contact.name,'sha256':sha(contact)},'manual_review_note':'No manual PASS is issued by this automated script.'}
    dest=out/'MAIN_FIGURES_STRUCTURAL_QA.json';dest.write_text(json.dumps(receipt,indent=2),encoding='utf8');print(json.dumps({'figures':[{k:r[k] for k in ['figure','automated_status','automated_issues']} for r in records],'output':str(dest)},indent=2))


if __name__=='__main__':main()
