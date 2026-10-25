from pathlib import Path
import argparse
import json,hashlib,datetime
from pdf2image import convert_from_path,pdfinfo_from_path
from PIL import Image,ImageDraw,ImageFont
from pypdf import PdfReader
parser=argparse.ArgumentParser(description='Render an external private PDF for actual page review.')
parser.add_argument('--pdf',type=Path,required=True)
parser.add_argument('--output-root',type=Path,required=True)
parser.add_argument('--poppler-path',type=Path,help='Directory containing pdftoppm and pdfinfo; otherwise use PATH')
parser.add_argument('--include-text',action='store_true',help='Also export private full page text; never commit that output')
args=parser.parse_args()
P=args.pdf.resolve();Q=args.output_root.resolve();O=Q
if P==Q or P in Q.parents or Q in P.parents:
    raise ValueError('Use a separate QA directory, not the input PDF or its ancestor.')
checkout=next((p for p in Path(__file__).resolve().parents if (p/'.git').exists()),None)
if checkout is not None and (Q==checkout or checkout in Q.parents):
    raise ValueError('Private document render outputs must be outside this publication checkout.')
Q.mkdir(parents=True,exist_ok=True)
poppler=str(args.poppler_path) if args.poppler_path else None
before=hashlib.sha256(P.read_bytes()).hexdigest()
info=pdfinfo_from_path(str(P),poppler_path=poppler);n=info['Pages']
pages=[]
reader=PdfReader(P)
for dpi in [96,144]:
    dest=Q/(str(dpi)+'dpi');dest.mkdir(exist_ok=True)
    result=convert_from_path(str(P),dpi=dpi,output_folder=str(dest),fmt='png',paths_only=True,poppler_path=poppler,thread_count=2)
    assert len(result)==n
    for i,p in enumerate(result,1):
        f=dest/('page-'+str(i).zfill(3)+'.png');Path(p).replace(f)
        with Image.open(f) as im:sz=list(im.size)
        pages.append({'page':i,'dpi':dpi,'path':f.relative_to(O).as_posix(),'pixels':sz,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
thumbs=[]
for i in range(1,n+1):
    with Image.open(Q/'96dpi'/('page-'+str(i).zfill(3)+'.png')) as im:
        im=im.convert('RGB');im.thumbnail((320,455));tile=Image.new('RGB',(340,490),'#e9edef');tile.paste(im,((340-im.width)//2,25));ImageDraw.Draw(tile).text((12,6),'Page '+str(i),fill='black');thumbs.append(tile)
sheet=Image.new('RGB',(340*4,490*((n+3)//4)),'white')
for i,im in enumerate(thumbs):sheet.paste(im,((i%4)*340,(i//4)*490))
sheet.save(Q/'MANUSCRIPT_PAGE_CONTACT_SHEET.png')
textpages=[{'page':i+1,'characters':len(p.extract_text() or ''),'text':p.extract_text() or ''} for i,p in enumerate(reader.pages)]
if args.include_text:
    (Q/'PDF_PAGE_TEXT.json').write_text(json.dumps(textpages,indent=2),encoding='utf8')
receipt={'status':'RENDERED_PENDING_ACTUAL_PAGE_VISUAL_REVIEW','pdf_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'pages':n,'engine':'External PDF; Poppler via pdf2image','renders':pages,'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
assert before==hashlib.sha256(P.read_bytes()).hexdigest(), 'Input PDF changed during rendering'
(Q/'PDF_RENDER_MANIFEST.json').write_text(json.dumps(receipt,indent=2),encoding='utf8')
print(json.dumps({'pages':n,'page_rasters':len(pages),'status':receipt['status']}))
