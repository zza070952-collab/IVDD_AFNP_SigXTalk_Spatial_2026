"""Optional external style normalization and separate Word-native PDF export.

Formatting preserves visible text and embedded media bytes. It does not insert
body text, captions or replacement content. Export requires Windows, Microsoft
Word and pywin32, and opens a separate invisible Word instance read-only.
"""
import argparse
import json
import os
from pathlib import Path
from docx import Document
from docx.shared import Pt, Mm, RGBColor
from docx.oxml.ns import qn
from authoring_checks import sha256,json_write,different_output,check_docx_preserved

def format_document(source,output,style_path,receipt):
    for path in [source,style_path]:
        different_output(path,output);different_output(path,receipt)
    different_output(output,receipt)
    spec=json.loads(Path(style_path).read_text(encoding='utf-8'))
    document=Document(source)
    for section in document.sections:
        for name,value in spec.get('section_mm',{}).items():
            if name not in ['top_margin','bottom_margin','left_margin','right_margin','page_width','page_height']:
                raise ValueError('Unsupported section property: '+name)
            setattr(section,name,Mm(value))
    for name,values in spec.get('styles',{}).items():
        style=document.styles[name]
        if 'font_name' in values:style.font.name=values['font_name']
        if 'font_pt' in values:style.font.size=Pt(values['font_pt'])
        if 'bold' in values:style.font.bold=values['bold']
        style.font.color.rgb=RGBColor(0,0,0)
        style.font.underline=False
        for key in ['space_before_pt','space_after_pt']:
            if key in values:setattr(style.paragraph_format,key.removesuffix('_pt'),Pt(values[key]))
        if 'line_spacing' in values:style.paragraph_format.line_spacing=values['line_spacing']
    heading_names={'Title','Subtitle'}|{s.name for s in document.styles if s.name.startswith('Heading')}
    for name in heading_names:
        if name not in document.styles:continue
        style=document.styles[name]
        style.font.color.rgb=RGBColor(0,0,0);style.font.underline=False
        for border in list(style.element.iter(qn('w:pBdr'))):border.getparent().remove(border)
    for paragraph in document.paragraphs:
        if paragraph.style.name in heading_names:
            for border in list(paragraph._p.iter(qn('w:pBdr'))):border.getparent().remove(border)
            for run in paragraph.runs:run.font.color.rgb=RGBColor(0,0,0);run.font.underline=False
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    document.save(output)
    check_docx_preserved(source,output,require_text=True)
    json_write(receipt,{'status':'FORMATTED_PENDING_RENDER_AND_VISUAL_REVIEW','input_sha256':sha256(source),
                        'output_sha256':sha256(output),'style_sha256':sha256(style_path),
                        'visible_text':'UNCHANGED','embedded_media_bytes':'UNCHANGED','visual_review':'NOT_PERFORMED_BY_THIS_COMMAND'})

def export_word(source,output,receipt):
    if os.name!='nt':raise RuntimeError('Native Word export requires Windows')
    import win32com.client
    source=Path(source).resolve();output=Path(output).resolve();different_output(source,output);different_output(source,receipt);different_output(output,receipt)
    original=sha256(source)
    output.parent.mkdir(parents=True,exist_ok=True)
    word=win32com.client.DispatchEx('Word.Application')
    word.Visible=False;word.DisplayAlerts=0
    document=None
    try:
        document=word.Documents.Open(str(source),ReadOnly=True,AddToRecentFiles=False,Visible=False)
        document.ExportAsFixedFormat(str(output),17)
    finally:
        if document is not None:document.Close(False)
        word.Quit()
    if sha256(source)!=original:raise ValueError('Input DOCX changed during export')
    if not output.is_file():raise ValueError('Word did not produce a PDF')
    json_write(receipt,{'status':'WORD_EXPORTED_PENDING_PDF_QA','input_sha256':original,'output_sha256':sha256(output),
                        'input_docx_unchanged':True,'visual_review':'NOT_PERFORMED_BY_THIS_COMMAND'})

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    f=sub.add_parser('format');f.add_argument('--styles',required=True,type=Path)
    w=sub.add_parser('word-pdf')
    for p in [f,w]:
        p.add_argument('--input',required=True,type=Path);p.add_argument('--output',required=True,type=Path);p.add_argument('--receipt',required=True,type=Path)
    args=parser.parse_args()
    if args.command=='format':format_document(args.input,args.output,args.styles,args.receipt)
    else:export_word(args.input,args.output,args.receipt)

if __name__=='__main__':main()
