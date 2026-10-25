"""Minimal synthetic checks; writes only to an ignored .smoke directory.

All text and shapes here are generic test fixtures, not study material. The
test deliberately leaves manual review pending. It never invents a review PASS.
"""
from pathlib import Path
import json
import fitz
from docx import Document
from docx.shared import Mm
from reportlab.pdfgen import canvas
from authoring_checks import sha256,json_write,binary_package_hashes,docx_text_parts
from assemble_pdf import assemble
from pdf_qa import audit,review_template,validate_review
from replace_references import apply
from format_docx import format_document

def rejected(call):
    try:call()
    except ValueError:return True
    raise AssertionError('Expected a rejected unsafe/stale operation')

def main():
    output=Path(__file__).resolve().parent/'.smoke'
    output.mkdir(exist_ok=True)
    (output/'page.md').write_text('# Example document\n\nThis paragraph contains a generic reference to Appendix A.\n\n## Example section\n\nThe fixture label remains visible.\n',encoding='utf-8')
    figure=output/'figure.pdf'
    can=canvas.Canvas(str(figure),pagesize=(180*72/25.4,75*72/25.4))
    can.setFont('Helvetica',10);can.drawString(20,180,'Example figure');can.rect(20,30,140,100);can.save()
    layout={'metadata':{'title':'Example document'},'minimum_font_pt':8,'pages':[{'blocks':[{'type':'markdown','path':'page.md'},{'type':'figure_pdf','path':'figure.pdf'}]}]}
    json_write(output/'layout.json',layout)
    assemble(output/'layout.json',output/'assembled.pdf',output/'build.json')
    report=audit(output/'assembled.pdf',output/'render',expected_pages=1)
    assert not report['errors'] and report['status']=='AUTOMATED_CHECKS_PASS_PENDING_MANUAL_REVIEW'
    review_template(report,output/'review_pending.json')
    checks={'pdf_assembly_and_automated_checks':True,'pending_review_rejected':rejected(lambda:validate_review(output/'assembled.pdf',output/'render/AUTOMATED_PDF_QA.json',output/'review_pending.json',output/'must_not_exist.json'))}
    plan={'rules':[{'id':'label','old':'Appendix A','new':'Appendix B','expected':1}]}
    json_write(output/'plan.json',plan)
    apply(output/'page.md',output/'revised.md',output/'plan.json',output/'md_audit.json')
    checks['markdown_literal_reference_replacement']=True
    # A figure-page raster becomes an independent synthetic DOCX media fixture.
    with fitz.open(figure) as doc:doc[0].get_pixmap(dpi=96).save(output/'toy.png')
    document=Document();document.add_paragraph('Example document',style='Title')
    paragraph=document.add_paragraph();paragraph.add_run('Appendix ');paragraph.add_run('A').bold=True
    paragraph.add_run(' retains value 12.5.')
    document.add_picture(str(output/'toy.png'),width=Mm(75));document.save(output/'input.docx')
    plan['rules'][0]['scope']={'paragraph_indices':[1]}
    plan['input_sha256']=sha256(output/'input.docx');json_write(output/'docx_plan.json',plan)
    apply(output/'input.docx',output/'revised.docx',output/'docx_plan.json',output/'docx_audit.json')
    assert binary_package_hashes(output/'input.docx')==binary_package_hashes(output/'revised.docx')
    checks['split_run_replacement_and_image_bytes_preserved']=True
    bad={'rules':[{'id':'value','old':'12.5','new':'99.0','expected':1}]}
    json_write(output/'bad_numeric_plan.json',bad)
    checks['changed_scientific_numeric_token_rejected']=rejected(lambda:apply(output/'input.docx',output/'invalid.docx',output/'bad_numeric_plan.json',output/'invalid_audit.json'))
    styles={'section_mm':{'left_margin':20,'right_margin':20},'styles':{'Normal':{'font_name':'Arial','font_pt':11},'Title':{'font_name':'Arial','font_pt':16}}}
    json_write(output/'styles.json',styles)
    format_document(output/'revised.docx',output/'formatted.docx',output/'styles.json',output/'format_audit.json')
    checks['formatting_preserves_all_visible_text_and_media_bytes']=docx_text_parts(output/'revised.docx')==docx_text_parts(output/'formatted.docx')
    too_tall={'pages':[{'blocks':[{'type':'figure_pdf','path':'figure.pdf','width_mm':200}]}]}
    json_write(output/'overflow.json',too_tall)
    checks['oversized_artwork_rejected']=rejected(lambda:assemble(output/'overflow.json',output/'overflow.pdf',output/'overflow_audit.json'))
    assert all(checks.values())
    json_write(output/'SMOKE_RECEIPT.json',{'status':'SYNTHETIC_AUTOMATED_CHECKS_PASS','checks':checks,'manual_visual_review':'PENDING; not manufactured by the smoke test','word_export':'OPTIONAL_NOT_EXERCISED'})
    print(json.dumps({'status':'SYNTHETIC_AUTOMATED_CHECKS_PASS','checks':checks,'manual_visual_review':'PENDING'}))

if __name__=='__main__':main()
