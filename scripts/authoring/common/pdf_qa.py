"""Separate automated PDF checks/renders from external manual review records."""
import argparse
import json
from pathlib import Path
import fitz
from authoring_checks import sha256, json_write

def audit(pdf, output_dir, minimum_font=8, expected_pages=None, dpis=(96,144)):
    pdf, output_dir = Path(pdf).resolve(), Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    errors, pages = [], []
    with fitz.open(pdf) as doc:
        if expected_pages is not None and len(doc) != expected_pages:
            errors.append('Unexpected page count')
        for index, page in enumerate(doc, 1):
            spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
            outside = [s['bbox'] for s in spans if s['bbox'][0] < -.2 or s['bbox'][1] < -.2 or s['bbox'][2] > page.rect.width+.2 or s['bbox'][3] > page.rect.height+.2]
            small = [s['size'] for s in spans if s['size']+.01 < minimum_font]
            if outside:
                errors.append(f'Page {index}: text outside page boundary')
            if small:
                errors.append(f'Page {index}: text below minimum font')
            blank = not spans and not page.get_images() and not page.get_drawings()
            if blank:
                errors.append(f'Page {index}: blank page')
            views=[]
            for dpi in dpis:
                relative=Path(f'{dpi}dpi')/f'page-{index:03d}.png'
                path=output_dir/relative
                path.parent.mkdir(exist_ok=True)
                page.get_pixmap(dpi=dpi, alpha=False).save(path)
                views.append({'dpi': dpi, 'file': relative.as_posix(), 'sha256': sha256(path)})
            pages.append({'page': index, 'width_pt': page.rect.width, 'height_pt': page.rect.height,
                          'text_span_count': len(spans), 'minimum_text_pt': min((s['size'] for s in spans),default=None),
                          'text_outside_page_bboxes': outside, 'below_minimum_font_sizes': small, 'blank': blank, 'renders': views})
    report={'status': 'AUTOMATED_CHECKS_FAIL' if errors else 'AUTOMATED_CHECKS_PASS_PENDING_MANUAL_REVIEW',
            'pdf_basename': pdf.name, 'pdf_sha256': sha256(pdf), 'page_count':len(pages), 'minimum_font_pt_required':minimum_font,
            'errors': errors, 'pages':pages, 'manual_review': 'NOT_PERFORMED_BY_THIS_COMMAND',
            'limits':'No automated claim about visual overlap, glyph quality, figure content or scientific correctness.'}
    json_write(output_dir/'AUTOMATED_PDF_QA.json', report)
    return report

def review_template(report, output):
    json_write(output, {'pdf_sha256':report['pdf_sha256'], 'reviewer':None, 'reviewed_utc':None,
                        'attestation':None, 'pages':[{'page':p['page'],'dpi':v['dpi'],'render_sha256':v['sha256'],'status':'PENDING','notes':''}
                        for p in report['pages'] for v in p['renders']]})

def validate_review(pdf, automated_path, manual_path, output):
    automatic=json.loads(Path(automated_path).read_text(encoding='utf-8'))
    manual=json.loads(Path(manual_path).read_text(encoding='utf-8'))
    if automatic['status']!='AUTOMATED_CHECKS_PASS_PENDING_MANUAL_REVIEW':
        raise ValueError('Automated checks have not passed')
    if automatic['pdf_sha256']!=sha256(pdf) or manual['pdf_sha256']!=automatic['pdf_sha256']:
        raise ValueError('Review is not bound to the current PDF hash')
    if not manual.get('reviewer') or not manual.get('reviewed_utc') or manual.get('attestation')!='I inspected these rendered page images.':
        raise ValueError('A named, dated external visual-review attestation is required')
    with fitz.open(pdf) as doc:
        if automatic['page_count']!=len(doc) or [p['page'] for p in automatic['pages']]!=list(range(1,len(doc)+1)):
            raise ValueError('Automated report does not cover every current PDF page')
    if any(not p['renders'] for p in automatic['pages']):
        raise ValueError('Every page must have at least one registered rendered view')
    required={(p['page'],v['dpi']):v for p in automatic['pages'] for v in p['renders']}
    supplied={(p['page'],p['dpi']):p for p in manual['pages']}
    if len(supplied)!=len(manual['pages']) or set(supplied)!=set(required):
        raise ValueError('Review must cover each page at every recorded DPI exactly once')
    for key,view in required.items():
        actual=Path(automated_path).parent/view['file']
        if sha256(actual)!=view['sha256'] or supplied[key].get('render_sha256')!=view['sha256'] or supplied[key].get('status')!='PASS':
            raise ValueError('Missing, stale or failing manual page review: '+str(key))
    json_write(output, {'status':'EXTERNAL_MANUAL_REVIEW_REGISTERED_PASS', 'pdf_sha256':automatic['pdf_sha256'],
                        'automated_report_sha256':sha256(automated_path),'manual_record_sha256':sha256(manual_path),
                        'reviewer':manual['reviewer'],'reviewed_utc':manual['reviewed_utc'],
                        'reviewed_views':len(required),'note':'This validates the external record and immutable render hashes; it does not perform the visual review.'})

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    a=sub.add_parser('audit-render'); a.add_argument('--pdf',required=True,type=Path); a.add_argument('--output-dir',required=True,type=Path)
    a.add_argument('--minimum-font',type=float,default=8); a.add_argument('--expected-pages',type=int); a.add_argument('--dpis',type=int,nargs='+',default=[96,144])
    a.add_argument('--review-template',type=Path)
    v=sub.add_parser('validate-review'); v.add_argument('--pdf',required=True,type=Path);v.add_argument('--automated-report',required=True,type=Path)
    v.add_argument('--manual-review',required=True,type=Path);v.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    if args.command=='audit-render':
        report=audit(args.pdf,args.output_dir,args.minimum_font,args.expected_pages,args.dpis)
        if args.review_template: review_template(report,args.review_template)
        if report['errors']: raise SystemExit(1)
    else: validate_review(args.pdf,args.automated_report,args.manual_review,args.output)

if __name__=='__main__':
    main()
