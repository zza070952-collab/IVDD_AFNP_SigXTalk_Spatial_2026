"""Assemble external Markdown text and existing vector figure PDFs.

The layout JSON is private input. This code contains no manuscript content.
Pages are explicitly planned; overflow is an error, never silent shrink-to-fit.
"""
import argparse
import json
import re
import tempfile
from pathlib import Path
from xml.sax.saxutils import escape
import fitz
from reportlab.pdfgen import canvas
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from authoring_checks import sha256, json_write, resolve_input, different_output

def markdown_inline(text):
    text = escape(text)
    text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'\1 (\2)', text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<i>\1</i>', text)
    text = re.sub(r'`([^`]+)`', r'\1', text)
    return text

def markdown_blocks(text, default_style):
    if '```' in text or re.search(r'^\s*\|', text, re.M):
        raise ValueError('Fenced code and Markdown tables require a separate, explicitly planned layout')
    for block in re.split(r'\n\s*\n', text.strip()):
        if not block.strip():
            continue
        match = re.match(r'^(#{1,6})\s+(.*)$', block, re.S)
        if match:
            yield ('title' if len(match[1]) == 1 else 'heading'), markdown_inline(match[2].replace('\n', ' '))
        else:
            yield default_style, markdown_inline(block.replace('\n', ' '))

def assemble(layout_path, output, receipt):
    layout_path = Path(layout_path).resolve()
    output = Path(output).resolve()
    different_output(layout_path, output)
    different_output(layout_path, receipt)
    different_output(output, receipt)
    config = json.loads(layout_path.read_text(encoding='utf-8'))
    root = layout_path.parent
    page = config.get('page', {})
    pw, ph = page.get('width_mm', 210) * mm, page.get('height_mm', 297) * mm
    margin = page.get('margin_mm', 15) * mm
    top, bottom = page.get('top_mm', 11) * mm, page.get('bottom_mm', 17) * mm
    width = pw - 2 * margin
    minimum = float(config.get('minimum_font_pt', 8))
    for alias, filename in config.get('font_files', {}).items():
        pdfmetrics.registerFont(TTFont(alias, str(resolve_input(root, filename))))
    defaults = {'title': ('Helvetica-Bold', 14, 17), 'heading': ('Helvetica-Bold', 11, 14),
                'body': ('Helvetica', 10.5, 13.2), 'caption': ('Helvetica', 9, 11), 'small': ('Helvetica', 9, 11.5)}
    styles = {}
    for name, (font, size, leading) in defaults.items():
        params = {'fontName': font, 'fontSize': size, 'leading': leading, 'textColor': '#000000'}
        params.update(config.get('styles', {}).get(name, {}))
        if params['fontSize'] < minimum:
            raise ValueError('Style smaller than declared minimum: ' + name)
        styles[name] = ParagraphStyle(name, **params)
    pages = config['pages']
    if not pages:
        raise ValueError('At least one explicit page is required')
    output.parent.mkdir(parents=True, exist_ok=True)
    inputs, placements, blocks = {}, [], []
    with tempfile.TemporaryDirectory(prefix='pdf-layout-', dir=output.parent) as temporary:
        base = Path(temporary) / 'base.pdf'
        can = canvas.Canvas(str(base), pagesize=(pw, ph))
        metadata = config.get('metadata', {})
        can.setTitle(metadata.get('title', ''))
        can.setAuthor(metadata.get('author', ''))
        for number, spec in enumerate(pages, 1):
            if number > 1:
                can.showPage()
            y = ph - top
            for item in spec['blocks']:
                source = resolve_input(root, item['path'])
                different_output(source, output)
                different_output(source, receipt)
                inputs[str(source)] = {'basename': source.name, 'sha256': sha256(source)}
                if item['type'] == 'markdown':
                    for style, text in markdown_blocks(source.read_text(encoding='utf-8'), item.get('style', 'body')):
                        paragraph = Paragraph(text, styles[style])
                        _, height = paragraph.wrap(width, ph)
                        if y - height < bottom:
                            raise ValueError(f'Text overflow on page {number}: {source.name}')
                        paragraph.drawOn(can, margin, y - height)
                        blocks.append({'page': number, 'source_basename': source.name, 'style': style, 'bbox_pt': [margin, ph-y, margin+width, ph-y+height]})
                        y -= height + item.get('gap_pt', 8)
                elif item['type'] == 'figure_pdf':
                    with fitz.open(source) as doc:
                        figure_page = int(item.get('page_index', 0))
                        rect = doc[figure_page].rect
                        target_width = item.get('width_mm', width/mm) * mm
                        scale = target_width / rect.width
                        height = rect.height * scale
                        font_sizes = [s['size'] * scale for b in doc[figure_page].get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
                        if font_sizes and min(font_sizes) + .01 < minimum:
                            raise ValueError('Figure text becomes smaller than the declared minimum: ' + source.name)
                        if target_width > width + .01 or y - height < bottom:
                            raise ValueError(f'Figure overflow on page {number}: {source.name}')
                        left = margin + (width - target_width) / 2
                        placements.append({'page_index': number-1, 'source': source, 'source_page_index': figure_page, 'rect_pt': [left, ph-y, left+target_width, ph-y+height]})
                        y -= height + item.get('gap_pt', 5)
                else:
                    raise ValueError('Unsupported block type: ' + str(item['type']))
            footer = config.get('footer', {})
            if footer:
                size = footer.get('font_pt', 8)
                if size < minimum:
                    raise ValueError('Footer font smaller than declared minimum')
                can.setFont(footer.get('font', 'Helvetica'), size)
                footer_y = footer.get('y_mm', 9) * mm
                if footer.get('text_path'):
                    source = resolve_input(root, footer['text_path'])
                    different_output(source, output)
                    different_output(source, receipt)
                    value = source.read_text(encoding='utf-8').strip()
                    if '\n' in value or pdfmetrics.stringWidth(value, footer.get('font', 'Helvetica'), size) > width:
                        raise ValueError('Footer must fit on one line')
                    can.drawString(margin, footer_y, value)
                    inputs[str(source)] = {'basename': source.name, 'sha256': sha256(source)}
                if footer.get('page_numbers', True):
                    can.drawRightString(pw-margin, footer_y, str(number))
        can.save()
        with fitz.open(base) as doc:
            for placement in placements:
                with fitz.open(placement['source']) as source:
                    doc[placement['page_index']].show_pdf_page(fitz.Rect(placement['rect_pt']), source, placement['source_page_index'])
            doc.save(output, garbage=4, deflate=True)
    json_write(receipt, {'status': 'ASSEMBLED_PENDING_AUTOMATED_AND_MANUAL_QA', 'pdf_basename': output.name,
                        'pdf_sha256': sha256(output), 'page_count': len(pages), 'layout_sha256': sha256(layout_path),
                        'inputs': list(inputs.values()), 'text_blocks': blocks,
                        'figure_placements': [{**p, 'source': p['source'].name} for p in placements]})

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--layout', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--receipt', required=True, type=Path)
    args = parser.parse_args()
    assemble(args.layout, args.output, args.receipt)

if __name__ == '__main__':
    main()
