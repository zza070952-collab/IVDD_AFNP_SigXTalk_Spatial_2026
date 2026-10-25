#!/usr/bin/env python3
"""Read-only DOCX inventory; explicit private output contains source text."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
from lxml import etree
from figure_references import extract
from output_guard import require_private_output

NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'a': 'http://schemas.openxmlformats.org/drawingml/2006/main',
      'wp': 'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}


def inventory(path):
    path = Path(path)
    before = path.read_bytes()
    records = []
    with zipfile.ZipFile(path) as package:
        root = etree.fromstring(package.read('word/document.xml'))
        for i, paragraph in enumerate(root.xpath('.//w:body/w:p', namespaces=NS)):
            text = ''.join(paragraph.xpath('.//w:t/text()', namespaces=NS))
            style = paragraph.xpath('./w:pPr/w:pStyle/@w:val', namespaces=NS)
            images = paragraph.xpath('.//a:blip/@r:embed', namespaces={**NS, 'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'})
            records.append({'index': i, 'style': style[0] if style else '', 'text': text,
                            'images': images, 'extents_emu': [dict(x.attrib) for x in paragraph.xpath('.//wp:extent', namespaces=NS)],
                            'main_figure_refs': [dict(text=r['text'], start=r['start'], end=r['end'], items=[x.label() for x in r['items']]) for r in extract(text)]})
        fields = root.xpath('.//w:instrText/text()', namespaces=NS)
        media = {name: {'bytes': len(package.read(name)), 'sha256': hashlib.sha256(package.read(name)).hexdigest()}
                 for name in package.namelist() if name.startswith('word/media/')}
        relationships = [dict(r.attrib) for r in etree.fromstring(package.read('word/_rels/document.xml.rels'))]
    if path.read_bytes() != before:
        raise ValueError('Read-only source changed')
    return {'sha256': hashlib.sha256(before).hexdigest(), 'paragraphs': records, 'fields': fields,
            'media': media, 'relationships': relationships, 'source_input_unchanged': True}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('docx', type=Path)
    parser.add_argument('--output', type=Path, required=True, help='Private JSON destination outside publication checkout')
    args = parser.parse_args()
    target = require_private_output(args.output)
    if target == args.docx.resolve() or target.exists():
        raise ValueError('Choose a fresh private output, separate from the source')
    result = inventory(args.docx)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps({'sha256': result['sha256'], 'paragraphs': len(result['paragraphs']), 'images': len(result['media'])}))
