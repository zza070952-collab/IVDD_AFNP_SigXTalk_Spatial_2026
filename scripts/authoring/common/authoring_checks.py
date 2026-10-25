"""Content-independent hashes, numeric guards and DOCX package checks."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
XML_SPACE = '{http://www.w3.org/XML/1998/namespace}space'
TEXT_PART = re.compile(r'word/(?:document|footnotes|endnotes|header\d+|footer\d+)\.xml$')
NUMBER = re.compile(r'(?<![\w.])[+\-−]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?')

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def json_write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

def different_output(source, output):
    if Path(source).resolve() == Path(output).resolve():
        raise ValueError('Input and output must be different files')

def numeric_tokens(text, reference_patterns=()):
    for pattern in reference_patterns:
        text = re.sub(pattern, '', text)
    return NUMBER.findall(text)

def check_numbers(before, after, patterns=()):
    if numeric_tokens(before, patterns) != numeric_tokens(after, patterns):
        raise ValueError('Protected numeric token sequence changed outside permitted reference-label patterns')

def paragraph_nodes(paragraph):
    """Exclude nested text-box paragraphs to avoid double replacements."""
    result = []
    def visit(node):
        for child in node:
            if child.tag == W + 'p':
                continue
            if child.tag == W + 't':
                result.append(child)
            visit(child)
    visit(paragraph)
    return result

def paragraph_text(paragraph):
    return ''.join(n.text or '' for n in paragraph_nodes(paragraph))

def binary_package_hashes(path):
    with zipfile.ZipFile(path) as archive:
        return {name: hashlib.sha256(archive.read(name)).hexdigest()
                for name in archive.namelist()
                if name.startswith(('word/media/', 'word/embeddings/'))}

def docx_text_parts(path):
    with zipfile.ZipFile(path) as archive:
        return {name: [paragraph_text(p) for p in ET.fromstring(archive.read(name)).iter(W + 'p')]
                for name in archive.namelist() if TEXT_PART.fullmatch(name)}

def check_docx_preserved(before, after, require_text=True):
    if binary_package_hashes(before) != binary_package_hashes(after):
        raise ValueError('Embedded media or object bytes changed')
    if require_text and docx_text_parts(before) != docx_text_parts(after):
        raise ValueError('Visible DOCX text changed during formatting')

def resolve_input(base, raw):
    path = Path(raw).expanduser()
    return (path if path.is_absolute() else Path(base) / path).resolve()
