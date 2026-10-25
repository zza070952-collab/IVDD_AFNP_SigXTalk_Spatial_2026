"""Extract same-paragraph oligo pairs without redistributing a private protocol."""
from __future__ import annotations
import re
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

def docx_paragraphs(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as archive:
        root=ET.fromstring(archive.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    result=[]
    for paragraph in root.findall('.//w:body/w:p',ns):
        values=[]
        for node in paragraph.iter():
            if node.tag.endswith('}t'):values.append(node.text or '')
            elif node.tag.endswith('}br'):values.append('\n')
            elif node.tag.endswith('}tab'):values.append('\t')
        result.append(''.join(values))
    return result

def parse_sequence_pairs(paragraphs: list[str]) -> list[dict]:
    """Pair S/AS or F/R within one paragraph only; never infer a missing strand."""
    result=[];keys=set()
    for i,text in enumerate(paragraphs,1):
        m=re.fullmatch(r'(si[A-Za-z0-9_]+)\s+S\s+([ACGU]+dTdT)\s+AS\s+([ACGU]+dTdT)\s*',text)
        if m:
            row={'record':m[1],'forward_5to3':m[2],'reverse_5to3':m[3],'record_type':'siRNA sense/antisense','primerbank_id':None,'amplicon_bp':None,'source_paragraph':i}
        else:
            m=re.fullmatch(r'([^（(\n]+)[（(]([0-9]+c[0-9]+)[，,]\s*([0-9]+)\s*bp[）)]\s+F\s+([ACGT]+)\s+R\s+([ACGT]+)\s*',text)
            if not m:continue
            row={'record':m[1].strip(),'forward_5to3':m[4],'reverse_5to3':m[5],'record_type':'qPCR primer','primerbank_id':m[2],'amplicon_bp':int(m[3]),'source_paragraph':i}
        key=(row['record_type'],row['record'])
        if key in keys:raise ValueError('Duplicate sequence record: '+str(key))
        keys.add(key);result.append(row)
    return result
