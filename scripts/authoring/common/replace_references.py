"""Apply an external literal replacement plan to Markdown or DOCX.

No replacement text is embedded here. Matching is simultaneous, exact and
scope-limited. Unchanged DOCX media/object members remain byte-identical.
"""
import argparse
import hashlib
import io
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
from authoring_checks import (W, XML_SPACE, TEXT_PART, sha256, json_write, different_output,
                              check_numbers, paragraph_nodes, paragraph_text, check_docx_preserved)

def permitted(rule, text, index, part, kind):
    scope=rule.get('scope',{})
    if part is not None and part not in scope.get('parts',['word/document.xml']):
        return False
    indices=scope.get('paragraph_indices' if kind=='docx' else 'line_indices')
    if indices is not None and index not in indices:
        return False
    pattern=scope.get('paragraph_regex' if kind=='docx' else 'line_regex')
    return pattern is None or re.search(pattern,text) is not None

def matches(text, rules, index, part, kind):
    available=[r for r in rules if permitted(r,text,index,part,kind)]
    by_old={r['old']:r for r in available}
    if len(by_old)!=len(available):
        raise ValueError('Ambiguous duplicate old strings in overlapping scopes')
    if not by_old:
        return []
    pattern=re.compile('|'.join(re.escape(x) for x in sorted(by_old,key=len,reverse=True)))
    return [(m.start(),m.end(),by_old[m[0]]) for m in pattern.finditer(text)]

def edit_text(text, found):
    for start,end,rule in reversed(found):
        text=text[:start]+rule['new']+text[end:]
    return text

def patch_nodes(nodes, start, end, replacement):
    offset=0
    touched=[]
    for node in nodes:
        value=node.text or ''
        stop=offset+len(value)
        if start<stop and end>offset:
            touched.append((node,max(0,start-offset),min(len(value),end-offset)))
        offset=stop
    if not touched:
        raise ValueError('Replacement span did not map to DOCX text nodes')
    for i,(node,a,b) in enumerate(touched):
        value=node.text or ''
        node.text=value[:a]+(replacement if i==0 else '')+value[b:]
        if node.text.startswith(' ') or node.text.endswith(' '):
            node.set(XML_SPACE,'preserve')

def apply(input_path, output_path, plan_path, receipt, dry_run=False):
    input_path, output_path, plan_path, receipt=map(Path,(input_path,output_path,plan_path,receipt))
    for source in [input_path,plan_path]:
        different_output(source,output_path); different_output(source,receipt)
    different_output(output_path,receipt)
    plan=json.loads(plan_path.read_text(encoding='utf-8'))
    if plan.get('input_sha256') and plan['input_sha256']!=sha256(input_path):
        raise ValueError('Replacement plan is bound to a different input hash')
    rules=plan['rules']
    if not rules or len({r['id'] for r in rules})!=len(rules):
        raise ValueError('Nonempty uniquely identified rules are required')
    if any(not r['old'] or not isinstance(r.get('expected'),int) or r['expected']<0 for r in rules):
        raise ValueError('Each rule needs a nonempty old string and exact expected match count')
    patterns=plan.get('reference_label_patterns',[])
    for pattern in patterns:
        re.compile(pattern)
    counts={r['id']:0 for r in rules}
    changes=[]
    replacements={}
    def record(before,after,found,index,part):
        check_numbers(before,after,patterns)
        for _,_,rule in found:counts[rule['id']]+=1
        changes.append({'part':part,'index':index,'rule_ids':[r['id'] for _,_,r in found],
                        'before_text_sha256':hashlib.sha256(before.encode()).hexdigest(),
                        'after_text_sha256':hashlib.sha256(after.encode()).hexdigest()})
    suffix=input_path.suffix.lower()
    if suffix in ['.md','.txt']:
        lines=input_path.read_bytes().decode('utf-8').splitlines(keepends=True)
        edited=[]
        for index,line in enumerate(lines):
            found=matches(line,rules,index,None,'md')
            after=edit_text(line,found)
            if found:record(line,after,found,index,None)
            edited.append(after)
        output_bytes=''.join(edited).encode('utf-8')
    elif suffix=='.docx':
        with zipfile.ZipFile(input_path) as archive:
            for name in archive.namelist():
                if not TEXT_PART.fullmatch(name):continue
                raw=archive.read(name)
                for _,ns in ET.iterparse(io.BytesIO(raw),events=['start-ns']):
                    prefix,uri=ns
                    try:ET.register_namespace(prefix,uri)
                    except ValueError:pass
                root=ET.fromstring(raw)
                changed=False
                for index,paragraph in enumerate(root.iter(W+'p')):
                    before=paragraph_text(paragraph)
                    found=matches(before,rules,index,name,'docx')
                    if not found:continue
                    if any(next(paragraph.iter(W+tag),None) is not None for tag in ['fldChar','instrText','fldSimple']):
                        raise ValueError('Matched paragraph contains live Word fields; flatten or edit fields explicitly in a separate authorized step')
                    after=edit_text(before,found)
                    record(before,after,found,index,name)
                    nodes=paragraph_nodes(paragraph)
                    for start,end,rule in reversed(found):patch_nodes(nodes,start,end,rule['new'])
                    if paragraph_text(paragraph)!=after:raise ValueError('DOCX text-node edit mismatch')
                    changed=True
                if changed:replacements[name]=ET.tostring(root,encoding='utf-8',xml_declaration=True)
    else:
        raise ValueError('Only .md, .txt and .docx are supported')
    for rule in rules:
        if counts[rule['id']]!=rule['expected']:
            raise ValueError(f"Rule {rule['id']}: expected {rule['expected']} matches, found {counts[rule['id']]}")
    if not dry_run:
        output_path.parent.mkdir(parents=True,exist_ok=True)
        if suffix=='.docx':
            with zipfile.ZipFile(input_path) as archive,zipfile.ZipFile(output_path,'w') as target:
                for item in archive.infolist():target.writestr(item,replacements.get(item.filename,archive.read(item.filename)))
            check_docx_preserved(input_path,output_path,require_text=False)
        else:output_path.write_bytes(output_bytes)
    json_write(receipt,{'status':'DRY_RUN_SCOPE_CHECKS_PASS' if dry_run else 'APPLIED_PENDING_RENDER_AND_VISUAL_REVIEW',
                        'input_basename':input_path.name,'input_sha256':sha256(input_path),'plan_sha256':sha256(plan_path),
                        'output_basename':output_path.name,'output_sha256':None if dry_run else sha256(output_path),
                        'rule_match_counts':counts,'changes':changes,'numeric_sequence_guard':'PASS',
                        'embedded_media_bytes':'UNCHANGED' if suffix=='.docx' and not dry_run else 'NOT_APPLICABLE_OR_DRY_RUN',
                        'visual_review':'NOT_PERFORMED_BY_THIS_COMMAND'})

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',required=True,type=Path);parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--plan',required=True,type=Path);parser.add_argument('--receipt',required=True,type=Path)
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    apply(args.input,args.output,args.plan,args.receipt,args.dry_run)

if __name__=='__main__':main()
