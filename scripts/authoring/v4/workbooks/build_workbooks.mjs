import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const payloadArg=process.argv.find(v=>v.startsWith('--payload='))?.slice('--payload='.length);
if(!payloadArg)throw new Error('Supply an external private --payload=<file>');
const payload=JSON.parse(await fs.readFile(payloadArg,'utf8'));
const root=payload.root;
const qa=path.join(root,'08_review','workbook_qa');
await fs.mkdir(qa,{recursive:true});
let audit=[];
try{audit=JSON.parse(await fs.readFile(path.join(qa,'WORKBOOK_ARTIFACT_BUILD.json'),'utf8'));}catch(e){if(e.code!=='ENOENT')throw e;}
const selected=process.argv.find(v=>v.startsWith('--books='))?.split('=')[1].split(',').map(Number);
const colName=(n)=>{let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
function fmt(col,type){
  if(type==='datetime')return 'yyyy-mm-dd hh:mm:ss';
  if(type==='integer')return '#,##0';
  if(type!=='number')return '@';
  if(/(^p$|^p_|_p$|p_adjusted|p_two|pvalue|p_value|fdr|q_value)/i.test(col))return '0.000E+00';
  if(/Ct|OD450|fold|ratio|fraction|score|aucell|ucell|PRS|specificity|fidelity|mean|median|CI95/i.test(col))return '0.0000';
  return '0.000000';
}
function width(header,rows,j,type){
  if(type==='number'||type==='integer')return Math.max(16,Math.min(24,header.length+3));
  const longest=Math.max(header.length,...rows.slice(0,200).map(r=>String(r[j]??'').length));
  if(/sha256/i.test(header))return 70;
  if(/source|provenance|notes|interpretation|limitation|rule|transform|missing|status/i.test(header))return Math.min(82,Math.max(28,longest+2));
  return Math.min(50,Math.max(18,longest+2));
}
function fill(sheet,headers,rows,types,tableName){
  const last=colName(headers.length-1),bottom=rows.length+1;
  // Set literal text before the block write. This prevents numeric-looking IDs,
  // mixed QC labels and ISO text from being auto-coerced by the writer.
  sheet.getRange(`A1:${last}${bottom}`).setNumberFormat('@');
  const values=rows.map(r=>r.map((v,j)=>v!==null&&types[j]==='datetime'?new Date(v):v));
  sheet.getRange(`A1:${last}${bottom}`).values=[headers,...values];
  const used=sheet.getRange(`A1:${last}${bottom}`);
  used.format.font={name:'Arial',size:10,color:'#20252B'};
  used.format.verticalAlignment='center';
  used.format.rowHeight=20;
  sheet.showGridLines=false;
  sheet.getRange(`A1:${last}1`).format={fill:'#34495E',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:48};
  for(let j=0;j<headers.length;j++){
    const c=colName(j),t=types[j];
    sheet.getRange(`${c}1:${c}${bottom}`).format.columnWidth=width(headers[j],rows,j,t);
    if(rows.length){
      const r=sheet.getRange(`${c}2:${c}${bottom}`);
      r.setNumberFormat(fmt(headers[j],t));
      r.format.horizontalAlignment=(t==='number'||t==='integer')?'right':'left';
      r.format.wrapText=t==='text';
    }
  }
  // Fit only the populated body. Short source rows stay compact; descriptions remain readable.
  if(rows.length)sheet.getRange(`A2:${last}${bottom}`).format.autofitRows();
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(1);
  const table=sheet.tables.add(`A1:${last}${bottom}`,true,tableName);
  table.showFilterButton=true;
  return {range:`A1:${last}${bottom}`,bottom,last};
}
for(let bi=0;bi<payload.books.length;bi++){
  if(selected&&!selected.includes(bi+1))continue;
  const b=payload.books[bi];
  const wb=Workbook.create();
  const idx=wb.worksheets.add('Index');
  idx.tabColor=b.internal?'#9C6500':'#34495E';
  const indexRows=b.sheets.map(s=>[s.name,s.panel,s.rows.length,s.headers.length,s.source,s.source_sha256,s.analysis_unit,s.independent_n,s.transformation,s.status,s.missing]);
  const indexHead=['Sheet','Figure_or_table','Data_rows','Fields','Source_file','Source_SHA256','Observation_unit','Independent_unit','Transformation','Release_status','Missing_values'];
  if(b.internal){
    const status=indexHead.splice(9,1)[0];indexHead.unshift(status);
    for(const r of indexRows){const status=r.splice(9,1)[0];r.unshift(status);}
  }
  const indexTypes=indexHead.map(h=>['Data_rows','Fields'].includes(h)?'integer':'text');
  fill(idx,indexHead,indexRows,indexTypes,`B${bi}Index`);
  const models=[{name:'Index',rows:indexRows.length,cols:indexHead.length}];
  const dict=[];
  for(let si=0;si<b.sheets.length;si++){
    const s=b.sheets[si];
    const sheet=wb.worksheets.add(s.name);
    fill(sheet,s.headers,s.rows,s.types,`B${bi}T${si}`);
    for(let j=0;j<s.headers.length;j++)dict.push([s.name,s.headers[j],s.types[j],s.units[j],s.missing,s.source]);
    models.push({name:s.name,rows:s.rows.length,cols:s.headers.length,source:s.source,source_sha256:s.source_sha256,typed_matrix_sha256:s.typed_matrix_sha256});
  }
  const ds=wb.worksheets.add('Data_dictionary');
  fill(ds,['Sheet','Field','Stored_type','Units_or_definition','Missing_values','Source_file'],dict,['text','text','text','text','text','text'],`B${bi}Dictionary`);
  models.push({name:'Data_dictionary',rows:dict.length,cols:6});
  wb.recalculate();
  const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:1500});
  const inspections=[];
  for(const m of models){
    const range=`A1:${colName(Math.min(m.cols,6)-1)}${Math.min(m.rows+1,7)}`;
    const preview=await wb.render({sheetName:m.name,range,scale:1.3,format:'png'});
    const png=`book${bi+1}_${m.name}.png`;
    await fs.writeFile(path.join(qa,png),new Uint8Array(await preview.arrayBuffer()));
    const check=await wb.inspect({kind:'table',range:`'${m.name}'!${range}`,include:'values,formulas',tableMaxRows:3,tableMaxCols:5,maxChars:600});
    inspections.push({sheet:m.name,range,preview:png,inspection:check.ndjson});
  }
  const out=path.join(root,b.file);
  await fs.mkdir(path.dirname(out),{recursive:true});
  await (await SpreadsheetFile.exportXlsx(wb)).save(out);
  try{await fs.rename(out+'.inspect.ndjson',path.join(qa,`book${bi+1}_export_inspect.ndjson`));}catch(e){if(e.code!=='ENOENT')throw e;}
  const bytes=await fs.readFile(out);
  const entry={file:b.file,bytes:bytes.length,sha256:crypto.createHash('sha256').update(bytes).digest('hex'),sheets:models.length,source_sheets:b.sheets.length,data_rows:b.sheets.reduce((a,s)=>a+s.rows.length,0),formula_error_scan:errors.ndjson,rendered_sheets:inspections.map(x=>({sheet:x.sheet,range:x.range,preview:x.preview})),internal_source_unresolved:b.internal,models};
  audit=audit.filter(e=>e.file!==b.file);audit.push(entry);audit.sort((a,b)=>a.file.localeCompare(b.file));
  await fs.writeFile(path.join(qa,`book${bi+1}_inspections.json`),JSON.stringify(inspections,null,2));
  await fs.writeFile(path.join(qa,'WORKBOOK_ARTIFACT_BUILD.json'),JSON.stringify(audit,null,2));
  console.log(JSON.stringify({file:b.file,bytes:bytes.length,sheets:models.length,rendered:inspections.length}));
}
