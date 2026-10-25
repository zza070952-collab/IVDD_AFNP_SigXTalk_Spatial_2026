import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';
const payloadArg=process.argv.find(x=>x.startsWith('--payload='))?.slice(10);
if(!payloadArg)throw new Error('Required --payload=/private/build/COMPACT_WORKBOOK_PAYLOAD.json');
const p=JSON.parse(await fs.readFile(payloadArg,'utf8'));
const outdir=p.output,qa=p.qa;
await fs.mkdir(outdir,{recursive:true});await fs.mkdir(qa,{recursive:true});
const col=(n)=>{let s='';for(n++;n>0;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;};
const fmt=(h,t)=>t==='integer'?'#,##0':t!=='number'?'@':/p_two|p_adj|p_raw|p_value|algorithm_p|(^P$)/i.test(h)?'0.000E+00':/Raw_Ct/i.test(h)?'0.000':/OD450/i.test(h)?'0.0000':'0.000000';
const audits=[];
const selected=process.argv.find(v=>v.startsWith('--book='))?.split('=')[1];
for(let bi=0;bi<p.books.length;bi++){
 const b=p.books[bi];if(selected&&String(bi+1)!==selected)continue;
 const wb=Workbook.create(),models=[],previews=[];
 for(let si=0;si<b.sheets.length;si++){
  const s=b.sheets[si],ws=wb.worksheets.add(s.name),last=col(s.headers.length-1),end=s.rows.length+5;
  ws.showGridLines=false;
  const data=ws.getRange(`A5:${last}${end}`);data.setNumberFormat('@');data.values=[s.headers,...s.rows];
  data.format.font={name:'Arial',size:10,color:'#20252B'};data.format.verticalAlignment='center';data.format.rowHeight=22;
  ws.getRange('A1').values=[[s.name.replaceAll('_',' ')]];ws.getRange('A1').format.font={name:'Arial',size:14,bold:true,color:'#20252B'};
  ws.getRange('A2').values=[['Internal review. Experimental source verification remains pending.']];ws.getRange('A2').format.font={name:'Arial',size:10,bold:true,color:'#8B4B00'};
  ws.getRange('A3').values=[['Static source values. Definitions and source hashes follow the data table.']];ws.getRange('A3').format.font={name:'Arial',size:10,color:'#444444'};
  ws.getRange(`A5:${last}5`).format={fill:'#34495E',font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:46};
  for(let j=0;j<s.headers.length;j++){
   const h=s.headers[j],t=s.types[j],c=col(j);
   let width=t==='number'||t==='integer'?Math.max(17,Math.min(25,h.length+2)):Math.min(42,Math.max(18,h.length+2,...s.rows.slice(0,100).map(r=>Math.min(42,String(r[j]??'').length+2))));
   if(/definition|notes|source_status|interpretation|limitation|rule|description/i.test(h))width=52;
   if(/sha256/i.test(h))width=70;
   ws.getRange(`${c}5:${c}${end}`).format.columnWidth=width;
   const body=ws.getRange(`${c}6:${c}${end}`);body.setNumberFormat(fmt(h,t));body.format.horizontalAlignment=t==='text'?'left':'right';body.format.wrapText=t==='text';
  }
  ws.getRange(`A6:${last}${end}`).format.autofitRows();
  const table=ws.tables.add(`A5:${last}${end}`,true,`CompactB${bi}T${si}`);table.showFilterButton=true;
  ws.freezePanes.freezeRows(5);ws.freezePanes.freezeColumns(1);
  const notes=[['Notes',s.description||''],['Source',s.source||''],['Source SHA256',s.source_sha256||'See record-specific source SHA256 fields'],['Missing values','Literal NA/NA reason means missing or inapplicable. Recorded zeros remain zeros. Raw undetermined Ct remains blank with original Cq_status.'],['Units and definitions','Column names and record-specific Units/Definition/Notes fields specify measurement units and inferential units. No formulas or hidden data.']];
  let nr=end+3;
  for(const [label,value] of notes){
   const lines=String(value).match(/.{1,105}(?:\s|$)|\S{1,105}/g)||[''];
   ws.getRange(`A${nr}`).values=[[label]];ws.getRange(`A${nr}`).format.font={name:'Arial',size:10,bold:true};nr++;
   for(const line of lines){ws.getRange(`A${nr}`).values=[[line.trim()]];ws.getRange(`A${nr}`).format.font={name:'Arial',size:10};nr++;}
   nr++;
  }
  const hash=crypto.createHash('sha256').update(JSON.stringify([s.headers,...s.rows])).digest('hex');
  models.push({name:s.name,rows:s.rows.length,columns:s.headers.length,range:`A5:${last}${end}`,typed_matrix_sha256:hash});
  for(let left=0;left<s.headers.length;left+=5){
   const range=`${col(left)}5:${col(Math.min(left+4,s.headers.length-1))}${Math.min(end,12)}`;
   const preview=await wb.render({sheetName:s.name,range,scale:1.5,format:'png'});
   const file=`book${bi+1}_${s.name}_cols${left+1}.png`;
   await fs.writeFile(path.join(qa,file),new Uint8Array(await preview.arrayBuffer()));previews.push({sheet:s.name,range,file});
  }
  const top=await wb.render({sheetName:s.name,range:'A1:E9',scale:1.2,format:'png'});await fs.writeFile(path.join(qa,`book${bi+1}_${s.name}_opening.png`),new Uint8Array(await top.arrayBuffer()));
 }
 wb.recalculate();
 const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:20},maxChars:1500});
 const out=path.join(outdir,b.filename);await (await SpreadsheetFile.exportXlsx(wb)).save(out);
 try{await fs.rename(out+'.inspect.ndjson',path.join(qa,`book${bi+1}_export_inspect.ndjson`));}catch(e){if(e.code!=='ENOENT')throw e;}
 const bytes=await fs.readFile(out);audits.push({filename:b.filename,bytes:bytes.length,sha256:crypto.createHash('sha256').update(bytes).digest('hex'),sheet_count:b.sheets.length,data_rows:b.sheets.reduce((n,s)=>n+s.rows.length,0),models,previews,error_scan:errors.ndjson});
 console.log(JSON.stringify({file:b.filename,sheets:b.sheets.length,rows:b.sheets.reduce((n,s)=>n+s.rows.length,0),bytes:bytes.length}));
}
await fs.writeFile(path.join(qa,selected?`ARTIFACT_WORKBOOK_BUILD_${selected}.json`:'ARTIFACT_WORKBOOK_BUILD.json'),JSON.stringify(audits,null,2));
