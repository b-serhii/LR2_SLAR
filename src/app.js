import {sizes, generateSystem,diagonalDominance,makeZeroPivotCase,parseImportedText} from './system.js';
import {gaussPivot,gaussWithoutPivot,jacobi,seidel,analyze} from './solvers.js';
import {drawConvergence,drawSolutions} from './charts.js';
const $=s=>document.querySelector(s);
const f=(x,p=3)=>x===null||x===undefined?'—':Number.isFinite(x)?Number(x.toFixed(p)).toLocaleString('uk-UA',{maximumFractionDigits:p,minimumFractionDigits:p}):'—';
const sci=x=>x===null||x===undefined?'—':Number.isFinite(x)?x.toExponential(2):'—';
const k=(n,d)=>`${n}|${d}`;
let selected=5,d=1,eps=1e-10,maxIterations=20000;
let custom=new Map(),results=new Map(),current=null,frameIndex=0,playToken=0,playing=false;

function getSystem(n,di=d){return custom.get(k(n,di))||generateSystem(n,di);}
function message(msg,type='ok'){
 const node=$('#notice');node.className='alert'+(type==='warn'?' warn':type==='error'?' error':'');
 node.innerHTML=`<strong>${type==='ok'?'✓':type==='warn'?'!':'×'}</strong><span></span>`;
 node.querySelector('span').textContent=msg;
}
function options(){return {iteration:{eps,maxIterations},gauss:{trace:true}};}
function model(n,di=d){return results.get(k(n,di));}
function computeAll(){
 stop();
 const m=Number($('#max-iters').value);
 if(!Number.isInteger(m)||m<20000||m>100000)throw new Error('Максимальна кількість ітерацій має бути від 20 000 до 100 000.');
 maxIterations=m;eps=Number($('#eps').value);
 const newResults=new Map();
 for(const di of [1,.1])for(const n of sizes)newResults.set(k(n,di),analyze(getSystem(n,di),options()));
 results=newResults;
 current=model(selected);
 updateAll();
}
function updateSizeButtons(){
 document.querySelectorAll('button[data-size]').forEach(b=>{const yes=Number(b.dataset.size)===selected;b.classList.toggle('active',yes);b.setAttribute('aria-pressed',String(yes));});
}
function renderMatrix(snapshot=null){
 const M=snapshot?.A??current.system.A,B=snapshot?.b??current.system.b,n=M.length;
 let html='<table class="matrix-grid"><thead><tr><th>i / j</th>';
 for(let j=0;j<n;j++)html+=`<th>${j+1}</th>`;
 html+='<th>bᵢ</th></tr></thead><tbody>';
 for(let i=0;i<n;i++){
  html+=`<tr><th>${i+1}</th>`;
  for(let j=0;j<n;j++){
   const value=M[i][j];
   let cls='';
   if(snapshot?.type==='pivot' && j===snapshot.k && i===snapshot.selectedRow)cls='choose';
   else if(snapshot && (snapshot.type==='swap'||snapshot.type==='eliminate') && i===snapshot.k && j===snapshot.k)cls='highlight';
   else if(snapshot?.type==='eliminate' && j===snapshot.k && i>snapshot.k)cls='focus';
   else if(Math.abs(value)<1e-12)cls='zero';
   else if(value<0)cls='negative';
   html+=`<td class="${cls}" title="a${i+1},${j+1}=${value}">${f(value,3)}</td>`;
  }
  html+=`<td class="rhs" title="b${i+1}=${B[i]}">${f(B[i],3)}</td></tr>`;
 }
 html+='</tbody></table>';
 $('#matrix').innerHTML=html;
}
function updateStats(){
 const {system,gauss,jacobi:j,seidel:s,dominance}=current;
 $('#st-n').textContent=`${selected} × ${selected}`;
 $('#st-d').textContent=custom.has(k(selected,d))?'Користувацька система':`d = ${d===1?'1':'0,1'}`;
 $('#st-g').textContent=sci(gauss.residual);
 $('#st-j').textContent=j.iterations;
 $('#st-jstatus').textContent=j.converged?'Досягнуто ε':'Ліміт / розбіжність';
 $('#st-s').textContent=s.iterations;
 $('#st-sstatus').textContent=s.converged?'Досягнуто ε':'Ліміт / розбіжність';
 const diag=diagonalDominance(system.A);
 if(diag.strict)message(`Матриця має строгу діагональну перевагу. Мінімальний запас: ${f(diag.minMargin,4)}. Якобі та Зейдель збігаються.`);
 else message('Строгої діагональної переваги немає: збіжність ітераційних методів не гарантована. Перевірте статус і результати.','warn');
}
function updatePlots(){
 for(const n of sizes){
  const r=model(n);
  const conv=$(`[data-card="conv-${n}"]`),sol=$(`[data-card="sol-${n}"]`);
  drawConvergence(conv.querySelector('.plot-wrap'),r.jacobi,r.seidel,{small:true});
  drawSolutions(sol.querySelector('.plot-wrap'),r,{small:true});
  conv.style.borderColor=sol.style.borderColor=n===selected?'#a99ff0':'#e4e9f2';
 }
}
function updateResults(){
 const arr=[];
 for(const n of sizes)for(const di of [1,.1]){
  const m=model(n,di);
  for(const r of [m.gauss,m.jacobi,m.seidel]){
   arr.push(`<tr class="${n===selected&&di===d?'focus':''}"><td>${n}</td><td>${di===1?'1':'0,1'}</td><td>${r.method}</td><td>${r.method==='Гаусс'?'—':r.iterations}</td><td>${sci(r.residual)}</td><td>${sci(r.error)}</td><td>${f(r.timeMs,4)}</td><td class="${r.converged?'ok':'bad'}">${r.converged?'✓':'Не збігся'}</td></tr>`);
  }
 }
 $('#results').innerHTML=arr.join('');
}
function methodTrace(){
 const method=$('#anim-method').value;
 if(method==='gauss')return current.gauss.steps;
 return (method==='jacobi'?current.jacobi:current.seidel).history;
}
function refreshFrame(){
 const method=$('#anim-method').value;
 const seq=methodTrace(),s=seq[Math.min(seq.length-1,frameIndex)];
 if(method==='gauss'){
  renderMatrix(s);
  $('#matrix-status').textContent=`Гаусс · ${s?.type??'початок'}`;
  $('#step-label').textContent=`Крок ${frameIndex+1}/${seq.length}: ${s?.note??''}`;
  updatePlots();
 }else {
  renderMatrix();
  $('#matrix-status').textContent=`${method==='jacobi'?'Якобі':'Зейдель'} · ітерація ${s?.iteration??0}`;
  const sol=s?.x||[];
  $('#step-label').textContent=`Ітерація ${s.iteration}; ‖Ax−b‖₂ = ${sci(s.residual)}; x = [${sol.slice(0,Math.min(4,sol.length)).map(x=>f(x,3)).join('; ')}${sol.length>4?'; …':''}]`;
  const r=method==='jacobi'?current.jacobi:current.seidel;
  const partial=r.history.slice(0,frameIndex+1);
  const partialJ=method==='jacobi'?partial:null,partialS=method==='seidel'?partial:null;
  // Автомасштабування під час анімації за поточними значеннями активних кривих.
  const c=$(`[data-card="conv-${selected}"] .plot-wrap`),v=$(`[data-card="sol-${selected}"] .plot-wrap`);
  drawConvergence(c,current.jacobi,current.seidel,{small:true,partialJ:partialJ??current.jacobi.history,partialS:partialS??current.seidel.history});
  drawSolutions(v,current,{small:true,shown:{[method]:sol}});
 }
}
function stop(){playToken++;playing=false;if($('#play'))$('#play').textContent='▶ Відтворити';}
function resetFrame(){stop();frameIndex=0;refreshFrame();}
function nextFrame(){
 const seq=methodTrace();
 if(frameIndex<seq.length-1)frameIndex++;
 refreshFrame();
}
function animate(){
 if(playing){stop();return;}
 const seq=methodTrace();if(frameIndex>=seq.length-1)frameIndex=0;
 playing=true;const mine=++playToken;$('#play').textContent='Ⅱ Пауза';
 let last=0;
 function tick(ts){
  if(mine!==playToken)return;
  if(!last||ts-last>=Number($('#speed').value)){
   nextFrame();last=ts;
   if(frameIndex>=seq.length-1){stop();return;}
  }
  requestAnimationFrame(tick);
 }
 requestAnimationFrame(tick);
}
function updateAll(){
 updateSizeButtons();updateStats();updatePlots();updateResults();frameIndex=0;refreshFrame();
}
$('#sizes').addEventListener('click',e=>{
 const b=e.target.closest('button[data-size]');if(!b)return;
 stop();selected=Number(b.dataset.size);current=model(selected);updateAll();
});
$('#dominance').addEventListener('change',e=>{d=Number(e.target.value);current=model(selected);updateAll();});
$('#run').addEventListener('click',()=>{try{computeAll();}catch(e){message(e.message,'error');}});
$('#generate').addEventListener('click',()=>{custom.clear();$('#json-input').value='';$('#json-error').textContent='';try{computeAll();}catch(e){message(e.message,'error');}});
$('#anim-method').addEventListener('change',resetFrame);
$('#play').addEventListener('click',animate);
$('#next').addEventListener('click',()=>{stop();nextFrame();});
$('#reset').addEventListener('click',resetFrame);
$('#fill-json').addEventListener('click',()=>{$('#json-input').value=JSON.stringify(getSystem(selected),null,2);});
$('#load-json').addEventListener('click',()=>{
  try{
   const data=parseImportedText($('#json-input').value);
   selected=data.n;
   custom.set(k(selected,d),data);
   $('#json-error').textContent='Власні дані застосовано.';
   computeAll();
  }catch(e){$('#json-error').textContent=`Помилка: ${e.message}`;}
});
$('#file-input').addEventListener('change',async e=>{
 const file=e.target.files?.[0];if(!file)return;
 try{$('#json-input').value=await file.text();$('#json-error').textContent=`Завантажено ${file.name}. Натисніть «Застосувати дані».`;}
 catch(err){$('#json-error').textContent=err.message;}
});
$('#zero-pivot').addEventListener('click',()=>{
 try{
  const sample=generateSystem(5,1),z=makeZeroPivotCase(sample);
  let err='';try{gaussWithoutPivot(z);}catch(e){err=e.message;}
  const fixed=gaussPivot(z);
  $('#zero-result').innerHTML=`<strong>До виправлення:</strong> ${err||'збою не виявлено'}<div class="divider"></div><strong>Після додавання вибору головного елемента:</strong> перестановок рядків = ${fixed.swaps}; ‖Ax−b‖₂ = ${sci(fixed.residual)}; ‖x−x*‖₂ = ${sci(fixed.error)}.`;
 }catch(e){$('#zero-result').textContent=`Помилка: ${e.message}`;}
});
$('#export').addEventListener('click',()=>{
 const lines=['n;d;Метод;Ітерації;Невязка_L2;Похибка_L2;Час_мс;Збіжність'];
 for(const n of sizes)for(const di of [1,.1]){
  const m=model(n,di);
  for(const r of [m.gauss,m.jacobi,m.seidel])lines.push([n,di,r.method,r.iterations,r.residual,r.error??'',r.timeMs,r.converged?'так':'ні'].join(';'));
 }
 const blob=new Blob(['\uFEFF'+lines.join('\n')],{type:'text/csv;charset=utf-8'});
 const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download='LR2_Slar_metrics_variant1.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1500);
});
try{computeAll();}catch(e){message(e.message,'error');}
