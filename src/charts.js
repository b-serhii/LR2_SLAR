/** SVG без сторонніх бібліотек. Усі межі визначаємо за поточними даними. */
const C={jacobi:'#6554D9',seidel:'#0D8D88',gauss:'#E0A24B',exact:'#24395B'};
const v=(x,p=2)=>Number.isFinite(x)?Number(x.toFixed(p)).toLocaleString('uk-UA'):'—';
function bound(nums,fraction=.12){
 const finite=nums.filter(Number.isFinite);
 if(!finite.length)return [-1,1];
 const min=Math.min(...finite),max=Math.max(...finite),w=Math.max(max-min,Math.max(1e-6,Math.abs(max)*.03));
 return [min-w*fraction,max+w*fraction];
}
function frame(SMALL){ const w=SMALL?570:1000,h=SMALL?255:385;
 return {w,h,L:SMALL?56:75,T:SMALL?21:32,R:SMALL?15:28,B:SMALL?37:48}; }
function render(svg,caption,view,body){
 svg.innerHTML=`<svg class="plot" viewBox="0 0 ${view.w} ${view.h}" role="img" aria-label="${caption}" xmlns="http://www.w3.org/2000/svg">${body}</svg>`;
}
function axis(vw,xmin,xmax,ymin,ymax,yFormatter,integerX=false){
 const W=vw.w-vw.L-vw.R,H=vw.h-vw.T-vw.B;
 const xp=x=>vw.L+(x-xmin)/(xmax-xmin)*W;
 const yp=y=>vw.T+(ymax-y)/(ymax-ymin)*H;
 let s='';
 for(let i=0;i<=4;i++){
  const y=vw.T+i*H/4,val=ymax+(ymin-ymax)*i/4;
  s+=`<line x1="${vw.L}" x2="${vw.w-vw.R}" y1="${y}" y2="${y}" stroke="#E6EBF3" stroke-dasharray="3 5"/><text class="tick" x="${vw.L-10}" y="${y+4}" text-anchor="end">${yFormatter(val)}</text>`;
 }
 for(let i=0;i<=4;i++){
  const x=vw.L+i*W/4,val=xmin+(xmax-xmin)*i/4;
  s+=`<line x1="${x}" x2="${x}" y1="${vw.T}" y2="${vw.h-vw.B}" stroke="#E6EBF3" stroke-dasharray="3 5"/><text class="tick" x="${x}" y="${vw.h-vw.B+20}" text-anchor="middle">${integerX?Math.round(val):v(val,1)}</text>`;
 }
 s+=`<line x1="${vw.L}" x2="${vw.w-vw.R}" y1="${vw.h-vw.B}" y2="${vw.h-vw.B}" stroke="#A0ADBF"/>`;
 s+=`<text class="axis" x="${vw.w-vw.R}" y="${vw.h-3}" text-anchor="end">${integerX?'Індекс компоненти i':'Номер ітерації'}</text>`;
 return {xp,yp,s};
}
function line(arr,xp,yp,color,dashed=false,width=3){
 let d='',open=false;
 for(const p of arr){
  const x=xp(p.x),y=yp(p.y);
  if(!Number.isFinite(x)||!Number.isFinite(y)){open=false;continue;}
  d+=`${open?'L':'M'}${x.toFixed(2)} ${y.toFixed(2)} `;open=true;
 }
 return `<path d="${d}" fill="none" stroke="${color}" stroke-width="${width}" stroke-linejoin="round" stroke-linecap="round" ${dashed?'stroke-dasharray="7 5"':''}/>`;
}
/** Y — логарифм (десятковий) нев'язки; координати відновлюємо підписами 10^k. */
export function drawConvergence(node,jacobi,seidel,{small=false,partialJ=null,partialS=null}={}) {
 const vw=frame(small);
 const j=partialJ ?? jacobi.history, s=partialS ?? seidel.history;
 const ys=[...j,...s].map(p=>Math.log10(Math.max(Number.MIN_VALUE,p.residual)));
 const [ymin,ymax]=bound(ys,.10);
 const xEnd=Math.max(1,...j.map(p=>p.iteration),...s.map(p=>p.iteration));
 const {xp,yp,s:axes}=axis(vw,0,xEnd,ymin,ymax,z=>`10^${v(z,1)}`,false);
 const project=arr=>arr.map(z=>({x:z.iteration,y:Math.log10(Math.max(Number.MIN_VALUE,z.residual))}));
 let html=axes+line(project(j),xp,yp,C.jacobi)+line(project(s),xp,yp,C.seidel);
 for(const [arr,color] of [[j,C.jacobi],[s,C.seidel]]){
  const last=arr.at(-1);if(last)html+=`<circle cx="${xp(last.iteration)}" cy="${yp(Math.log10(Math.max(Number.MIN_VALUE,last.residual)))}" r="${small?4:6}" fill="${color}" stroke="white" stroke-width="2"><title>Ітерацій: ${last.iteration}, ||Ax−b||=${last.residual.toExponential(3)}</title></circle>`;
 }
 render(node,'Збіжність Якобі та Зейделя, логарифмічна шкала',vw,html);
}
/** Чотири криві, кожна за індексом компоненти i=1..n. */
export function drawSolutions(node,model,{small=false,shown=null}={}) {
 const vw=frame(small),{system,gauss,jacobi,seidel}=model,n=system.n;
 const j=shown?.jacobi??jacobi.x,s=shown?.seidel??seidel.x;
 const curves=[{key:'exact',vector:system.xStar},{key:'gauss',vector:gauss.x},{key:'jacobi',vector:j},{key:'seidel',vector:s}].filter(item=>item.vector);
 const yy=curves.flatMap(k=>k.vector);
 const [ymin,ymax]=bound(yy,.11);
 const {xp,yp,s:axes}=axis(vw,1,n,ymin,ymax,z=>v(z,1),true);
 let html=axes;
 for(const c of curves){
  const points=c.vector.map((y,i)=>({x:i+1,y}));
  html+=line(points,xp,yp,C[c.key],c.key==='exact',c.key==='exact'?small?3.0:4:small?2:2.4);
  points.forEach((p,i)=>{html+=`<circle cx="${xp(p.x)}" cy="${yp(p.y)}" r="${c.key==='exact'?(small?3.5:4.5):(small?2.0:2.6)}" fill="${C[c.key]}" stroke="white" stroke-width="${c.key==='exact'?1.6:0.8}"><title>${c.key}, x${i+1}=${v(p.y,5)}</title></circle>`;});
 }
 render(node,'Точний розв’язок і результати Гаусса, Якобі, Зейделя',vw,html);
}
export const chartColors=C;
