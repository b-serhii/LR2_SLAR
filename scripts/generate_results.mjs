import {generateSystem,makeZeroPivotCase} from '../src/system.js';
import {analyze,gaussWithoutPivot,gaussPivot} from '../src/solvers.js';
import {writeFileSync} from 'node:fs';
const t=[];const md=[];
for(const n of [5,10,20])for(const d of [1,.1]){
 const system=generateSystem(n,d),{gauss,jacobi,seidel,dominance}=analyze(system,{iteration:{eps:1e-10,maxIterations:20000}});
 const count={};for(const r of [gauss,jacobi,seidel]){
  // 100 repeated runs provide more reproducible timing compared with a single micro-call.
  // Standard clock includes calculations; history and outcomes are the full specified algorithm.
  const times=[];
  for(let i=0;i<40;i++) {
   const rr=analyze(system,{iteration:{eps:1e-10,maxIterations:20000}})[r.method==='Гаусс'?'gauss':r.method==='Якобі'?'jacobi':'seidel'];
   times.push(rr.timeMs);
  }
  times.sort((a,b)=>a-b);const median=(times[19]+times[20])/2;
  const row={n,d,method:r.method,iterations:r.iterations,residual:r.residual,error:r.error,medianTimeMs:median,converged:r.converged,margin:dominance.minMargin};
  t.push(row);count[r.method]=r.iterations;
 }
 console.log('n',n,'d',d,'iter',count);
}
const cols=['n','d','method','iterations','residual','error','medianTimeMs','converged','margin'];
const csv=[cols.join(';'),...t.map(r=>cols.map(c=>r[c]).join(';'))].join('\n');
writeFileSync(new URL('../docs/analytical_results.csv',import.meta.url),'\uFEFF'+csv,'utf8');
writeFileSync(new URL('../docs/analytical_results.json',import.meta.url),JSON.stringify(t,null,2),'utf8');
const s=makeZeroPivotCase(generateSystem(5,1));let error='';try{gaussWithoutPivot(s);}catch(e){error=e.message;}
const fx=gaussPivot(s,{trace:true});
writeFileSync(new URL('../docs/zero_pivot_demonstration.txt',import.meta.url),`Нульовий головний елемент a11 = ${s.A[0][0]}\nГаусс без вибору: ${error}\nГаусс із вибором: перестановок ${fx.swaps}, нев’язка ${fx.residual}, похибка ${fx.error}\nДля sin(i+j) переставлення перших двох рядків не гарантує нуль a11, тому застосовано еквівалентне перетворення R1 <- R1 - (a11/a21) R2 та b1 аналогічно.\n`);
