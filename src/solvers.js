import {diagonalDominance,errorNorm,multiply,norm2,residual,validateSystem} from './system.js';
const clock=()=>performance.now();
const DEFAULT_EPS=1e-10;
const DEFAULT_MAX=20000;

export function gaussWithoutPivot(input) {
  const {A,b,n}=validateSystem(input);const M=A.map(row=>row.slice());const y=b.slice();
  for(let k=0;k<n;k++) {
    if(Math.abs(M[k][k])<1e-14) throw new Error(`Метод Гаусса без вибору: нульовий головний елемент у рядку ${k+1}.`);
    for(let i=k+1;i<n;i++) {
      const f=M[i][k]/M[k][k];
      for(let j=k+1;j<n;j++) M[i][j]-=f*M[k][j];
      M[i][k]=0; y[i]-=f*y[k];
    }
  }
  return backSubstitute(M,y);
}
function backSubstitute(A,b,tolerance=1e-13) {
  const n=b.length, x=Array(n).fill(0);
  for(let i=n-1;i>=0;i--) {
    const scale=Math.max(1,...A[i].map(Math.abs));
    if(Math.abs(A[i][i])<=tolerance*scale) throw new Error('Вироджена або чисельно вироджена матриця: неможливий зворотний хід.');
    let z=b[i];for(let j=i+1;j<n;j++)z-=A[i][j]*x[j];
    x[i]=z/A[i][i];
  }
  return x;
}
/** Частковий вибір головного елемента за стовпцем з перевіркою масштабу. */
export function gaussPivot(input, {trace=false}={}) {
  const {A,b,xStar,n}=validateSystem(input);const start=clock();
  const M=A.map(r=>r.slice());const rhs=b.slice(); const steps=[];
  let swaps=0;let operations=0;
  function save(type,k,selectedRow=null,note=''){
    if(trace)steps.push({type,k,selectedRow, note, A:M.map(r=>r.slice()), b:rhs.slice()});
  }
  save('start',0,null,'Початкова розширена матриця [A|b]');
  for(let k=0;k<n;k++) {
    let pivot=k;
    for(let i=k+1;i<n;i++)if(Math.abs(M[i][k])>Math.abs(M[pivot][k]))pivot=i;
    save('pivot',k,pivot,`Пошук максимального |aᵢ,${k+1}| у стовпці ${k+1}`);
    const maxCell=Math.abs(M[pivot][k]);
    const scale=Math.max(...M.map(row=>Math.max(...row.map(Math.abs))));
    if(maxCell<=1e-13*Math.max(scale,Number.MIN_VALUE))throw new Error(`Вироджена або майже вироджена матриця: стовпець ${k+1}.`);
    if(pivot!==k) {
      [M[k],M[pivot]]=[M[pivot],M[k]];
      [rhs[k],rhs[pivot]]=[rhs[pivot],rhs[k]];
      swaps++;
      save('swap',k,k,`Переставляю рядки ${k+1} та ${pivot+1}`);
    }
    for(let i=k+1;i<n;i++) {
      const factor=M[i][k]/M[k][k];
      if(factor===0)continue;
      M[i][k]=0;
      for(let j=k+1;j<n;j++) M[i][j]-=factor*M[k][j];
      rhs[i]-=factor*rhs[k];operations++;
    }
    save('eliminate',k,k,`Занулюю елементи під головним у стовпці ${k+1}`);
  }
  const x=backSubstitute(M,rhs);
  const timeMs=clock()-start;
  save('result',n-1,n-1,'Зворотний хід завершено');
  return {method:'Гаусс',x,iterations:0,steps,timeMs,swaps,operations,
    residual:residual(A,b,x), error:errorNorm(x,xStar),converged:true,
    reason:'Розв’язок знайдено прямим методом',history:[]};
}

export function iterativeSolve(input, method, {eps=DEFAULT_EPS,maxIterations=DEFAULT_MAX, traceEvery=1}={}) {
  const {A,b,xStar,n}=validateSystem(input);
  if(!['jacobi','seidel'].includes(method)) throw new Error('Невідомий ітераційний метод.');
  if(!(eps>0) || !Number.isFinite(eps))throw new Error('ε має бути додатним скінченним числом.');
  if(!Number.isInteger(maxIterations)||maxIterations<1||maxIterations>100000)throw new Error('Кількість ітерацій має бути в межах 1–100000.');
  if(A.some((r,i)=>Math.abs(r[i])<1e-14))throw new Error(`Нульовий діагональний елемент: ${method==='jacobi'?'Якобі':'Зейдель'} неможливий без перестановки рівнянь.`);
  const start=clock(), name=method==='jacobi'?'Якобі':'Зейдель';
  const dominance=diagonalDominance(A);
  const history=[{iteration:0,residual:residual(A,b,Array(n).fill(0)),x:Array(n).fill(0)}];
  let x=Array(n).fill(0), iterations=0, r=history[0].residual;
  let reason='Досягнуто ліміту ітерацій', converged=r<=eps;
  if(converged)reason='Початкове наближення задовольняє критерій';
  for(let k=1; k<=maxIterations && !converged; k++) {
    const old=x.slice();
    if(method==='jacobi') {
      const next=Array(n).fill(0);
      for(let i=0;i<n;i++) {
        let s=b[i];for(let j=0;j<n;j++)if(j!==i)s-=A[i][j]*old[j];
        next[i]=s/A[i][i];
      }
      x=next;
    } else {
      for(let i=0;i<n;i++) {
        let s=b[i];
        for(let j=0;j<n;j++)if(j!==i)s-=A[i][j]*(j<i?x[j]:old[j]);
        x[i]=s/A[i][i];
      }
    }
    r=residual(A,b,x);iterations=k;
    if(k%Math.max(1,traceEvery)===0 || r<=eps || k===maxIterations)history.push({iteration:k,residual:r,x:x.slice()});
    if(!Number.isFinite(r) || x.some(v=>!Number.isFinite(v)) || r>1e100) {
      reason='Розбіжність або переповнення';break;
    }
    if(r<=eps){converged=true;reason='Норма нев’язки ≤ ε';break;}
    if(k>=5000 && Math.abs(old.reduce((s,v,i)=>s+Math.abs(v-x[i]),0))===0){reason='Стагнація чисел double';break;}
  }
  const timeMs=clock()-start;
  return {method:name,x,iterations,timeMs,converged,reason,residual:r,error:errorNorm(x,xStar),history,dominance};
}
export const jacobi=(input,opts={})=>iterativeSolve(input,'jacobi',opts);
export const seidel=(input,opts={})=>iterativeSolve(input,'seidel',opts);
export function analyze(input,opts={}) {
 const data=validateSystem(input);
 const dominance=diagonalDominance(data.A);
 const gauss=gaussPivot(data,opts.gauss||{});
 const j=jacobi(data,opts.iteration||{}),s=seidel(data,opts.iteration||{});
 return {gauss,jacobi:j,seidel:s,dominance,system:data};
}
