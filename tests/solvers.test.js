import test from 'node:test';
import assert from 'node:assert/strict';
import {generateSystem,diagonalDominance,makeZeroPivotCase,multiply,residual,validateSystem,parseImportedText} from '../src/system.js';
import {gaussPivot,gaussWithoutPivot,jacobi,seidel,analyze} from '../src/solvers.js';
const close=(a,b,tol=1e-7)=>assert.ok(Math.abs(a-b)<tol,`Очікувалося ${b}, отримано ${a}`);
for (const n of [5,10,20]) {
 for(const d of [1,.1]) {
  test(`Варіант 1, n=${n}, d=${d}: генератор та строга діагональна перевага`,()=>{
   const s=generateSystem(n,d);
   for(let i=0;i<n;i++){
    close(s.xStar[i],n-i,1e-12);
    for(let j=0;j<n;j++)if(i!==j)close(s.A[i][j],Math.sin(i+j+2),1e-12);
   }
   const diag=diagonalDominance(s.A);
   assert.equal(diag.strict,true);
   close(diag.minMargin,d,1e-11);
   close(residual(s.A,s.b,s.xStar),0,1e-11);
  });
  test(`Гаусс з вибором головного: n=${n}, d=${d}`,()=>{
   const s=generateSystem(n,d),r=gaussPivot(s);
   assert.ok(r.converged);assert.ok(r.residual<1e-10);
   assert.ok(r.error<1e-9); assert.equal(r.iterations,0);
  });
  test(`Якобі: n=${n}, d=${d}`,()=>{
   const s=generateSystem(n,d),r=jacobi(s,{eps:1e-10,maxIterations:20000});
   assert.equal(r.converged,true,r.reason);
   assert.ok(r.residual<=1.01e-10);assert.ok(r.error<1e-8);
   assert.ok(r.history.length>1);
  });
  test(`Зейдель: n=${n}, d=${d}`,()=>{
   const s=generateSystem(n,d),r=seidel(s,{eps:1e-10,maxIterations:20000});
   assert.equal(r.converged,true,r.reason);
   assert.ok(r.residual<=1.01e-10);assert.ok(r.error<1e-8);
  });
 }
}
test('Перевірка конкретної матриці з a11=0: Гаусс без вибору падає, з вибором працює',()=>{
 const orig=generateSystem(5,1);const z=makeZeroPivotCase(orig);
 assert.equal(z.A[0][0],0);
 assert.throws(()=>gaussWithoutPivot(z),/нульовий головний/);
 const fixed=gaussPivot(z,{trace:true});
 assert.ok(fixed.swaps>=1);
 assert.ok(fixed.error<1e-9);
 assert.ok(fixed.steps.some(x=>x.type==='pivot') && fixed.steps.some(x=>x.type==='swap'));
});
test('Вироджена матриця: повторювані рядки призводять до повідомлення',()=>{
 const s=generateSystem(5,1);s.A[1]=s.A[0].slice();s.b[1]=s.b[0];
 assert.throws(()=>gaussPivot(s),/Вироджена/);
});
test('Відсутність діагональної переваги виявляється',()=>{
 const s=generateSystem(5,1);s.A[0][0]=0.5;
 assert.equal(diagonalDominance(s.A).strict,false);
});
test('Нульовий діагональний елемент: ітераційні методи повідомляють про помилку',()=>{
 const s=generateSystem(5,1);s.A[0][0]=0;
 assert.throws(()=>jacobi(s),/Нульовий діагональний/);
 assert.throws(()=>seidel(s),/Нульовий діагональний/);
});
test('Надто малий maxIterations: повертається неконвергентний результат',()=>{
 const s=generateSystem(20,0.1),r=jacobi(s,{maxIterations:1,eps:1e-14});
 assert.equal(r.converged,false);assert.equal(r.iterations,1);
});
test('Валідація розмірів і скінченності',()=>{
 assert.throws(()=>validateSystem({A:[[1]],b:[1]}),/5×5/);
 const s=generateSystem(5,1);s.A[0][0]=Infinity;
 assert.throws(()=>validateSystem(s),/скінченними/);
});
test('JSON імпорт працює і зберігає відсутність точного розв’язку',()=>{
 const s=generateSystem(5,1);delete s.xStar;
 const v=parseImportedText(JSON.stringify(s));
 assert.equal(v.xStar,null);
});
test('Властивість: зміна одного рядка, що зберігає СЛАР, не змінює x*',()=>{
 const s=generateSystem(10,1),z=makeZeroPivotCase(s);
 const r=multiply(z.A,s.xStar).map((v,i)=>v-z.b[i]);
 assert.ok(Math.hypot(...r)<1e-10);
});
test('Неоптимальний критерій зупинки за зміною x може приховати велику нев’язку',()=>{
 const n=5,A=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1e12:1e11));
 const b=Array(n).fill(1e7),input={A,b};
 const first=Array(n).fill(1e-5);
 const dx=Math.max(...first.map(Math.abs));
 assert.ok(dx<1e-4,'Наївний критерій за приростом оголосив би збіжність');
 assert.ok(residual(A,b,first)>1e6,'А фактична нев’язка все ще дуже велика');
 const r=jacobi(input,{eps:1e-4,maxIterations:20000});
 assert.ok(r.converged);
 assert.ok(r.residual<=1e-4);
 assert.ok(r.iterations>1);
});
