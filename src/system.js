/** Варіант №1, індекси i,j у формулах починаються з одиниці. */
export const sizes = [5, 10, 20];
export const dValues = [1, 0.1];

export function knownSolution(n) {
  return Array.from({length:n}, (_, z) => n-z); // x*_i = n-i+1, i=z+1
}
export function multiply(A, x) {
  return A.map(row => row.reduce((s, a, j) => s + a*x[j], 0));
}
export function generateSystem(n, d=1) {
  if (!sizes.includes(n)) throw new Error('Доступні розміри: 5, 10 або 20.');
  if (!(d > 0) || !Number.isFinite(d)) throw new Error('Параметр d має бути додатним.');
  const A = Array.from({length:n}, (_, i) => Array.from({length:n}, (_, j) => i===j ? 0 : Math.sin(i+j+2)));
  for (let i=0;i<n;i++) A[i][i]=A[i].reduce((s,a,j)=>s+(i===j?0:Math.abs(a)),0)+d;
  const xStar=knownSolution(n);
  return {A,b:multiply(A,xStar),xStar,n,d,generated:true};
}
export function validateSystem(input) {
  if (!input || !Array.isArray(input.A) || !Array.isArray(input.b)) throw new Error('Потрібен об’єкт із матрицею A та вектором b.');
  const n=input.A.length;
  if (!sizes.includes(n)) throw new Error('Матриця має розмір 5×5, 10×10 або 20×20.');
  if (input.b.length!==n || input.A.some(row=>!Array.isArray(row)||row.length!==n)) throw new Error('Несумісні розміри A та b.');
  if (input.A.some(row=>row.some(v=>typeof v!=='number'||!Number.isFinite(v))) || input.b.some(v=>typeof v!=='number'||!Number.isFinite(v))) throw new Error('Усі елементи A і b повинні бути скінченними числами.');
  const xStar = input.xStar===undefined||input.xStar===null ? null : input.xStar;
  if(xStar && (!Array.isArray(xStar)||xStar.length!==n||xStar.some(v=>typeof v!=='number'||!Number.isFinite(v)))) throw new Error('Вектор xStar має містити рівно n чисел.');
  return { A:input.A.map(r=>r.slice()), b:input.b.slice(), xStar:xStar ? xStar.slice():null, n, generated:false, d:null };
}
export function diagonalDominance(A, tolerance=1e-12) {
  const margins = A.map((row,i)=>Math.abs(row[i])-row.reduce((s,v,j)=>s+(j===i?0:Math.abs(v)),0));
  return {margins, strict:margins.every(v=>v>tolerance), weak:margins.every(v=>v>=-tolerance), minMargin:Math.min(...margins)};
}
export function norm2(vector) { return Math.hypot(...vector); }
export function residual(A,b,x) {
  return norm2(multiply(A,x).map((v,i)=>v-b[i]));
}
export function errorNorm(x,xStar) {return xStar ? norm2(x.map((v,i)=>v-xStar[i])):null;}

/**
 * Доказова тестова система з a11=0. Операція R1 <- R1 - c·R2
 * (одночасно для A та b) не змінює розв’язок системи.
 * Просте переставлення перших двох рядків не занулює a11 для sin(i+j).
 */
export function makeZeroPivotCase(base) {
  const data=validateSystem(base);
  if(Math.abs(data.A[1][0])<1e-14) throw new Error('Другий рядок не підходить для створення нульового головного елемента.');
  const c=data.A[0][0]/data.A[1][0];
  data.A[0]=data.A[0].map((v,j)=>v-c*data.A[1][j]);
  data.A[0][0]=0; // логічно нульова різниця, усуваємо дрібний шум double
  data.b[0] -= c*data.b[1];
  data.xStar = base.xStar ? base.xStar.slice():null;
  return data;
}

export function parseImportedText(text) {
  let data;
  try { data=JSON.parse(text); }
  catch { throw new Error('Файл має містити коректний JSON: {"A": [...], "b": [...], "xStar": [...]}'); }
  return validateSystem(data);
}
