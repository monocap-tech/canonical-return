#!/usr/bin/env node
const assert=require("node:assert/strict"),fs=require("node:fs");
const S=10n**340n;
const floor=(n,d)=>{assert(d>0n);return n/d-(n<0n&&n%d!==0n?1n:0n)};
const ceil=(n,d)=>-floor(-n,d),iv=(l,h=l)=>[l,h],add=(x,y)=>iv(x[0]+y[0],x[1]+y[1]),neg=x=>iv(-x[1],-x[0]),mul=(x,y)=>{let z=x.flatMap(a=>y.map(b=>a*b));return iv(floor(z.reduce((a,b)=>a<b?a:b),S),ceil(z.reduce((a,b)=>a>b?a:b),S))},divn=(x,d)=>iv(floor(x[0],d),ceil(x[1],d)),scale=(x,k)=>k>=0n?iv(x[0]*k,x[1]*k):iv(x[1]*k,x[0]*k);
const cases=[["Atan5","atan",1,5],["Atan239","atan",1,239],["ExpA","exp",11,10],["ExpHalfA","exp",11,20]],res={};
let terms=0,mutations=0;
for(const [key,kind,p0,q0] of cases){
 const p=BigInt(p0),q=BigInt(q0),x=iv(floor(p*S,q),ceil(p*S,q));
 assert(kind==="atan"?3n*x[1]<=S:x[0]>=0n&&x[1]<3n*S);
 let t=kind==="atan"?x:iv(S),v=neg(mul(x,x)),z=iv(0n),fac=1n,rows=[];
 for(let i=0;i<400;i++){
  const j=BigInt(i);if(kind==="exp"&&i){fac*=j;t=divn(mul(t,x),j)}
  const d=kind==="atan"?2n*j+1n:j,num=kind==="atan"?(-1n)**j*p**d:p**j,den=kind==="atan"?q**d*d:q**j*fac,term=kind==="atan"?divn(t,d):t;
  assert(den>0n&&term[0]*den<=num*S&&num*S<=term[1]*den,key+" "+i);terms++;
  assert((floor(num*S,den)+1n)*den>num*S);assert((ceil(num*S,den)-1n)*den<num*S);mutations+=2;
  rows.push(term);z=add(z,term);if(kind==="atan")t=mul(t,v);
 }
 res[key]={kind,p:p0,q:q0,terms:400,rows:rows.map(r=>r.map(String)),endpoints:[String(z[0]-1n),String(z[1]+1n)],width:String(z[1]-z[0]+2n)};
}
const pi=add(scale(res.Atan5.endpoints.map(BigInt),16n),scale(res.Atan239.endpoints.map(BigInt),-4n));
assert(pi[0]>3n*S&&pi[1]*7n<22n*S);res.Pi={endpoints:pi.map(String),width:String(pi[1]-pi[0])};
assert(terms===1600&&mutations===3200);
const widths={Atan5:"398",Atan239:"402",ExpA:"402",ExpHalfA:"400",Pi:"7976"};
for(const [k,v] of Object.entries(res))assert.equal(v.width,widths[k]);
let matched=false;
if(process.argv[2]){
 const t=JSON.parse(fs.readFileSync(process.argv[2],"utf8"));assert.equal(t.scale,String(S));assert.deepEqual(Object.keys(t.cases),Object.keys(res));
 for(const [k,v] of Object.entries(res))for(const [f,x] of Object.entries(v))assert.deepEqual(t.cases[k][f],x,k+" "+f);
 matched=true;
}
console.log(JSON.stringify({status:"INDEPENDENT_CONSTANT_REPLAY_PASS",terms,mutations,matched,cases:res},null,2));
