#!/usr/bin/env node
const assert=require("node:assert/strict"),fs=require("node:fs");
const S=10n**340n;
const pi=JSON.parse(fs.readFileSync(process.argv[2],"utf8")).cases.Pi.endpoints.map(BigInt);
const lg=JSON.parse(fs.readFileSync(process.argv[3],"utf8")).cases["9"].endpoints.map(BigInt);
const floor=(n,d)=>{assert(d>0n);return n/d-(n<0n&&n%d!==0n?1n:0n)},ceil=(n,d)=>-floor(-n,d),point=n=>[n*S,n*S],add=(a,b)=>[a[0]+b[0],a[1]+b[1]],neg=a=>[-a[1],-a[0]],min=z=>z.reduce((a,b)=>a<b?a:b),max=z=>z.reduce((a,b)=>a>b?a:b);
const mul=(a,b)=>{const z=a.flatMap(x=>b.map(y=>x*y));return [floor(min(z),S),ceil(max(z),S)]};
const div=(a,b)=>{assert(b[0]>0n);const z=a.flatMap(x=>b.map(y=>[x*S,y]));return [min(z.map(([n,d])=>floor(n,d))),max(z.map(([n,d])=>ceil(n,d)))]};
const mok=(a,b,c)=>a.every(x=>b.every(y=>c[0]*S<=x*y&&x*y<=c[1]*S)),dok=(a,b,c)=>b.every(y=>y>0n)&&a.every(x=>b.every(y=>c[0]*y<=x*S&&x*S<=c[1]*y));
let ops=0,mutations=0;const check=(f,a,b,c)=>{assert(f(a,b,c));assert(!f(a,b,[c[0]+1n,c[1]]));assert(!f(a,b,[c[0],c[1]-1n]));ops++;mutations+=2};
const M=(a,b)=>{const c=mul(a,b);check(mok,a,b,c);return c},D=(a,b)=>{const c=div(a,b);check(dok,a,b,c);return c};
const collar=D(point(11n),point(5n)),nu0=D(pi,collar),nu=M(point(85n),nu0),arg=M(nu,lg),two=mul(point(2n),pi),num=arg[0]+arg[1],den=two[0]+two[1];
let period=num/den;const rem=num%den;if(2n*rem>den||2n*rem===den&&period%2n!==0n)period++;
const periodic=M(point(2n*period),pi),reduced=add(arg,neg(periodic)),q=neg(M(reduced,reduced)),halfNu=D(nu,point(2n)),shift=D(point(2049n),point(4n)),atanArg=D(halfNu,shift),v=neg(M(atanArg,atanArg));
assert(reduced[0]>0n&&reduced[1]<4n*S&&atanArg[0]>0n&&3n*atanArg[1]<=S);assert.equal(period,42n);
const end=a=>a.map(String),traces={};
for(const [name,e,init] of [["Sin",1,reduced],["Cos",0,point(1n)],["Atan",-1,atanArg]]){
 let t=init,z=e<0?point(0n):t,ts=[end(t)],ms=[];
 for(let i=0;i<(e<0?400:399);i++){
  if(e<0){const term=D(t,point(BigInt(2*i+1)));ms.push(end(term));z=add(z,term);t=M(t,v)}
  else{const mid=M(t,q);ms.push(end(mid));t=D(mid,point(BigInt((2*i+e+1)*(2*i+e+2))));z=add(z,t)}
  ts.push(end(t));
 }
 traces[name]={t:ts,mid:ms,out:end([z[0]-1n,z[1]+1n]),width:String(z[1]-z[0]+2n)};
}
assert.equal(ops,2406);assert.equal(mutations,4812);
assert.deepEqual(Object.fromEntries(Object.entries(traces).map(([k,v])=>[k,v.width])),{Sin:"14402668",Cos:"14297535",Atan:"707"});
const data={scale:String(S),period:Number(period),pi:end(pi),log9:end(lg),collar:end(collar),nu0:end(nu0),nu:end(nu),arg:end(arg),periodic:end(periodic),reduced:end(reduced),q:end(q),halfNu:end(halfNu),shift:end(shift),atanArg:end(atanArg),v:end(v),traces};
let matched=false;
if(process.argv[4]){const a=JSON.parse(fs.readFileSync(process.argv[4],"utf8"));for(const [k,v] of Object.entries(data))assert.deepEqual(a[k],v,k);matched=true}
console.log(JSON.stringify({status:"INDEPENDENT_SIGNED_PILOT_PASS",ops,mutations,matched,...data},null,2));
