#!/usr/bin/env node
const assert=require("node:assert/strict"),fs=require("node:fs");
const S=10n**340n,ns=[2,3,4,5,7,8,9];
const gcd=(a,b)=>b?gcd(b,a%b):a;
const floor=(a,b)=>{assert(b>0n);return a/b-(a<0n&&a%b!==0n?1n:0n)};
const ceil=(a,b)=>-floor(-a,b);
const mul=(x,y)=>{const z=x.flatMap(a=>y.map(b=>a*b));return [floor(z.reduce((a,b)=>a<b?a:b),S),ceil(z.reduce((a,b)=>a>b?a:b),S)]};
const reduced={},cases={};let terms=0,mutations=0;
for(const n of ns){
 let k=0,dy=1n;while(BigInt(n)>2n*dy){dy*=2n;k++}
 let g=gcd(BigInt(n),dy),a=BigInt(n)/g,b=dy/g;g=gcd(a-b,a+b);
 let p=(a-b)/g,q=(a+b)/g;
 const key=b===1n?String(a):String(a)+"/"+String(b);
 if(!reduced[key]){
  let u=[floor(p*S,q),ceil(p*S,q)];assert(20n*u[1]<7n*S&&u[0]*q<=p*S&&p*S<=u[1]*q);
  const v=mul(u,u);let t=u,l=0n,h=0n,rows=[];
  for(let i=0;i<400;i++){
   const d=2n*BigInt(i)+1n,num=p**d,den=q**d*d,x=floor(t[0],d),y=ceil(t[1],d);
   assert(den>0n&&x*den<=num*S&&num*S<=y*den,n+" "+i);
   assert((floor(num*S,den)+1n)*den>num*S);
   assert((ceil(num*S,den)-1n)*den<num*S);
   rows.push([String(x),String(y)]);l+=x;h+=y;terms++;mutations+=2;t=mul(t,v);
  }
  reduced[key]={y:[Number(a),Number(b)],p:Number(p),q:Number(q),rows,endpoints:[String(2n*l-1n),String(2n*h+1n)],width:String(2n*(h-l)+2n)};
 }
 const [l,u]=reduced[key].endpoints.map(BigInt),[l2,u2]=reduced["2"].endpoints.map(BigInt),L=l+BigInt(k)*l2,U=u+BigInt(k)*u2;
 cases[n]={k,reduced:key,endpoints:[String(L),String(U)],width:String(U-L)};
}
assert.equal(terms,2000);assert.equal(mutations,4000);
assert.deepEqual(Object.fromEntries(Object.entries(cases).map(([k,v])=>[k,v.width])),{"2":"802","3":"1596","4":"1604","5":"2406","7":"2406","8":"2406","9":"3210"});
let matched=false,log2Matched=false;
if(process.argv[2]){
 const t=JSON.parse(fs.readFileSync(process.argv[2],"utf8"));
 assert.equal(t.scale,String(S));assert.deepEqual(Object.keys(t.reduced),Object.keys(reduced));
 for(const [k,v] of Object.entries(reduced))for(const [f,x] of Object.entries(v))assert.deepEqual(t.reduced[k][f],x,k+" "+f);
 assert.deepEqual(t.cases,cases);matched=true;
}
if(process.argv[3]){
 const t=JSON.parse(fs.readFileSync(process.argv[3],"utf8"));
 assert.deepEqual(t.term_bounds,reduced["2"].rows);assert.deepEqual(t.actual_endpoints,reduced["2"].endpoints);log2Matched=true;
}
console.log(JSON.stringify({status:"INDEPENDENT_PRIME_LOG_REPLAY_PASS",terms,mutations,newTerms:1600,matched,log2Matched,reduced,cases},null,2));
