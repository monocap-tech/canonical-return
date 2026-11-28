#!/usr/bin/env node
const assert=require("node:assert/strict"),fs=require("node:fs");
const S=10n**340n;
const pi=JSON.parse(fs.readFileSync(process.argv[2],"utf8")).cases.Pi.endpoints.map(BigInt);
const log2=JSON.parse(fs.readFileSync(process.argv[3],"utf8")).actual_endpoints.map(BigInt);
const floor=(n,d)=>{assert(d>0n);return n/d-(n<0n&&n%d!==0n?1n:0n)},ceil=(n,d)=>-floor(-n,d);
const point=x=>[x*S,x*S],plus=(a,b)=>[a[0]+b[0],a[1]+b[1]],minus=(a,b)=>[a[0]-b[1],a[1]-b[0]];
const min=z=>z.reduce((a,b)=>a<b?a:b),max=z=>z.reduce((a,b)=>a>b?a:b);
const mul=(a,b)=>{const z=a.flatMap(x=>b.map(y=>x*y));return [floor(min(z),S),ceil(max(z),S)]};
const div=(a,b)=>{assert(b[0]>0n);const z=a.flatMap(x=>b.map(y=>[x*S,y]));return [min(z.map(([n,d])=>floor(n,d))),max(z.map(([n,d])=>ceil(n,d)))]};
const mok=(a,b,c)=>a.every(x=>b.every(y=>c[0]*S<=x*y&&x*y<=c[1]*S));
const dok=(a,b,c)=>b[0]>0n&&b[1]>0n&&a.every(x=>b.every(y=>c[0]*y<=x*S&&x*S<=c[1]*y));
let mutations=0,steps=0;
const check=(fn,a,b,c)=>{assert(fn(a,b,c));assert(!fn(a,b,[c[0]+1n,c[1]]));assert(!fn(a,b,[c[0],c[1]-1n]));mutations+=2;steps++};
assert(pi[0]+pi[1]>4n*S);
const x=div(pi,point(2n));assert(x[0]+x[1]>=2n*S&&x[0]+x[1]<=4n*S);
const u=div(minus(x,point(1n)),plus(x,point(1n))),v=mul(u,u);
assert(20n*u[1]<7n*S&&u[0]>0n);
check(dok,pi,point(2n),x);check(dok,minus(x,point(1n)),plus(x,point(1n)),u);check(mok,u,u,v);
let t=u,z=[0n,0n],powers=[t],terms=[];
for(let i=0;i<400;i++){
 const d=point(BigInt(2*i+1)),term=div(t,d),nt=mul(t,v);
 check(dok,t,d,term);check(mok,t,v,nt);z=plus(z,term);terms.push(term);powers.push(nt);t=nt;
}
const pilotOps=steps,red=[2n*z[0]-1n,2n*z[1]+1n],full=plus(red,log2);
const ints=[-7n,-3n,0n,2n,5n].map(x=>x*S/3n);
let badDen=0;
for(let i=0;i<ints.length;i++)for(let j=i;j<ints.length;j++)for(let k=0;k<ints.length;k++)for(let l=k;l<ints.length;l++){
 const a=[ints[i],ints[j]],b=[ints[k],ints[l]],c=mul(a,b);
 check(mok,a,b,c);
 if(b[0]>0n)check(dok,a,b,div(a,b));
 else{assert(!dok(a,b,[-100n*S,100n*S]));badDen++}
}
assert.equal(pilotOps,803);assert.equal(steps-pilotOps,270);assert.equal(mutations,2146);
assert.equal(full[1]-full[0],5588n);assert.equal(red[1]-red[0],4786n);
const end=a=>a.map(String);
const fields={scale:String(S),pi:end(pi),log2:end(log2),k:1,x:end(x),u:end(u),v:end(v),powers:powers.map(end),terms:terms.map(end),reduced:end(red),actual:end(full),width:String(full[1]-full[0])};
let matched=false;
if(process.argv[4]){
 const a=JSON.parse(fs.readFileSync(process.argv[4],"utf8"));
 for(const [k,v] of Object.entries(fields))assert.deepEqual(a[k],v,k);
 matched=true;
}
console.log(JSON.stringify({status:"INDEPENDENT_LOG_PI_REPLAY_PASS",pilotOps,syntheticOps:steps-pilotOps,mutations,badDen,matched,...fields},null,2));
