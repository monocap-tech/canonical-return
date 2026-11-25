// Exact rational fourth-moment feasibility and decision boundaries; no new integral.
function audit(old,second,third){
 const abs=x=>x<0n?-x:x;
 function gcd(a,b){a=abs(a);b=abs(b);while(b){const t=a%b;a=b;b=t;}return a||1n;}
 function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){const v=a.split("/");return q(BigInt(v[0]),BigInt(v[1]));}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}const c=gcd(a,b);return [a/c,b/c];}
 const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>q(-a[0],a[1]),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),inv=a=>q(a[1],a[0]),cmp=(a,b)=>a[0]*b[1]<b[0]*a[1]?-1:a[0]*b[1]>b[0]*a[1]?1:0;
 const I=(a,b=a)=>[q(a),q(b)],plus=(a,b)=>[add(a[0],b[0]),add(a[1],b[1])],minus=(a,b)=>plus(a,[neg(b[1]),neg(b[0])]),times=(a,b)=>{const v=a.flatMap(x=>b.map(y=>mul(x,y))).sort(cmp);return [v[0],v[3]];},divide=(a,b)=>{if(b[0][0]<=0n)throw Error("nonpositive denominator");return times(a,[inv(b[1]),inv(b[0])]);},sq=a=>{const p=mul(a[0],a[0]),r=mul(a[1],a[1]);return [a[0][0]<=0n&&a[1][0]>=0n?q(0):cmp(p,r)<0?p:r,cmp(p,r)>0?p:r];};
 const g=I("72/25"),s=I("21/8"),ell=I("51/200"),u=I("1101/200"),gs=times(g,s),ss=sq(s);
 function solve(G,M,N,T,K){
 const A=plus(plus(times(gs,G),times(u,M)),N),B=plus(plus(times(gs,M),times(u,N)),T),C0=plus(times(gs,N),times(u,T)),b0=plus(times(s,G),M),b1=plus(times(s,M),N);
 const D=minus(minus(times(gs,G),times(ell,M)),N),E=minus(minus(times(gs,M),times(ell,N)),T),F0=minus(times(gs,N),times(ell,T)),d0=minus(times(s,G),M),d1=minus(times(s,M),N);
 const U=minus(divide(G,ell),divide(sq(b0),times(ell,A))),V=plus(divide(G,u),divide(sq(d0),times(u,D)));
 const Ares=minus(b1,divide(times(B,b0),A)),Dres=minus(d1,divide(times(E,d0),D));
 const acc=plus(minus(divide(sq(B),A),C0),divide(sq(Ares),times(ell,minus(U,K))));
 const rej=minus(minus(F0,divide(sq(E),D)),divide(sq(Dres),times(u,minus(K,V))));
 const low=plus(divide(sq(N),G),divide(sq(minus(T,divide(times(M,N),G))),minus(N,divide(sq(M),G))));
 const high=minus(times(ss,N),divide(sq(minus(times(ss,M),T)),minus(times(ss,G),N)));
 return {accept_boundary:acc,reject_boundary:rej,support_lower_boundary:low,support_upper_boundary:high,upper_schur_zero:minus(divide(sq(B),A),C0),lower_schur_zero:minus(F0,divide(sq(E),D))};
 }
 const out={};for(const k of ["1","2"]){const G=I(...old[k].true_G),M=I(...old[k].true_M),N=I(...second[k].true_N),T=I(...third[k].true_T),K=I(...old[k].K);const r=solve(G,M,N,T,K);
 const render=a=>a.map(x=>x[0]+"/"+x[1]);
 out[k]=Object.fromEntries(Object.entries(r).map(([name,v])=>[name,render(v)]));
 }
 return out;
}
const fs=require('fs');const out=audit(...process.argv.slice(2,5).map(p=>JSON.parse(fs.readFileSync(p,'utf8'))));if(process.argv[5])fs.writeFileSync(process.argv[5],JSON.stringify(out,null,2)+'\n');else console.log(JSON.stringify(out,null,2));
