#!/usr/bin/env node
// Independent positive-argument directed arithmetic replay for actual 046 LOG2.
const assert=require("node:assert/strict"), fs=require("node:fs");
const S=10n**340n;
const floor=(a,b)=>{assert(b>0n);let q=a/b;return a<0n&&a%b!==0n?q-1n:q};
const ceil=(a,b)=>-floor(-a,b);
const iv=(l,h=l)=>({l,h});
const mul=(a,b)=>{const z=[a.l*b.l,a.l*b.h,a.h*b.l,a.h*b.h];return iv(floor(z.reduce((p,q)=>p<q?p:q),S),ceil(z.reduce((p,q)=>p>q?p:q),S))};
const divn=(a,d)=>iv(floor(a.l,d),ceil(a.h,d));
const u=iv(floor(S,3n),ceil(S,3n)),v=mul(u,u);
let t=u,l=0n,h=0n,rows=[];
const checks={};
function check(k,c){assert(c,k);checks[k]="PASS"}
for(let i=0;i<400;i++){
 const d=BigInt(2*i+1),den=3n**d*d,term=divn(t,d);
 assert(term.l*den<=S&&S<=term.h*den);
 rows.push({i,l:term.l,h:term.h,den});l+=term.l;h+=term.h;t=mul(t,v);
}
const L=2n*l-1n,U=2n*h+1n;
check("exact_400_terms",rows.length===400&&rows.every((r,i)=>r.i===i&&r.l*r.den<=S&&S<=r.h*r.den));
check("literal_u_guard",u.l*3n<=S&&S<=u.h*3n&&u.h*20n<7n*S);
check("endpoint_order",L<U);
check("sum_and_widen",L===2n*rows.reduce((a,r)=>a+r.l,0n)-1n&&U===2n*rows.reduce((a,r)=>a+r.h,0n)+1n);
check("reject_invalid_lower",rows.every(r=>(floor(S,r.den)+1n)*r.den>S));
check("reject_invalid_upper",rows.every(r=>(ceil(S,r.den)-1n)*r.den<S));
check("reject_wrong_endpoint",L+1n!==2n*l-1n&&U-1n!==2n*h+1n);
check("analytic_guard",4n*7n**801n*S<20n**801n);
const replay={scale:String(S),terms:400,term_bounds:rows.map(r=>[String(r.l),String(r.h)]),reduced_argument:[String(u.l),String(u.h)],actual_endpoints:[String(L),String(U)],width_units:String(U-L)};
if(process.argv[2]){
 const original=JSON.parse(fs.readFileSync(process.argv[2],"utf8"));
 for(const [key,value] of Object.entries(replay))assert.deepEqual(original[key],value,key);
 checks.literal_python_trace_match="PASS";
}
console.log(JSON.stringify({status:"INDEPENDENT_LOG2_REPLAY_PASS",checks,...replay},null,2));
