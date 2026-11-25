// Independent exact rational sector replay; no physical quadrature regeneration.
function audit(old,second,targets,logs,k){
const abs=x=>x<0n?-x:x;function gcd(a,b){a=abs(a);b=abs(b);while(b){const t=a%b;a=b;b=t;}return a||1n;}
function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){const v=a.split("/");return q(BigInt(v[0]),BigInt(v[1]));}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}const c=gcd(a,b);return [a/c,b/c];}
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>q(-a[0],a[1]),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),inv=a=>q(a[1],a[0]),cmp=(a,b)=>a[0]*b[1]<b[0]*a[1]?-1:a[0]*b[1]>b[0]*a[1]?1:0;
const I=(a,b=a)=>[q(a),q(b)],pt=a=>[a,a],plus=(a,b)=>[add(a[0],b[0]),add(a[1],b[1])],minus=(a,b)=>plus(a,[neg(b[1]),neg(b[0])]),times=(a,b)=>{const v=a.flatMap(x=>b.map(y=>mul(x,y))).sort(cmp);return [v[0],v[3]];},divide=(a,b)=>{if(b[0][0]<=0n)throw Error("nonpositive denominator");return times(a,[inv(b[1]),inv(b[0])]);},sq=a=>{const p=mul(a[0],a[0]),q=mul(a[1],a[1]);return [a[0][0]<=0n&&a[1][0]>=0n?[0n,1n]:cmp(p,q)<0?p:q,cmp(p,q)>0?p:q];};
const check=(b,msg)=>{if(!b)throw Error(msg);},lt=(a,b)=>cmp(a,b)<0,le=(a,b)=>cmp(a,b)<=0,overlap=(a,b)=>le(a[0],b[1])&&le(b[0],a[1]);
const rows={};for(const line of logs.trim().split("\n")){const [name,l,h,sep,den]=line.split(" ");check(sep==="/"&&den==="1e15"&&!rows[name],"log format");rows[name]=[q(l,10n**15n),q(h,10n**15n)];check(le(rows[name][0],rows[name][1]),"endpoint order");}
check(Object.keys(rows).length===14,"row count");
// independent rational weights log(p)/sqrt(n): 96-term atanh and integer sqrt
function isqrt(n){let x=n,y=(x+1n)/2n;while(y<x){x=y;y=(x+n/x)/2n;}return x;}
function logbounds(p){const z=q(p-1,p+1),zz=mul(z,z);let v=z,sum=q(0);for(let i=0;i<96;i++){sum=add(sum,mul(q(2,2*i+1),v));v=mul(v,zz);}const tail=mul(mul(q(2,193),v),inv(add(q(1),neg(zz))));return [sum,add(sum,tail)];}
const weights=[2,3,4,5,7,8,9].map((n,i)=>{const p=[2,3,2,5,7,2,3][i],D=10n**40n,b=isqrt(BigInt(n)*D*D);return divide(logbounds(p),[q(b,D),q(b+1n,D)]);});
const ws=weights.reduce(plus,I(0));check(lt(mul(q(2),ws[1]),q(8)),"weight bound");
const recomputed=weights.reduce((a,w,i)=>plus(a,times(I(-2),times(w,rows["third_corr"+i]))),I(0));check(overlap(rows.T_proxy,recomputed),"prime correlation overlap");
const s=I("21/8"),g=I("72/25"),ell=I("51/200"),u=I("1101/200"),G=I(...old[k].true_G),M=I(...old[k].true_M),N=I(...second[k].true_N),K=I(...old[k].K),eta=pt(rows.projection_error[1]),eps=I(old[k].proxy_error_upper);
check(le(q(0),rows.N_proxy[0])&&lt(rows.N_proxy[1],q("121/400"))&&overlap(rows.N_proxy,I(...second[k].proxy_N)),"second moment consistency");
check(le(rows.projection_error[0],q(0))&&le(q(0),eta[0])&&lt(eta[1],q("1/10")),"projection uncertainty");
const shift=plus(times(I(2n*10n**22n),I("1/"+10n**70n)),times(I(2n*10n**19n),I("32/"+10n**35n)));check(lt(times(I(8),shift)[1],q("1/1000000000000")),"shift budget");
const delta=plus(plus(times(s,eps),I("1/10000000000000")),eta),norm=plus(I("11/20"),eta),err=plus(times(s,plus(times(I(2),times(norm,delta)),sq(delta))),I("1/1000000000000"));check(lt(norm[1],q(1)),"proxy norm");
const T=[add(rows.T_proxy[0],neg(err[1])),add(rows.T_proxy[1],err[1])];
const h=minus(times(sq(s),G),N),c=minus(plus(times(g,h),times(sq(s),M)),T),m=plus(times(g,G),M),d=plus(times(g,M),N),e=plus(times(g,N),T),schur=minus(e,divide(sq(d),m));
check(h[0][0]>0n&&c[0][0]>0n&&schur[0][0]>0n,"nondegenerate moments");
const U=divide(minus(minus(times(g,G),M),divide(sq(h),c)),times(ell,u)),V=plus(divide(sq(G),m),divide(sq(minus(M,divide(times(d,G),m))),schur));
const status=lt(V[1],K[0])&&lt(K[1],U[0])?"FOUR_SCALAR_MOMENTS_INSUFFICIENT":lt(U[1],K[0])?"PASS_WITNESS":lt(K[1],V[0])?"FAIL_COMPARISON_WITNESS":"UNDECIDED_INTERVAL_WIDTH";
const render=a=>a.map(x=>x[0]+"/"+x[1]);
return {sector:k,status,weighted_correlation_consistency:"PASS",transport:render(err),true_T:render(T),sharp_lower:render(V),sharp_upper:render(U),K:render(K)};
}
const fs=require('fs');const [p,q,t,l,k,out]=process.argv.slice(2);const read=x=>JSON.parse(fs.readFileSync(x,'utf8'));const result=audit(read(p),read(q),read(t),fs.readFileSync(l,'utf8'),k);console.log(JSON.stringify(result,null,2));if(out)fs.writeFileSync(out,JSON.stringify(result,null,2)+'\n');
