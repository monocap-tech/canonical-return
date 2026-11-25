// CR-ALP-064 fifth signed-prime moment: exact frozen trials, no new integral.
function audit(old,second,third,fourth){
const abs=x=>x<0n?-x:x;function gcd(a,b){a=abs(a);b=abs(b);while(b){let c=a%b;a=b;b=c;}return a||1n;}
function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){let [x,y]=a.split("/");return q(BigInt(x),BigInt(y));}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}let c=gcd(a,b);return [a/c,b/c];}
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>q(-a[0],a[1]),sub=(a,b)=>add(a,neg(b)),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]),sq=a=>mul(a,a),cmp=(a,b)=>a[0]*b[1]<b[0]*a[1]?-1:a[0]*b[1]>b[0]*a[1]?1:0;
const sum=(...a)=>a.reduce(add,q(0)),scale=(n,a)=>mul(q(n),a),render=x=>x[0]+"/"+x[1],dec=x=>Number(x[0])/Number(x[1]),dot=(a,b)=>sum(...a.map((x,i)=>mul(x,b[i]))),mid=a=>div(add(...a),q(2)),lo=a=>a.sort(cmp)[0],hi=a=>a.sort(cmp).at(-1);
const g=q("72/25"),s=q("21/8"),l=q("51/200"),u=q("1101/200"),lu=mul(l,u),ss=sq(s),zero=q(0);
function solve2(A,b){const [[a,c],[_,d]]=A,det=sub(mul(a,d),sq(c));if(cmp(a,zero)<=0||cmp(det,zero)<=0)throw Error("nonpositive Gram");return [div(sub(mul(d,b[0]),mul(c,b[1])),det),div(sub(mul(a,b[1]),mul(c,b[0])),det)];}
function build(v){
 const [G,M,N,T,S,K]=v,h=sub(mul(ss,G),N),j=sub(mul(ss,M),T),k=sub(mul(ss,N),S),r=add(mul(g,h),j),b=add(mul(g,j),k),c=add(mul(g,k),sub(mul(ss,T),zero));
 const A=[[add(mul(g,G),M),add(mul(g,M),N)],[add(mul(g,M),N),add(mul(g,N),T)]],d=[add(mul(g,N),T),add(mul(g,T),S)],e=mul(g,S);
 return {G,M,N,T,S,K,h,j,k,r,b,c,A,d,e};
}
const k="2",box=[old[k].true_G,old[k].true_M,second[k].true_N,third[k].true_T,fourth[k].true_S,old[k].K].map(a=>a.map(x=>q(x)));
const v=box.map(mid),a=build(v),{G,M,N,T,S,K,h,j,r,b,c,A,d,e}=a;
const f=sub(sub(mul(g,G),M),mul(lu,K)),base=div(sq(h),r),res=sub(j,div(mul(b,h),r)),gap=sub(f,base);if(cmp(gap,zero)<=0)throw Error("already pass");
const yu2=div(gap,res),yu1=div(sub(h,mul(b,yu2)),r);
const alpha=solve2(A,[G,M]),ad=solve2(A,d),baseV=dot([G,M],alpha),resV=sub(N,dot(d,alpha)),gapV=sub(K,baseV);if(cmp(gapV,zero)<=0)throw Error("already fail");
const z3=div(gapV,resV),z12=alpha.map((x,i)=>sub(x,mul(ad[i],z3)));
const freeze=x=>q(Math.round(dec(x)*100000),100000);
const Y=[yu1,yu2].map(freeze),Z=[...z12,z3].map(freeze);
const ac=[],re=[],supportlo=[],supporthi=[];
for(let i=0;i<64;i++){
 const a=build(box.map((x,j)=>x[(i>>j)&1])),{G,M,N,T,S,K,h,j,r,b,c,A,d,e}=a;
 const [x,y]=Y,[z0,z1,z2]=Z;
 const top=sum(sub(sub(mul(g,G),M),mul(lu,K)),neg(scale(2,sum(mul(h,x),mul(j,y)))),mul(r,sq(x)),scale(2,mul(b,mul(x,y))),mul(c,sq(y)));
 ac.push(div(top,sq(y)));
 const bottom=sub(sub(scale(2,dot([G,M,N],Z)),dot([z0,z1],A.map(row=>dot(row,[z0,z1])))),sum(scale(2,mul(z2,dot(d,[z0,z1]))),mul(e,sq(z2)),K));
 re.push(div(bottom,sq(z2)));
 // Positive weighted moment matrix for (s +/- p)dmu.
 for(const sign of [1,-1]){
   const m=[G,M,N,T,S],w=Array.from({length:4},(_,j)=>add(mul(s,m[j]),scale(sign,m[j+1])));
   const D=[[w[0],w[1]],[w[1],w[2]]],d=[w[2],w[3]],c=dot(d,solve2(D,d));
   if(sign===1)supportlo.push(sub(c,mul(s,S)));else supporthi.push(sub(mul(s,S),c));
 }
}
const acceptance=hi(ac),rejection=lo(re),support=[lo(supportlo),hi(supporthi)];
return {sector:"2",trial_upper:Y.map(render),trial_lower:Z.map(render),sufficient_acceptance_strict_above:render(acceptance),sufficient_rejection_strict_below:render(rejection),support_outer:support.map(render),decimal:{acceptance:dec(acceptance),rejection:dec(rejection),support:support.map(dec)},status:"FIFTH_MOMENT_NOT_EVALUATED"};
}
const fs=require('fs');const args=process.argv.slice(2);const r=audit(...args.slice(0,4).map(p=>JSON.parse(fs.readFileSync(p,'utf8'))));if(args[4])fs.writeFileSync(args[4],JSON.stringify(r,null,2)+'\n');else console.log(JSON.stringify(r,null,2));
