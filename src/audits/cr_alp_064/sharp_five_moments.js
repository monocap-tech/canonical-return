// CR-ALP-064: independent exact transport and optimal five-moment discrimination.
function audit(old,second,third,fourth,targets,support,logs){
const abs=x=>x<0n?-x:x;function gcd(a,b){a=abs(a);b=abs(b);while(b){let c=a%b;a=b;b=c;}return a||1n;}
function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){let [x,y]=a.split("/");return q(BigInt(x),BigInt(y));}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}let c=gcd(a,b);return [a/c,b/c];}
const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>q(-a[0],a[1]),sub=(a,b)=>add(a,neg(b)),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]),sq=a=>mul(a,a),cmp=(a,b)=>a[0]*b[1]<b[0]*a[1]?-1:a[0]*b[1]>b[0]*a[1]?1:0;
const lt=(a,b)=>cmp(a,b)<0,le=(a,b)=>cmp(a,b)<=0,mi=(...a)=>a.sort(cmp)[0],ma=(...a)=>a.sort(cmp).at(-1),sum=(...x)=>x.reduce(add,q(0)),scale=(n,x)=>mul(q(n),x),render=x=>x[0]+"/"+x[1],pair=a=>a.map(render),check=(b,msg)=>{if(!b)throw Error(msg);},overlap=(a,b)=>le(a[0],b[1])&&le(b[0],a[1]),same=(a,b)=>cmp(a,b)===0;
const zero=q(0),one=q(1),g=q("72/25"),s=q("21/8"),l=q("51/200"),u=q("1101/200"),gs=mul(g,s);
function isqrt(n){if(n<2n)return n;let x=n,y=(x+1n)/2n;while(y<x){x=y;y=(x+n/x)/2n;}return x;}
function sqrtup(x){const D=10n**30n;return q(isqrt(x[0]*D*D/x[1])+1n,D);}
function gram(v,side){const [G,M,N,T,S]=v;return side===0?{a:sum(mul(gs,G),mul(u,M),N),b:sum(mul(gs,M),mul(u,N),T),c:sum(mul(gs,N),mul(u,T),S),d:add(mul(s,G),M),e:add(mul(s,M),N),endpoint:l,sign:-1}:{a:sub(sub(mul(gs,G),mul(l,M)),N),b:sub(sub(mul(gs,M),mul(l,N)),T),c:sub(sub(mul(gs,N),mul(l,T)),S),d:sub(mul(s,G),M),e:sub(mul(s,M),N),endpoint:u,sign:1};}
function trial(v,side,y){const {a,b,c,d,e,endpoint,sign}=gram(v,side),[x,z]=y;const corr=sub(scale(2,add(mul(d,x),mul(e,z))),sum(mul(a,sq(x)),scale(2,mul(b,mul(x,z))),mul(c,sq(z))));return div(add(v[0],scale(sign,corr)),endpoint);}
function optimal(v,side){const {a,b,c,d,e}=gram(v,side);const det=sub(mul(a,c),sq(b));check(lt(zero,a)&&lt(zero,det),"nonpositive Gram corner");const y=[div(sub(mul(c,d),mul(b,e)),det),div(sub(mul(a,e),mul(b,d)),det)];return {value:trial(v,side,y),y};}
function corners(box){return Array.from({length:32},(_,i)=>box.map((v,j)=>v[(i>>j)&1]));}
// Exact abstract Radau endpoint equalities, arbitrary positive weights.
for(const side of [0,1]){const ts=[side===0?l:u,q(1),q(4)],ws=[q(2),q(3),q(5)];const v=Array.from({length:5},(_,k)=>sum(...ts.map((t,i)=>{let x=ws[i];for(let j=0;j<k;j++)x=mul(x,sub(t,g));return x;})));const j=sum(...ts.map((t,i)=>div(ws[i],t)));check(same(optimal(v,side).value,j),"Radau equality");check(le(optimal(v,1).value,j)&&le(j,optimal(v,0).value),"abstract bound");}
check(lt(scale(8,add(mul(q(2n*10n**22n),q(1n,10n**70n)),mul(q(2n*10n**19n),q(32n,10n**35n)))),q(1n,10n**12n)),"shift");
const result={};for(const k of ["1","2"]){
const saved=fourth[k],rows={};for(const line of logs[Number(k)-1].trim().split("\n")){const [name,a,b,sep,den]=line.split(" ");check(!rows[name]&&sep==="/"&&den==="1e15","format");rows[name]=[q(a,10n**15n),q(b,10n**15n)];check(le(...rows[name]),"interval");}check(Object.keys(rows).length===10,"row count");
const [w,a,b,p]=["W2","sine4","bridge4","S_proxy"].map(n=>rows[n]);
check(overlap(p,[ma(zero,sub(sub(w[0],a[1]),b[1])),sub(sub(w[1],a[0]),b[0])]),"projection identity");
check(le(zero,rows.N_proxy[0])&&lt(rows.N_proxy[1],q("121/400"))&&overlap(rows.N_proxy,second[k].proxy_N.map(x=>q(x))),"N replay");
const eta=rows.projection_error[1];check(le(rows.projection_error[0],zero)&&le(zero,eta)&&lt(eta,q("1/10"))&&same(eta,q(third[k].projection_error)),"eta replay");
const delta=sum(mul(s,q(old[k].proxy_error_upper)),q("1/10000000000000"),eta),delta2=add(mul(s,delta),q("1/1000000000000")),B=sqrtup(p[1]),err=add(scale(2,mul(B,delta2)),sq(delta2));
const admissible=[q(support[k].support_lower_boundary[0]),q(support[k].support_upper_boundary[1])],raw=[ma(zero,sub(p[0],err)),add(p[1],err)];check(overlap(admissible,raw),"support");
const S=[ma(raw[0],admissible[0]),mi(raw[1],admissible[1])];
check(pair(S).join()===saved.true_S.join()&&render(err)===saved.transport&&render(delta2)===saved.second_application_error&&render(B)===saved.proxy_norm_bound,"native exact transport mismatch");
const box=[old[k].true_G,old[k].true_M,second[k].true_N,third[k].true_T].map(v=>v.map(x=>q(x))).concat([S]),K=old[k].K.map(x=>q(x));const cs=corners(box);
const fixedU=mi(ma(...cs.map(v=>trial(v,0,targets[k].trial_upper.map(x=>q(x))))),q(third[k].inverse_bounds[1])),fixedV=ma(mi(...cs.map(v=>trial(v,1,targets[k].trial_lower.map(x=>q(x))))),q(third[k].inverse_bounds[0]));
check(pair([fixedV,fixedU]).join()===saved.inverse_bounds.join(),"native fixed trial mismatch");
const mid=box.map(v=>div(add(...v),q(2))),yu=optimal(mid,0).y,yv=optimal(mid,1).y;
// U concave (inf affine), V convex (sup affine). Entire box Gram PD by convexity.
const lowerU=mi(...cs.map(v=>optimal(v,0).value)),upperV=ma(...cs.map(v=>optimal(v,1).value));
const upperU=ma(...cs.map(v=>trial(v,0,yu))),lowerV=mi(...cs.map(v=>trial(v,1,yv)));
check(le(lowerU,upperU)&&le(lowerV,upperV),"extrema enclosure");
const status=lt(upperU,K[0])?"PASS_WITNESS":lt(K[1],lowerV)?"FAIL_COMPARISON_WITNESS":lt(upperV,K[0])&&lt(K[1],lowerU)?"FIVE_SCALAR_MOMENTS_INSUFFICIENT":"UNDECIDED";
result[k]={status,native_transport_and_trials:"EXACT_MATCH",sharp_upper:pair([lowerU,upperU]),sharp_lower:pair([lowerV,upperV]),midpoint_upper_trial:pair(yu),midpoint_lower_trial:pair(yv),K:pair(K),true_S:pair(S),native_inverse:pair([fixedV,fixedU])};
}
return result;
}
const fs=require('fs');const args=process.argv.slice(2);const out=audit(...args.slice(0,6).map(p=>JSON.parse(fs.readFileSync(p,'utf8'))),args.slice(6,8).map(p=>fs.readFileSync(p,'utf8')));if(args[8])fs.writeFileSync(args[8],JSON.stringify(out,null,2)+'\n');else console.log(JSON.stringify(out,null,2));
