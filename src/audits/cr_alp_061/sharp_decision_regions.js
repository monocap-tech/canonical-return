// Exact rational third-moment decision regions from frozen 058/059 JSON.
function audit(old,second){
const abs=x=>x<0n?-x:x;
function gcd(a,b){a=abs(a);b=abs(b);while(b){const t=a%b;a=b;b=t;}return a||1n;}
function rat(a,b=1n){if(typeof a==="string"&&a.includes("/")){const v=a.split("/");return rat(BigInt(v[0]),BigInt(v[1]));}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}const c=gcd(a,b);return [a/c,b/c];}
const add=(a,b)=>rat(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>rat(-a[0],a[1]),mul=(a,b)=>rat(a[0]*b[0],a[1]*b[1]),inv=a=>rat(a[1],a[0]),cmp=(a,b)=>a[0]*b[1]<b[0]*a[1]?-1:a[0]*b[1]>b[0]*a[1]?1:0;
const I=(a,b=a)=>[rat(a),rat(b)],pt=a=>[a,a],plus=(a,b)=>[add(a[0],b[0]),add(a[1],b[1])],minus=(a,b)=>plus(a,[neg(b[1]),neg(b[0])]),times=(a,b)=>{const v=a.flatMap(x=>b.map(y=>mul(x,y))).sort(cmp);return [v[0],v[3]];},divide=(a,b)=>{if(b[0][0]<=0n)throw Error("nonpositive denominator");return times(a,[inv(b[1]),inv(b[0])]);},sq=a=>{const p=mul(a[0],a[0]),q=mul(a[1],a[1]);return [a[0][0]<=0n&&a[1][0]>=0n?rat(0):cmp(p,q)<0?p:q,cmp(p,q)>0?p:q];};
const s=I("21/8"),g=I("72/25"),ell=I("51/200"),u=I("1101/200"),su=times(ell,u),ss=sq(s);
const res={};for(const k of ["1","2"]){
const G=I(...old[k].true_G),M=I(...old[k].true_M),N=I(...second[k].true_N),K=I(...old[k].K);
const h=minus(times(ss,G),N),f=minus(minus(times(g,G),M),times(su,K));
const m=plus(times(g,G),M),d=plus(times(g,M),N),j0=divide(sq(G),m),gap=minus(K,j0);
const accept=minus(plus(times(g,h),times(ss,M)),divide(sq(h),f));
const reject=plus(minus(divide(sq(d),m),times(g,N)),divide(sq(minus(M,divide(times(d,G),m))),gap));
const tmin=minus(divide(sq(plus(times(s,M),N)),plus(times(s,G),M)),times(s,N));
const tmax=minus(times(s,N),divide(sq(minus(times(s,M),N)),minus(times(s,G),M)));
const render=a=>a.map(x=>x[0]+"/"+x[1]); const decimal=a=>a.map(x=>Number(x[0])/Number(x[1]));
if(cmp(reject[1],accept[0])>=0)throw Error("overlapping decision thresholds");
res[k]={accept_boundary:render(accept),reject_boundary:render(reject),support_lower_boundary:render(tmin),support_upper_boundary:render(tmax)};
console.log({sector:k,accept:decimal(accept),reject:decimal(reject),support_min:decimal(tmin),support_max:decimal(tmax)});
}
return res;
}
const fs=require('fs'); const output=audit(...process.argv.slice(2,4).map(p=>JSON.parse(fs.readFileSync(p,'utf8')))); if(process.argv[4])fs.writeFileSync(process.argv[4],JSON.stringify(output,null,2)+'\n');
