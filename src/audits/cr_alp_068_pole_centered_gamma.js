// Exact arithmetic for CR-ALP-068 all-resolvent pole centering and split-tail bound.
function audit(prior){
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){let c=a%b;a=b;b=c;}return a||1n;}
function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){let [n,d]=a.split("/");return q(BigInt(n),BigInt(d));}a=BigInt(a);b=BigInt(b);if(b<=0n)throw Error("den");const c=gcd(a,b);return[a/c,b/c];}
const a=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]),s=(x,y)=>a(x,q(-y[0],y[1])),m=(x,y)=>q(x[0]*y[0],x[1]*y[1]),d=(x,y)=>q(x[0]*y[1],x[1]*y[0]),cmp=(x,y)=>x[0]*y[1]<y[0]*x[1]?-1:x[0]*y[1]>y[0]*x[1]?1:0,pow=(x,n)=>q(x[0]**BigInt(n),x[1]**BigInt(n)),min=(x,y)=>cmp(x,y)<0?x:y,str=x=>x[0]+"/"+x[1],fl=x=>x[0]/x[1];
function isqrt(n){if(n<0n)throw Error("sqrt");if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));while(true){const y=(x+n/x)/2n;if(y>=x)return x;x=y;}}
function sqrt(x,upper){const k=10n**12n,v=isqrt(x[0]*k*k/x[1]);return q(v+(upper?1n:0n),k);}
function loglo(x){const u=d(s(x,q(1)),a(x,q(1)));let z=q(0);for(let i=0;i<24;i++)z=a(z,d(pow(u,2*i+1),q(2*i+1)));return m(q(2),z);}
const l2=loglo(q(2));function logint(k){let b=0;while((1n<<BigInt(b+1))<=k)b++;return q(fl(m(a(m(q(b),l2),loglo(q(k,1n<<BigInt(b)))),q(10**9))),10**9);}
const aa=q(11,10),v0=q(120),piL=q(314159,100000),piU=q(22,7),sqrt2L=q(7,5),S=q(21,8),gamma=q(72,25);
if(cmp(d(m(q(86),piL),m(q(2),aa)),v0)<=0||cmp(pow(sqrt2L,2),q(2))>=0)throw Error("geometry");
const C=m(m(d(q(2),aa),q(101,100)),a(d(aa,sqrt2L),q(1,32))),Imax=d(piU,m(q(2),sqrt2L));
const small=q(prior.small_pole_centered),oldTotal=q(prior.total_penalty);
if(cmp(a(small,m(C,Imax)),oldTotal)!==0||cmp(oldTotal,q(17,10))>=0)throw Error("prior constants");
const starts=[];for(let j=86;j<512;j+=4)starts.push(j);starts.push(512,601,701,851,961,1251,1537,2049,4097,8193,16385);
const out=[];let global=null;
for(let i=0;i<starts.length;i++){const j=starts[i],end=starts[i+1]?starts[i+1]-1:null,k=fl(d(m(q(j),piL),m(q(2),aa))),g=s(s(logint(k),q(37,20)),q(2,k));let pen=oldTotal,chosen=null;
if(end!==null){const vmax=d(m(q(end),piU),m(q(2),aa)),ratio=s(q(1),d(pow(v0,2),pow(vmax,2))),v32=m(vmax,sqrt(vmax,true));
if(cmp(ratio,q(0))<0||cmp(ratio,q(1))>=0)throw Error("ratio");
let best=oldTotal;
for(const B of [64,96,128,160,192,224,256,320,384,512,640,768,1024,1536,2048,4096]){
 const b=q(B),alpha=pow(d(pow(b,2),a(pow(v0,2),pow(b,2))),2);
 // First tail beta>B has preceding interval [beta-2,beta], all above B-2.
 const bm=q(B-2),tail=min(Imax,d(v32,m(q(3),m(bm,sqrt(bm,false)))));
 const candidate=a(small,m(C,a(m(alpha,Imax),m(s(q(1),alpha),tail))));
 if(cmp(candidate,best)<0){best=candidate;chosen=B;}
}
pen=m(pow(ratio,2),best);
}
const D=s(g,pen),coerc=s(D,S);if(global===null||cmp(coerc,global)<0)global=coerc;
out.push({start:j,end,k:String(k),split_beta:chosen,penalty_upper:str(pen),gamma_lower:str(D),weil_lower:str(coerc),gain_over_scalar_weil_floor:str(s(D,gamma))});
}
if(cmp(global,q(0))<=0)throw Error("no global coercivity");
return {reference_frequency:"120",boundary_constant:str(oldTotal),global_weil_floor:str(global),bands:out};
}
const fs=require('fs');const result=audit(JSON.parse(fs.readFileSync(process.argv[2],'utf8')));if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(result,null,2)+'\n');else console.log(JSON.stringify(result,null,2));
