// Exact arithmetic for CR-ALP-066 full-frequency Gamma comparison.
function audit(){
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){const c=a%b;a=b;b=c;}return a||1n;}
function q(n,d=1n){n=BigInt(n);d=BigInt(d);if(d<0n){n=-n;d=-d;}const g=gcd(n,d);return[n/g,d/g];}
const a=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]),s=(x,y)=>a(x,q(-y[0],y[1])),m=(x,y)=>q(x[0]*y[0],x[1]*y[1]),d=(x,y)=>q(x[0]*y[1],x[1]*y[0]),lt=(x,y)=>x[0]*y[1]<y[0]*x[1],pow=(x,n)=>q(x[0]**BigInt(n),x[1]**BigInt(n)),sum=x=>x.reduce(a,q(0)),str=x=>x[0]+"/"+x[1],floor=x=>x[0]/x[1];
function atan(x,n){return sum(Array.from({length:n},(_,i)=>m(q(i%2?-1:1,2*i+1),pow(x,2*i+1))));}
const piL=s(m(q(16),atan(q(1,5),8)),m(q(4),atan(q(1,239),7))),piU=s(m(q(16),atan(q(1,5),7)),m(q(4),atan(q(1,239),8)));
const A=q(11,10),v0=q(120),rt0=q(109,10),spacing=q(57,20),r2=q(7,5);
if(!lt(q(314159,100000),piL)||!lt(piU,q(22,7))||!lt(v0,d(m(q(86),piL),m(q(2),A)))||!lt(spacing,d(piL,A))||!lt(pow(rt0,2),v0)||!lt(pow(r2,2),q(2)))throw Error("geometry");
let fact=1n;const exp5=sum(Array.from({length:7},(_,k)=>{if(k)fact*=BigInt(k);return q(5n**BigInt(k),fact);}));
if(!lt(q(100),exp5)||!lt(q(2),a(q(1),A)))throw Error("exponential");
const small=m(m(d(q(2),A),q(24)),m(d(a(d(q(1),v0),d(q(2),spacing)),rt0),d(q(1),rt0)));
const centered=m(small,pow(d(q(32),v0),4));
const tail=m(m(d(q(2),A),q(101,100)),m(a(d(A,r2),q(1,32)),d(q(22,7),m(q(2),r2))));
const total=a(centered,tail);if(!lt(total,q(7,4)))throw Error("penalty");
function loglo(x){const v=d(s(x,q(1)),a(x,q(1)));return m(q(2),sum(Array.from({length:24},(_,i)=>d(pow(v,2*i+1),q(2*i+1)))));}
const l2=loglo(q(2));
function logint(k){let h=BigInt(k),b=0;while((1n<<BigInt(b+1))<=h)b++;return q(floor(m(a(m(q(b),l2),loglo(q(h,1n<<BigInt(b)))),q(10**9))),10**9);}
const gamma=q(72,25),sP=q(21,8),theta=q(1,10);
const bands=[86,124,201,321,481,601,701,851,961,1251,2049,4097,8193].map(j=>{const k=floor(d(m(q(j),q(314159,100000)),m(q(2),A)));const g=s(s(logint(k),q(37,20)),q(2,k));const D=s(g,q(7,4)),gain=s(D,gamma),blend=s(a(m(s(q(1),theta),gamma),m(theta,D)),sP);return {j,k:String(k),gamma_diagonal_lower:str(D),gain_over_gamma:str(gain),blended_weil_diagonal_lower:str(blend)};});
if(!lt(q(0),q(...bands[0].blended_weil_diagonal_lower.split("/")))||!lt(q(0),q(...bands.find(x=>x.j===701).gain_over_gamma.split("/"))))throw Error("band");
return {small_uncentered:str(small),small_pole_centered:str(centered),large_beta:str(tail),total_penalty:str(total),certified_penalty:"7/4",theta:"1/10",bands};
}
const result=audit();const fs=require('fs');if(process.argv[2])fs.writeFileSync(process.argv[2],JSON.stringify(result,null,2)+'\n');else console.log(JSON.stringify(result,null,2));
