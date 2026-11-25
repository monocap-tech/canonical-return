// Frozen linear weighted inverse trials, exact endpoint checks.
function audit(old,second,third){
 const abs=x=>x<0n?-x:x;function gcd(a,b){a=abs(a);b=abs(b);while(b){let c=a%b;a=b;b=c;}return a||1n;}
 function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){let [x,y]=a.split("/");return q(BigInt(x),BigInt(y));}a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero");if(b<0n){a=-a;b=-b;}let c=gcd(a,b);return [a/c,b/c];}
 const a=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]),neg=x=>q(-x[0],x[1]),sub=(x,y)=>a(x,neg(y)),m=(x,y)=>q(x[0]*y[0],x[1]*y[1]),d=(x,y)=>q(x[0]*y[1],x[1]*y[0]),sq=x=>m(x,x),cmp=(x,y)=>x[0]*y[1]<y[0]*x[1]?-1:x[0]*y[1]>y[0]*x[1]?1:0,dec=x=>Number(x[0])/Number(x[1]),mid=x=>(dec(x[0])+dec(x[1]))/2;
 const g=q("72/25"),s=q("21/8"),ell=q("51/200"),u=q("1101/200"),gs=m(g,s);
 const out={};for(const k of ["1","2"]){
 const X=[old[k].true_G,old[k].true_M,second[k].true_N,third[k].true_T,old[k].K].map(v=>v.map(x=>q(x)));
 const [G,M,N,T,K]=X.map(mid),gg=dec(g),ss=dec(s),ll=dec(ell),uu=dec(u);
 const A=gg*ss*G+uu*M+N,B=gg*ss*M+uu*N+T,C0=gg*ss*N+uu*T,b0=ss*G+M,b1=ss*M+N;
 const D=gg*ss*G-ll*M-N,E=gg*ss*M-ll*N-T,F0=gg*ss*N-ll*T,d0=ss*G-M,d1=ss*M-N;
 const U=G/ll-b0*b0/(ll*A),V=G/uu+d0*d0/(uu*D),ar=b1-B*b0/A,dr=d1-E*d0/D;
 const beta=ar/(ar*ar/(ll*(U-K))),alpha=(b0-B*beta)/A,delta=dr/(dr*dr/(uu*(K-V))),zeta=(d0-E*delta)/D;
 const trial=[alpha,beta,zeta,delta].map(x=>q(Math.round(10000*x),10000));
 const [y,b,z,c]=trial;const up=[],lo=[];
 for(let mask=0;mask<32;mask++){
 const [G,M,N,T,K]=X.map((v,j)=>v[(mask>>j)&1]);
 const A=a(a(m(gs,G),m(u,M)),N),B=a(a(m(gs,M),m(u,N)),T),C0=a(m(gs,N),m(u,T)),b0=a(m(s,G),M),b1=a(m(s,M),N);
 const D=sub(sub(m(gs,G),m(ell,M)),N),E=sub(sub(m(gs,M),m(ell,N)),T),F0=sub(m(gs,N),m(ell,T)),d0=sub(m(s,G),M),d1=sub(m(s,M),N);
 const acc=d(sub(sub(sub(a(sub(m(ell,K),G),m(q(2),a(m(b0,y),m(b1,b)))),m(A,sq(y))),m(q(2),m(B,m(y,b)))),m(C0,sq(b))),sq(b));
 const rej=d(a(a(a(sub(sub(m(u,K),G),m(q(2),a(m(d0,z),m(d1,c)))),m(D,sq(z))),m(q(2),m(E,m(z,c)))),m(F0,sq(c))),sq(c));
 up.push(acc);lo.push(rej);
 }
 up.sort(cmp);lo.sort(cmp);
 const render=x=>x[0]+"/"+x[1];
 out[k]={trial_upper:[render(y),render(b)],trial_lower:[render(z),render(c)],sufficient_acceptance_strict_below:render(up[0]),sufficient_rejection_strict_above:render(lo.at(-1))};
 }
 return out;
}
const fs=require('fs');const r=audit(...process.argv.slice(2,5).map(p=>JSON.parse(fs.readFileSync(p,'utf8'))));if(process.argv[5])fs.writeFileSync(process.argv[5],JSON.stringify(r,null,2)+'\n');else console.log(JSON.stringify(r,null,2));
