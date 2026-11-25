// CR-ALP-062 candidate: exact BigInt rational matrix regression.
function audit(){
  const abs=x=>x<0n?-x:x;
  function gcd(a,b){a=abs(a);b=abs(b);while(b){let t=a%b;a=b;b=t;}return a||1n;}
  function q(a,b=1n){a=BigInt(a);b=BigInt(b);if(!b)throw Error("zero denominator");if(b<0n){a=-a;b=-b;}let c=gcd(a,b);return [a/c,b/c];}
  const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]), neg=a=>q(-a[0],a[1]),sub=(a,b)=>add(a,neg(b)),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]),eq=(a,b)=>a[0]===b[0]&&a[1]===b[1];
  const tr=A=>A[0].map((_,i)=>A.map(r=>r[i])),plus=(A,B)=>A.map((r,i)=>r.map((x,j)=>add(x,B[i][j]))),minus=(A,B)=>A.map((r,i)=>r.map((x,j)=>sub(x,B[i][j]))),scale=(k,A)=>A.map(r=>r.map(x=>mul(k,x))),prod=(A,B)=>A.map(r=>tr(B).map(c=>r.reduce((s,x,i)=>add(s,mul(x,c[i])),q(0)))),diag=xs=>xs.map((x,i)=>xs.map((_,j)=>i===j?x:q(0))),gram=(R,D)=>prod(tr(R),prod(D,R)),eye=n=>diag(Array.from({length:n},()=>q(1))),mat=A=>A.map(r=>r.map(x=>q(x))),same=(A,B)=>A.every((r,i)=>r.every((x,j)=>eq(x,B[i][j])));
  const psd=A=>eq(A[0][1],A[1][0])&&A[0][0][0]>=0n&&A[1][1][0]>=0n&&sub(mul(A[0][0],A[1][1]),mul(A[0][1],A[1][0]))[0]>=0n;
  const check=(b,msg)=>{if(!b)throw Error(msg);};
  const ell=q(51,200),u=q(1101,200),g=q(72,25),s=q(21,8);
  function test(ts,R,Y,A,B,which){
    const p=ts.map(t=>sub(t,g)),G=gram(R,eye(ts.length)),M=gram(R,diag(p)),N=gram(R,diag(p.map(x=>mul(x,x)))),T=gram(R,diag(p.map(x=>mul(mul(x,x),x)))),J=gram(R,diag(ts.map(t=>div(q(1),t))));
    const H0=minus(scale(mul(s,s),G),N),H1=minus(plus(scale(g,H0),scale(mul(s,s),M)),T),HY=prod(H0,Y);
    const U=scale(div(q(1),mul(ell,u)),plus(minus(minus(scale(g,G),M),plus(HY,tr(HY))),gram(Y,H1)));
    const D=plus(prod(G,A),prod(M,B)),X=plus(prod(R,A),prod(prod(diag(p),R),B)),V=minus(plus(D,tr(D)),gram(X,diag(ts)));
    check(psd(minus(U,J))&&psd(minus(J,V)),"Loewner comparison");
    if(which==="lower")check(same(V,J),"Gaussian lower equality");
    if(which==="upper")check(same(U,J),"Lobatto upper equality");
    if(which==="both")check(same(U,J)&&same(V,J),"endpoint equality");
  }
  const R=mat([[1,2],[3,-1],[2,1],[1,-2]]),Y=mat([[1,0],[1,2]]),A=mat([[1,2],[0,1]]),B=scale(q(-1,7),mat([[1,0],[0,2]]));
  test([q(1,2),q(2),q(4),q(5)],R,Y,A,B,"none");
  test([q(1),q(4)],mat([[1,1],[0,1]]),Y,scale(div(sub(q(5),g),q(4)),eye(2)),scale(q(-1,4),eye(2)),"lower");
  test([ell,q(2),u],mat([[1,2],[3,-1],[2,1]]),scale(q(1,2),eye(2)),A,B,"upper");
  test([ell,u],mat([[1,1],[0,1]]),Y,scale(div(g,mul(ell,u)),eye(2)),scale(div(q(-1),mul(ell,u)),eye(2)),"both");
  return "PASS: noncommuting rational trials, Gaussian lower equality, Lobatto upper equality, singular endpoint case; no actual Weil moment evaluated.";
}
console.log(audit());
