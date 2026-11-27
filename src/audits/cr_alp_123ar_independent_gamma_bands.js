// CR-ALP-123AR: independent exact Gamma band audit.
function audit(b){
  function gcd(a,b){a=a<0n?-a:a;while(b){[a,b]=[b,a%b]}return a||1n}
  function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){let [n,d]=a.split("/");return q(n,d)}a=BigInt(a);b=BigInt(b);if(b<0n){a=-a;b=-b}if(b===0n)throw Error("zero");let d=gcd(a,b);return [a/d,b/d]}
  const add=(a,b)=>q(a[0]*b[1]+b[0]*a[1],a[1]*b[1]),neg=a=>q(-a[0],a[1]),sub=(a,b)=>add(a,neg(b)),mul=(a,b)=>q(a[0]*b[0],a[1]*b[1]),div=(a,b)=>q(a[0]*b[1],a[1]*b[0]),pow=(a,n)=>q(a[0]**BigInt(n),a[1]**BigInt(n)),cmp=(a,b)=>a[0]*b[1]<b[0]*a[1]?-1:a[0]*b[1]>b[0]*a[1]?1:0,min=(a,b)=>cmp(a,b)<0?a:b,fl=a=>a[0]/a[1],st=a=>a.join("/");
  function check(t,msg){if(!t)throw Error(msg)}
  function sum(xs){return xs.reduce(add,q(0))}
  function isqrt(n){if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));for(;;){let y=(x+n/x)/2n;if(y>=x)return x;x=y}}
  function sqr(a,upper){const g=10n**18n;return q(isqrt(a[0]*g*g/a[1])+(upper?1n:0n),g)}
  function loglo(x){let k=0;while(cmp(x,q(2))>=0){x=div(x,q(2));k++}function ls(z){let t=div(sub(z,q(1)),add(z,q(1)));return mul(q(2),sum(Array.from({length:32},(_,i)=>div(pow(t,2*i+1),q(2*i+1)))))}return add(ls(x),mul(q(k),ls(q(2))))}
  function expLo(x,n){let t=q(1),s=t;for(let k=1;k<=n;k++){t=div(mul(t,x),q(k));s=add(s,t)}return s}
  function atan(x,n){return sum(Array.from({length:n},(_,i)=>mul(q(i%2?-1:1,2*i+1),pow(x,2*i+1))))}
  let piLo=sub(mul(q(16),atan(q(1,5),10)),mul(q(4),atan(q(1,239),9))),piHi=sub(mul(q(16),atan(q(1,5),9)),mul(q(4),atan(q(1,239),10)));
  check(cmp(piLo,q(314159,100000))>0&&cmp(piHi,q(22,7))<0,"pi");
  check(cmp(expLo(q(46,25),16),q(44,7))>0,"log 2pi");
  check(cmp(expLo(q(11,10),6),q(2))>0&&cmp(expLo(q(5),10),q(100))>0,"exp");
  let A=q(11,10),V=q(120),h=q(57,20),r0=q(109,10),r2=q(7,5);
  check(cmp(pow(r0,2),V)<0&&cmp(pow(r2,2),q(2))<0&&cmp(div(piLo,A),h)>0&&cmp(div(mul(q(86),piLo),mul(q(2),A)),V)>0,"geometry");
  let sm=mul(mul(div(q(2),A),q(24)),mul(div(add(div(q(1),V),div(q(2),h)),pow(r0,2)),pow(div(q(32),V),4)));
  let C=mul(mul(div(q(2),A),q(101,100)),add(div(A,r2),q(1,32))),I=div(q(22,7),mul(q(2),r2)),total=add(sm,mul(C,I));
  check(cmp(total,q(b.boundary_constant))===0,"total");
  let mD=null,mSaved=null,mTail=null,minSymbolReserve=null,rows=[];
  for(let [i,z] of b.bands.entries()){
    check(z.start===(i===0?86:b.bands[i-1].end+1),"coverage");
    const k=fl(div(mul(q(z.start),q(314159,100000)),mul(q(2),A)));
    check(k===BigInt(z.k)&&k>=120n,"k");
    let g=sub(sub(loglo(q(k)),q(46,25)),q(22,7n*k));
    let gOld=sub(sub(q(fl(mul(loglo(q(k)),q(10n**9n))),10n**9n),q(37,20)),q(2,k));
    check(cmp(g,gOld)>0,"direct g");
    let reserve=sub(g,gOld);minSymbolReserve=minSymbolReserve===null?reserve:min(minSymbolReserve,reserve);
    let p=total;
    if(z.end!==null){
      let v=div(mul(q(z.end),q(22,7)),mul(q(2),A));let ratio=sub(q(1),div(pow(V,2),pow(v,2)));
      check(z.end>=z.start&&cmp(ratio,q(0))>=0&&cmp(ratio,q(1))<0,"band geometry");
      if(z.split_beta!==null){
        let B=q(z.split_beta);check(cmp(B,q(64))>=0,"split");let B2=pow(B,2),alpha=pow(div(B2,add(pow(V,2),B2)),2),bm=sub(B,q(2)),tail=min(I,div(mul(v,sqr(v,true)),mul(q(3),mul(bm,sqr(bm,false)))));
        p=add(sm,mul(C,add(mul(alpha,I),mul(sub(q(1),alpha),tail))))
      }
      p=mul(pow(ratio,2),min(total,p))
    }
    check(cmp(p,q(z.penalty_upper))<=0,`penalty ${z.start}`);
    let D=sub(g,p),saved=q(z.gamma_lower);
    check(cmp(D,saved)>=0,`gamma ${z.start}`);
    mD=mD===null?D:min(mD,D);mSaved=mSaved===null?saved:min(mSaved,saved);if(z.start>=2049)mTail=mTail===null?saved:min(mTail,saved);
    rows.push({start:z.start,end:z.end,direct_lower:st(D),saved_lower:st(saved),penalty_upper:st(p),classification:"PASS"});
  }
  check(rows.length===118&&b.bands.at(-1).end===null&&cmp(mSaved,q(72,25))>0&&cmp(mTail,q("4447938456/1000000000"))>=0,"terminal");
  return {status:"CONDITIONAL_GAMMA_COMPARATOR_AUDIT_PASS",small:st(sm),total:st(total),min_saved_gamma:st(mSaved),min_saved_gamma_1e9:String(fl(mul(mSaved,q(10n**9n)))),min_direct_gamma_1e9:String(fl(mul(mD,q(10n**9n)))),min_symbol_reserve_1e9:String(fl(mul(minSymbolReserve,q(10n**9n)))),min_tail_gamma_1e9:String(fl(mul(mTail,q(10n**9n)))),bands:rows,scope:"Original resolvent identity and prime norm premise retained; independent rational band check and direct digamma symbol lower."}
}
const fs=require('fs'),crypto=require('crypto');
const raw=fs.readFileSync(process.argv[2]);
if(crypto.createHash('sha1').update(Buffer.from('blob '+raw.length+'\0')).update(raw).digest('hex')!=='1a7fe8629daa8b8e6cec185a79b6746d6b28975d')throw Error('input blob');
const result=audit(JSON.parse(raw));
const output=JSON.stringify(result,null,2)+'\n';
if(process.argv[3])fs.writeFileSync(process.argv[3],output);else process.stdout.write(output);
