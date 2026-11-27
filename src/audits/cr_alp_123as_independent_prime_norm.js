// CR-ALP-123AS: independently recompute exact prime translation norm.
function primeAudit(){
  function gcd(a,b){while(b){[a,b]=[b,a%b]}return a}
  function q(a,b=1n){a=BigInt(a);b=BigInt(b);if(b<=0n||a<0n)throw Error("nonnegative rational required");let g=gcd(a,b);return [a/g,b/g]}
  const add=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]),mul=(x,y)=>q(x[0]*y[0],x[1]*y[1]),div=(x,y)=>q(x[0]*y[1],x[1]*y[0]),sub=(x,y)=>q(x[0]*y[1]-y[0]*x[1],x[1]*y[1]),cmp=(x,y)=>x[0]*y[1]<y[0]*x[1]?-1:x[0]*y[1]>y[0]*x[1]?1:0,key=x=>x.join("/"),ceil=x=>(x[0]+x[1]-1n)/x[1];
  function check(x,s){if(!x)throw Error(s)}
  function sqrtfloor(n){let lo=0n,hi=1n;while(hi*hi<=n)hi*=2n;while(hi-lo>1n){let m=(lo+hi)/2n;if(m*m<=n)lo=m;else hi=m}return lo}
  function log12(x){check(cmp(x,q(1))>=0&&cmp(x,q(2))<=0,"log range");let u=div(sub(x,q(1)),add(x,q(1))),u2=mul(u,u),term=u,s=q(0);for(let i=0;i<48;i++){s=add(s,div(mul(q(2),term),q(2*i+1)));term=mul(term,u2)}return [s,add(s,div(mul(q(2),term),mul(q(97),sub(q(1),u2))))]}
  function log(n){let p=1n,k=0n;while(2n*p<=n){p*=2n;k++}let x=log12(q(n,p)),l=log12(q(2));return x.map((v,i)=>add(v,mul(q(k),l[i])))}
  check(cmp(log(9n)[1],q(11,5))<0&&cmp(q(11,5),log(10n)[0])<0,"collar");
  const N=[2n,3n,4n,5n,7n,8n,9n],P=[2n,3n,2n,5n,7n,2n,3n],G=10n**9n,S=10n**18n;
  const W=N.map((n,i)=>{let r=q(sqrtfloor(n*S*S),S);check(cmp(mul(r,r),q(n))<=0,"root");return ceil(mul(q(G),div(log(P[i])[1],r)))});
  const old=[490130n,634285n,346574n,719763n,735485n,245065n,366205n];
  check(W.every((w,i)=>w<=old[i]*1000n),"weight refinement");
  const inside=x=>cmp(q(1),x)<0&&cmp(x,q(10))<0;
  let bp=[q(1),q(10)],v=[1n],steps=[];
  for(let k=1;k<=10;k++){
    let map=new Map(bp.map(x=>[key(x),x]));
    for(let x of bp)for(let n of N)for(let y of [mul(x,q(n)),div(x,q(n))])if(inside(y))map.set(key(y),y);
    let next=[...map.values()].sort(cmp),values=[];
    for(let i=0;i<next.length-1;i++){
      const l=next[i],h=next[i+1];let val=0n;
      for(let j=0;j<N.length;j++)for(let s of [q(N[j]),q(1,N[j])]){
        const L=mul(l,s),H=mul(h,s);
        if(cmp(H,q(1))<=0||cmp(L,q(10))>=0)continue;
        check(cmp(L,q(1))>=0&&cmp(H,q(10))<=0,"support crossing");
        let a=0,b=bp.length-1;
        while(b-a>1){let m=(a+b)>>1;if(cmp(bp[m],L)<=0)a=m;else b=m}
        check(cmp(bp[a],L)<=0&&cmp(H,bp[a+1])<=0,"old cell crossing");
        val+=W[j]*v[a];
      }
      values.push(val);
    }
    bp=next;v=values;let m=v.reduce((a,b)=>a>b?a:b,0n);
    steps.push({power:k,cells:v.length,max_row_mass_numerator:String(m),denominator:String(G**BigInt(k))});
  }
  check(v.length===2189,"cell count");let m=BigInt(steps[9].max_row_mass_numerator),d=G**10n;
  check(m*8n**10n<d*21n**10n,"norm");
  let up=2600n;while(m*1000n**10n>=d*up**10n)up++;
  return {status:"INDEPENDENT_PRIME_NORM_PASS",weights:Object.fromEntries(N.map((n,i)=>[String(n),String(W[i])])),weight_grid:String(G),steps,verified_norm_upper:"21/8",refined_millesimal_upper:String(up)+"/1000",scope:"Original seven paired translations; enlarged interval (1,10) only for norm domination. No new source or support claim."}
}
const fs=require('fs'); const output=JSON.stringify(primeAudit(),null,2)+'\n';
if(process.argv[2])fs.writeFileSync(process.argv[2],output);else process.stdout.write(output);
