// Exact finite diagnostics for CROSS-PW-004. The polynomial model is not the infinite-convolution detector.
function replay(){
 const assert=(p,m)=>{if(!p)throw Error(m)};
 const gcd=(a,b)=>b===0n?(a<0n?-a:a):gcd(b,a%b);
 const add=([a,b],[c,d])=>{const p=a*d+b*c,q=b*d,g=gcd(p,q);return [p/g,q/g]};
 const deriv=a=>a.slice(1).map((v,i)=>v*BigInt(i+1));
 const powD=(a,k)=>{for(let i=0;i<k;i++)a=deriv(a);return a};
 let h=[1n];
 for(let k=0;k<20;k++){const b=Array(h.length+2).fill(0n);h.forEach((v,i)=>{b[i]+=v;b[i+2]-=v});h=b}
 const d16=powD(h,16),d18=powD(h,18);
 const f=d16.map((v,i)=>4n*(d18[i]??0n)-v);
 for(let k=0;k<16;k++){
  let s=[0n,1n];f.forEach((v,i)=>{if((i+k)%2===0)s=add(s,[2n*v,BigInt(i+k+1)])});
  assert(s[0]===0n,"moment "+k);
 }
 assert(f.some(v=>v!==0n)&&f.every((v,i)=>i%2===0||v===0n),"nonzero even polynomial");
 for(let n=1;n<=128;n++){let sum=0n;for(let j=1;j<=n;j++)sum+=1n<<BigInt(n-j);assert(sum===(1n<<BigInt(n))-1n,"radius sum");}
 // numerator of 4*w^18-w^16 at rational w=p/q, common q^18
 const symbol=(p,q)=>p**16n*(4n*p*p-q*q);
 assert(symbol(1n,2n)===0n&&symbol(-1n,2n)===0n,"pole factors");
 return {radius_partial_sums:128,polynomial_model_moments:16,pole_symbol_checks:2,even_nonzero_model:true,status:"PASS"}
}
console.log(JSON.stringify(replay()));
