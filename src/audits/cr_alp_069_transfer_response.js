// Exact outward fixed-point interval response transfer; no quadrature.
function audit(trial,diag,log,newJSON,oldJSON,oldBandsJSON){
 const S=10n**30n,zero=[0n,0n],point=x=>[x,x],one=point(S),min=x=>x.reduce((a,b)=>a<b?a:b),max=x=>x.reduce((a,b)=>a>b?a:b);
 function fl(n,d){if(d<=0n)throw Error("bad denominator");return n>=0n?n/d:-((-n+d-1n)/d);}
 const ce=(n,d)=>-fl(-n,d);
 function rat(s){const t=String(s);if(t.includes("/")){const [n,d]=t.split("/").map(BigInt);return[fl(n*S,d),ce(n*S,d)];}if(!/^-?\d+(\.\d+)?$/.test(t))throw Error("number "+t);const neg=t.startsWith("-"),[a,b=""]=t.replace("-","").split(".");if(b.length>30)throw Error("decimal");const n=(BigInt(a)*S+BigInt((b+"0".repeat(30)).slice(0,30)))*(neg?-1n:1n);return point(n);}
 const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],sub=(a,b)=>[a[0]-b[1],a[1]-b[0]],mul=(a,b)=>{const p=a.flatMap(x=>b.map(y=>x*y));return[fl(min(p),S),ce(max(p),S)]},div=(a,b)=>{if(b[0]<=0n)throw Error("divide nonpositive");const p=a.flatMap(x=>b.map(y=>fl(x*S,y))),q=a.flatMap(x=>b.map(y=>ce(x*S,y)));return[min(p),max(q)]},sq=a=>{const z=a.map(x=>x*x);return[a[0]<=0n&&a[1]>=0n?0n:fl(min(z),S),ce(max(z),S)]},sum=x=>x.reduce(add,zero),str=x=>{const neg=x<0n,n=neg?-x:x;return(neg?"-":"")+(n/S)+"."+(n%S).toString().padStart(30,"0")},pair=x=>x.map(str),overlap=(a,b)=>a[0]<=b[1]&&b[0]<=a[1];
 function isqrt(n){if(n<0n)throw Error("sqrt");if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));while(true){const y=(x+n/x)/2n;if(y>=x)return x;x=y;}}
 const root=x=>{if(x[0]<0n)throw Error("negative root");return[isqrt(x[0]*S),isqrt(x[1]*S)+1n]},up=a=>point(a[1]),down=a=>point(a[0]);
 const bands=JSON.parse(newJSON).bands,oldBands=JSON.parse(oldBandsJSON).bands,old=JSON.parse(oldJSON)["2"],jids=Array.from({length:982},(_,i)=>86+2*i),all=Array.from({length:1024},(_,i)=>2+2*i);
 const bank=t=>{const a=t.trim().split("\n").map(x=>x.split(" "));if(a.length!==982||a.some((x,i)=>Number(x[0])!==jids[i]||x[1]!==x[2]))throw Error("bank");return a.map(x=>rat(x[1]));};
 const x0=bank(trial),d0=bank(diag),r={};for(const line of log.trim().split("\n")){const [k,l,h,slash,sc]=line.split(" ");if(slash!=="/"||sc!=="1e15"||k in r)throw Error("row");r[k]=[BigInt(l)*10n**15n,BigInt(h)*10n**15n];if(r[k][0]>r[k][1])throw Error("interval");}
 if(Object.keys(r).length!==2054)throw Error("row count");
 const gamma=rat("72/25"),theta=rat("1/10"),prime=rat("21/8"),grid=10n**21n;
 const lower9=a=>point(fl(a[0],grid)*grid);
 const newD=b=>lower9(rat(b.gamma_lower)),oldD=b=>lower9(add(mul(sub(one,theta),gamma),mul(theta,rat(b.gamma_diagonal_lower))));
 if(bands[0].start!==86||bands.at(-1).end!==null||bands.some((b,i)=>i&&b.start!==bands[i-1].end+1))throw Error("coverage");
 const db=bands.map(b=>sub(newD(b),prime)),dmin=point(min(db.map(x=>x[0]))),dtail=point(min(bands.map((b,i)=>[b,i]).filter(([b])=>b.end===null||b.end>2048).map(([b,i])=>db[i][0])));
 if(dmin[0]<=0n||dtail[0]<dmin[0])throw Error("coercivity");
 const d1=jids.map((j,i)=>{const old=oldBands.filter(b=>b.j<=j).at(-1);if(!overlap(oldD(old),d0[i]))throw Error("old diagonal "+j);return newD(bands.find(b=>j>=b.start&&(b.end===null||j<=b.end)));});
 const poles=jids.map(j=>r["pole"+j]),cp=sum(x0.map((x,i)=>mul(x,poles[i]))),pp=sum(poles.map(sq)),lam=div(cp,pp),x=x0.map((v,i)=>sub(v,mul(lam,poles[i]))),X2=sum(x.map(sq));
 if(!overlap(X2,r.X2)||!overlap(sum(x.map((v,i)=>mul(v,poles[i]))),zero))throw Error("trial reproduction");
 const Delta=x.map((v,i)=>mul(v,sub(d1[i],d0[i]))),cross=sum(Delta.map((v,i)=>mul(v,r["y"+jids[i]]))),Delta2=sum(Delta.map(sq)),xd=sum(x.map((v,i)=>mul(v,Delta[i]))),bd=sum(Delta.map((v,i)=>mul(v,poles[i])));
 const E=sub(r.E,xd),Y2=add(sub(r.Y2,mul(rat("2"),cross)),Delta2),Yb=sub(r.Yb,bd);
 const low=all.filter(j=>j<86),B=sub(Yb,sum(low.map(j=>mul(r["y"+j],r["pole"+j])))),Q=sub(sub(Y2,sum(low.map(j=>sq(r["y"+j])))),div(sq(B),r.pole_st)),projection=div(B,r.pole_st);
 if(Q[1]<0n)throw Error("Q norm");
 const coeff=jids.map((j,i)=>sub(sub(r["y"+j],Delta[i]),mul(projection,poles[i]))),finite=sum(coeff.map(sq)),tail=sub(Q,finite);if(tail[1]<0n)throw Error("tail");
 const A=add(sum(coeff.map((c,i)=>div(up(sq(c)),sub(d1[i],prime)))),div(up(tail),dtail));
 const shift=rat("1/10000000000000"),eps=rat(old.proxy_error_upper),invroot=root(div(one,dmin)),P=add(add(up(E),mul(root(up(X2)),shift)),sq(add(root(up(A)),mul(shift,invroot))));
 if(P[1]<0n)throw Error("negative penalty");
 const T=sq(add(root(up(P)),mul(eps,invroot))),K=rat(old.K[0]),margin=sub(K,T);
 return {status:margin[0]>0n?"PASS_WITNESS":"UNDECIDED",N:2048,trial:"original weighted 067",dmin:pair(dmin),dtail:pair(dtail),X2:pair(X2),delta_norm2:pair(Delta2),energy:pair(E),Q_residual_norm2:pair(Q),finite_norm2:pair(finite),tail_norm2:pair(tail),weighted_residual_upper:str(A[1]),proxy_inverse_upper:str(P[1]),true_inverse_upper:str(T[1]),K_lower:str(K[0]),margin_lower:str(margin[0])};
}
const fs=require('fs');const args=process.argv.slice(2);const result=audit(...args.slice(0,6).map(p=>fs.readFileSync(p,'utf8')));if(args[6])fs.writeFileSync(args[6],JSON.stringify(result,null,2)+'\n');else console.log(JSON.stringify(result,null,2));
