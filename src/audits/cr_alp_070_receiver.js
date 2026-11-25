// CR-ALP-070 exact receiver arithmetic; original correction source authentication is external.
function audit(trial,diag,log,newJSON,oldJSON,oldBandsJSON,originalTrial){
 const S=10n**30n,zero=[0n,0n],point=x=>[x,x],one=point(S),min=x=>x.reduce((a,b)=>a<b?a:b),max=x=>x.reduce((a,b)=>a>b?a:b);
 function fl(n,d){if(d<=0n)throw Error("bad denominator");return n>=0n?n/d:-((-n+d-1n)/d);}
 const ce=(n,d)=>-fl(-n,d);
 function rat(s){const t=String(s);if(t.includes("/")){const [n,d]=t.split("/").map(BigInt);return[fl(n*S,d),ce(n*S,d)];}const match=t.match(/^([+-]?)(\d+)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$/);if(!match)throw Error("number "+t);const [,sign,a,b="",e="0"]=match;const exp=Number(e)-b.length;if(!Number.isSafeInteger(exp)||Math.abs(exp)>1000)throw Error("exponent");const n=BigInt(a+b)*(sign==="-"?-1n:1n);const p=exp+30;return p>=0?point(n*10n**BigInt(p)):[fl(n,10n**BigInt(-p)),ce(n,10n**BigInt(-p))];}
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

 // Original correction X: interval coefficients in the ordinary normalized sine basis
 // plus Db times the complete matching pole function. Source linkage is a separate
 // prerequisite; this checker cannot certify that a supplied bank encloses X.
 const words=originalTrial.trim().split(/\s+/);let cursor=0;
 const take=()=>{if(cursor>=words.length)throw Error("short original correction");return words[cursor++];};
 const takeI=()=>{const l=rat(take()),h=rat(take());if(l[0]>h[1])throw Error("reversed input");return[l[0],h[1]];};
 if(take()!=="2")throw Error("wrong sector");const Db=takeI();if(take()!=="425")throw Error("wrong original band");
 const rawX={};for(let j=2;j<=850;j+=2){if(Number(take())!==j)throw Error("original index");rawX[j]=takeI();}
 if(cursor!==words.length)throw Error("trailing original correction data");
 const low=all.filter(j=>j<86),pole=j=>r["pole"+j];
 const X=j=>add(rawX[j]||zero,mul(Db,pole(j)));
 // These are necessary interval consistency checks, not a proof of source identity.
 if(!low.every(j=>overlap(X(j),zero)))throw Error("low sine constraint not enclosed");
 const bnorm=add(r.pole_st,sum(low.map(j=>sq(pole(j)))));
 const poleMoment=add(sum(Object.entries(rawX).map(([j,v])=>mul(v,pole(j)))),mul(Db,bnorm));
 if(!overlap(poleMoment,zero))throw Error("pole constraint not enclosed");
 const gain=d1.map(D=>sub(D,gamma));if(gain.some(a=>a[0]<0n))throw Error("new diagonal does not dominate");
 const tailBands=bands.filter(b=>b.end===null||b.end>2048);
 const tailGain=tailBands.map(b=>sub(newD(b),gamma));
 if(tailGain.some(a=>a[0]<0n))throw Error("negative tail gain");
 const nonneg=a=>{if(a[1]<0n)throw Error("negative squared norm");return[max([0n,a[0]]),a[1]];};
 const B=sub(r.Yb,sum(low.map(j=>mul(r["y"+j],pole(j))))),projection=div(B,r.pole_st);
 const e0=jids.map(j=>sub(r["y"+j],mul(projection,pole(j))));
 const Q0=nonneg(sub(sub(r.Y2,sum(low.map(j=>sq(r["y"+j])))),div(sq(B),r.pole_st)));
 const T0=nonneg(sub(Q0,sum(e0.map(sq))));
 const cTail=nonneg(sub(r.pole_st,sum(poles.map(sq))));
 const finiteGain=sum(x.map((v,i)=>mul(gain[i],sq(add(v,X(jids[i]))))));
 const tailGainLower=mul(mul(point(min(tailGain.map(a=>a[0]))),sq(Db)),cTail);
 const totalGainLower=down(add(finiteGain,tailGainLower));
 const oldAdjustment=sum(x.map((v,i)=>mul(sub(d0[i],gamma),sq(v))));
 const K=rat(old.K[0]),eps=rat(old.proxy_error_upper),eta=rat("1/10000000000000"),normX=root(up(X2));
 const energyError=mul(add(mul(rat("2"),eps),eta),normX);
 const available=sub(add(sub(K,r.E),totalGainLower),add(oldAdjustment,energyError));
 const h=x.map((v,i)=>add(mul(sub(d1[i],d0[i]),v),mul(gain[i],X(jids[i]))));
 const z=e0.map((v,i)=>sub(v,h[i])),df=d1.map(D=>sub(D,prime));
 const mid=a=>point(fl(a[0]+a[1],2n));
 // Any fixed rational t is legal. Use a finite deterministic candidate bank;
 // selecting the least certified upper bound does not rely on optimization.
 const tFinite=mid(div(sum(z.map((v,i)=>div(mul(v,poles[i]),df[i]))),sum(poles.map((v,i)=>div(sq(v),df[i])))));
 const centers=[zero,tFinite,...tailGain.map(a=>mid(sub(zero,mul(a,Db))))];
 const lo=min(centers.map(c=>c[0])),hi=max(centers.map(c=>c[1]));
 for(let i=0;i<=32;i++)centers.push(point(lo+fl((hi-lo)*BigInt(i),32n)));
 const candidates=[...new Set(centers.map(a=>a[0].toString()))].map(s=>point(BigInt(s)));
 const invroot=root(div(one,dmin)),baseTail=root(div(up(T0),dtail));
 let best=null;
 for(const t of candidates){
  const finite=sum(z.map((v,i)=>div(up(sq(sub(v,mul(t,poles[i])))),df[i])));
  const pTail=point(max(tailBands.map((b,i)=>div(sq(add(mul(tailGain[i],Db),t)),sub(newD(b),prime))[1])));
  const tail=sq(add(baseTail,root(mul(up(cTail),pTail))));
  const U=up(add(finite,tail)),truePenalty=sq(add(root(U),mul(add(eps,eta),invroot)));
  const margin=sub(available,truePenalty);
  if(best===null||margin[0]>best.margin[0])best={t,finite,tail,U,truePenalty,margin};
 }
 return {status:best.margin[0]>0n?"CONDITIONAL_RECEIVER_PASS":"UNDECIDED",source_linkage:"REQUIRES_INDEPENDENT_AUTHENTICATION_OF_ORIGINAL_X",N:2048,candidates:candidates.length,dmin:pair(dmin),dtail:pair(dtail),original_pole_coefficient:pair(Db),original_pole_constraint:pair(poleMoment),old_projected_residual_norm2:pair(Q0),old_residual_tail_norm2:pair(T0),pole_tail_norm2:pair(cTail),finite_gain:pair(finiteGain),tail_gain_lower:str(tailGainLower[0]),total_gain_lower:str(totalGainLower[0]),old_diagonal_adjustment:pair(oldAdjustment),energy_transport_upper:str(energyError[1]),available_energy_lower:str(available[0]),pole_subtraction:str(best.t[0]),finite_weighted_upper:str(best.finite[1]),infinite_weighted_upper:str(best.tail[1]),true_inverse_penalty_upper:str(best.truePenalty[1]),margin_lower:str(best.margin[0])};
}

const fs=require('fs');const args=process.argv.slice(2);const result=audit(...args.slice(0,7).map(p=>fs.readFileSync(p,'utf8')));if(args[7])fs.writeFileSync(args[7],JSON.stringify(result,null,2)+'\n');else console.log(JSON.stringify(result,null,2));
