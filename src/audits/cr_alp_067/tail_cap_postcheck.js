// Exact diagnostic: impossibility of passing with fixed response and capped diagonal, not original Weil negativity.
function audit(j,g,l){
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){let c=a%b;a=b;b=c;}return a||1n;}
function q(a,b=1n){if(typeof a==="string"&&a.includes("/")){let [n,d]=a.split("/");return q(BigInt(n),BigInt(d));}a=BigInt(a);b=BigInt(b);if(b<0n){a=-a;b=-b;}const c=gcd(a,b);return[a/c,b/c];}
const a=(x,y)=>q(x[0]*y[1]+y[0]*x[1],x[1]*y[1]),s=(x,y)=>a(x,q(-y[0],y[1])),m=(x,y)=>q(x[0]*y[0],x[1]*y[1]),d=(x,y)=>q(x[0]*y[1],x[1]*y[0]),cmp=(x,y)=>x[0]*y[1]<y[0]*x[1]?-1:x[0]*y[1]>y[0]*x[1]?1:0,str=x=>x[0]+"/"+x[1],floor=x=>x[0]/x[1];
const rows={};for(const line of l.trim().split("\n")){const [key,lo,hi,slash,scale]=line.split(" ");if(slash!=="/"||scale!=="1e15"||key in rows)throw Error("format");rows[key]=[q(lo,10n**15n),q(hi,10n**15n)];}
if(j.theta!=="1/10"||j.N!==2048)throw Error("wrong trial");
const K=q(j.K_lower),E=rows.E[0],Q=q(j.Q_residual_norm2[0]),theta=q(1,10),gamma=q(72,25),sp=q(21,8);
const diag=b=>s(q(floor(m(a(m(s(q(1),theta),gamma),m(theta,q(b.gamma_diagonal_lower))),q(10**9))),10**9),sp);
const ds=g.bands.map(b=>diag(b));if(ds.some((x,i)=>i&&cmp(x,ds[i-1])<0))throw Error("nonmonotone");
const cap=ds.at(-1),eta=q(1,10**13);
// Allow a loose norm bound sqrt(X2)<10 and sqrt(Q)<1 for the rational-shift transport.
// ||e_exact||² >= ||e_rat||² - 2 eta (since ||e_rat||<1).
if(cmp(q(j.Q_residual_norm2[1]),q(1))>=0||cmp(rows.X2[1],q(100))>=0)throw Error("transport");
const Qlo=s(Q,m(q(2),eta)),Elo=s(E,m(q(10),eta)),budget=s(K,Elo);
if(cmp(budget,q(0))<=0)throw Error("budget");
const inf=s(a(Elo,d(Qlo,cap)),K);
if(cmp(inf,q(0))<=0)throw Error("not obstructed");
return {theta:"1/10",cap:str(cap),E_lower_after_shift:str(Elo),Q_norm2_lower_after_shift:str(Qlo),remaining_budget_upper:str(budget),minimum_weighted_certificate_value:str(a(Elo,d(Qlo,cap))),minimum_excess_over_K:str(inf),necessary_effective_floor:str(d(Qlo,budget)),status:"FIXED_TRIAL_CAPPED_DIAGONAL_CERTIFICATE_IMPOSSIBLE",max_theta_strict_coercivity_from_first_band:str(d(s(gamma,sp),s(gamma,q(g.bands[0].gamma_diagonal_lower)))),band_necessary_mass_diagnostics:g.bands.slice(1).map((b,i)=>({through:b.j-1,weight_upper:str(ds[i]),max_mass_below_cutoff_if_pass:str(m(budget,ds[i]))}))};
}
const fs=require('fs');const [j,g,l,out]=process.argv.slice(2);const result=audit(JSON.parse(fs.readFileSync(j,'utf8')),JSON.parse(fs.readFileSync(g,'utf8')),fs.readFileSync(l,'utf8'));if(out)fs.writeFileSync(out,JSON.stringify(result,null,2)+'\n');else console.log(JSON.stringify(result,null,2));
