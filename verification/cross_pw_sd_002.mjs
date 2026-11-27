// CROSS-PW/SD-002: exact model diagnostics, no actual zeta or prime data.
const q=100n,p=1n;
// rank lower fraction ((1-e)^2/(1+e)) with e=p/q.
const num=(q-p)**2n, den=q*(q+p);
// For a rank <= 10 fixed subspace, target dimension 11 is impossible.
if(num*11n<=den*10n)throw Error("rank obstruction failed");
// Exact finite exponential-mode analogue: roots r=2,1/2 of J_{n+1}-(5/2)J_n+J_{n-1}=0.
const J=[1n,2n,4n,8n,16n,32n];
for(let i=1;i+1<J.length;i++)if(2n*J[i+1]-5n*J[i]+2n*J[i-1]!==0n)throw Error("mode residual");
console.log(JSON.stringify({rank_lower_numerator:String(num),rank_lower_denominator:String(den),target_dimension:11,rank_cap:10,mode_residuals:4,status:"PASS"}));
