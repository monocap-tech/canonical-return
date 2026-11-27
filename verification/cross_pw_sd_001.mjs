// Exact integer checks for CROSS-PW/SD-001; run with node.
// These are finite model tests, not a prime-discrepancy or RH computation.
const s=100n;
const gap=(s-2n)**2n*(17n*s)-2n*(s+1n)*s**2n;
if(gap<=0n) throw Error("SD scalar normalization obstruction failed");
function k4num(p,q){return 25n*q**6n-147n*p**2n*q**4n+123n*p**4n*q**2n-p**6n;}
const k0=k4num(0n,1n),k9=k4num(9n,10n);
if(k0<=0n||k9>=0n) throw Error("PW model kernel sign check failed");
console.log(JSON.stringify({scalar_gap_numerator:String(gap),k4_at_zero:String(k0),k4_at_9_10_numerator:String(k9),k4_at_9_10_denominator:"1000000",status:"PASS"}));
