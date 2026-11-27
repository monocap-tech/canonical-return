// CROSS-PW/SD-003 finite symbolic identities with exact integer arithmetic.
// Not evaluation of the actual zero spectrum or a positivity certificate.
for(let p=-20n;p<=20n;p++){const q=13n;if(q*q-(2n*p-q)**2n!==4n*p*(q-p))throw Error("spectral denominator");}
for(let n=1n;n<=128n;n++){if(1n-(-4n*n-1n)**2n!==-8n*n*(2n*n+1n))throw Error("trivial denominator");}
// A real dilation of 2 on one channel cannot preserve a positive metric on that eigenline.
if(2n**2n-1n!==3n)throw Error("metric obstruction");
console.log(JSON.stringify({spectral_denominator_checks:41,trivial_denominator_checks:128,metric_defect:3,status:"PASS"}));
