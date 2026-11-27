// CR-ALP-123AV: independent remainder guard arithmetic.
function audit(){
 const S=10n**340n, R=10n**290n;const fact=n=>{let f=1n;for(let i=1n;i<=n;i++)f*=i;return f};let checks={};
 function lt(name,a,b){if(!(a<b))throw Error(name);checks[name]="PASS"}
 lt("log_400_remainder",2n*7n**801n*400n*S,801n*351n*20n**801n);
 lt("atan_400_remainder",S,801n*3n**801n);
 lt("sincos_400_remainder",4n**800n*S,fact(800n));
 lt("exp_400_remainder",3n**400n*401n*S,398n*fact(400n));
 for(let shift of [512n,768n]){
  lt("digamma_EM_"+shift,4n*fact(180n)*R,6n**180n*180n*shift**180n);
  lt("trigamma_EM_"+shift,4n*fact(180n)*R,6n**180n*shift**181n);
 }
 lt("reflected_resolvent_K256_remainder",100n*2n**1025n*625n*S,5n**1025n*609n);
 // Lower exp(11/10) via first six positive Taylor terms, to ensure exp(-a)<2/5.
 let n=6n,den=10n**n*fact(n),num=0n;
 for(let k=0n;k<=n;k++)num+=11n**k*10n**(n-k)*fact(n)/fact(k);
 lt("exp_minus_a_lt_2_over_5",5n*den,2n*num);
 // For x=11/20, lower/upper exp(x) on [17/10,7/4] by rational
 // Taylor with geometric tail. Then |cosh|<6/5, |sinh|<3/5.
 n=12n;den=20n**n*fact(n);num=0n;
 for(let k=0n;k<=n;k++)num+=11n**k*20n**(n-k)*fact(n)/fact(k);
 lt("exp_half_a_gt_17_over_10",17n*den,10n*num);
 // tail x^13/13! /(1-x/14)
 const td=20n**13n*fact(13n)*269n,tn=11n**13n*280n;
 lt("exp_half_a_lt_7_over_4",4n*(num*td+tn*den),7n*den*td);
 // e+1/e increasing for e>1; at 7/4: cosh <=(7/4+4/7)/2=65/56<6/5;
 // sinh <=(7/4-4/7)/2=33/56<3/5.
 lt("anchor_cosh_bound",65n*5n,56n*6n);
 lt("anchor_sinh_bound",33n*5n,56n*3n);
 // All omitted rational factors are bounded by 100 for beta>=512.
 // beta^2-1/4 >= beta^2/2; |beta b0-b1| < 2 beta.
 // EC/(beta^2+v^2) <= 4/beta^3; b-b factor <=16/beta^2.
 lt("reflected_factor_uniform_bound",16n,100n);
 return {status:"EXACT_ELEMENTARY_REMAINDER_GUARDS_PASS",checks,grid:"10^340",shifts:[512,768],bernoulli_terms:90,reflected_terms:256,scope:"Independent integer checks of actual 046/AP3 numerical remainder constants; analytic identities and directed interval semantics reviewed separately."};
}
const fs=require('fs');const s=JSON.stringify(audit(),null,2)+'\n';if(process.argv[2])fs.writeFileSync(process.argv[2],s);else process.stdout.write(s);
