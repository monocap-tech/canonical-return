// Exact original correction attachment from recovered contracted intervals.
function audit(mid,exact){
const S=10n**30n;function fl(n,d){if(d<=0n)throw Error("denominator");return n>=0n?n/d:-((-n+d-1n)/d);}
function enc(x){let n,d;if(x.includes("/")){[n,d]=x.split("/").map(BigInt);}else{const neg=x.startsWith("-"),[a,b=""]=x.replace("-","").split(".");n=BigInt(a+b)*(neg?-1n:1n);d=10n**BigInt(b.length);}return[fl(n*S,d),-fl(-n*S,d)];}
function row(r){let a=enc(r[0]),b=enc(r[1]);if(a[0]>b[1])throw Error("order");return[a[0],b[1]];}
const abs=x=>x<0n?-x:x,max=(a,b)=>a>b?a:b;
const error=(a,b)=>max(abs(a[0]-b[1]),abs(a[1]-b[0]));
function read(x){const r=x.trim().split("\n").map(y=>y.split(" "));if(r.length!==428||r[0][0]!=="2"||r[2][0]!=="425")throw Error("headers");return {pole:row(r[1]),sine:r.slice(3).map((v,i)=>{if(Number(v[0])!==2*(i+1)||v.length!==3)throw Error("id");return row(v.slice(1));})};}
function sqrt(n){if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));for(;;){const y=(x+n/x)/2n;if(y>=x)return x;x=y;}}
const m=read(mid),e=read(exact),all=m.sine.map((v,i)=>error(v,e.sine[i])),high=all.slice(42),normUpper=sqrt(high.reduce((a,b)=>a+b*b,0n))+1n,allUpper=sqrt(all.reduce((a,b)=>a+b*b,0n))+1n,pole=error(m.pole,e.pole);
if(normUpper>=S/10n**10n)throw Error("attachment tolerance");
return {method:"Q85 projection contractivity on original and archived finite high-sine coefficient differences; pole component killed exactly",all_sine_coordinates:425,high_sine_coordinates:383,scale:String(S),all_sine_l2_error_upper_numerator:String(allUpper),pole_error_upper_numerator:String(pole),projected_correction_l2_error_upper_numerator:String(normUpper),receiver_tolerance:"1/10000000000",status:"PASS_ATTACHMENT"};
}
const fs=require('fs');const [mid,exact,out]=process.argv.slice(2);const result=audit(fs.readFileSync(mid,'utf8'),fs.readFileSync(exact,'utf8'));if(out)fs.writeFileSync(out,JSON.stringify(result,null,2)+'\n');else console.log(JSON.stringify(result,null,2));
