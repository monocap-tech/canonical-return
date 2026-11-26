"""RETURN101-B: independent Arb raw endpoint coefficient regeneration.
Use the physical integration formulas of RETURN098-C, not historical
mpmath Gamma or raw coefficient evaluators. Strict containment required.
"""
import sys,json,hashlib
from fractions import Fraction as F
from pathlib import Path
from flint import arb,acb,fmpq,ctx
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
ctx.prec=384;Q=10**55
def A(v):
 v=F(v);return arb(fmpq(v.numerator,v.denominator))
def load(p,h):
 b=(root/p).read_bytes();assert hashlib.sha256(b).hexdigest()==h,p
 return json.loads(b)
def enc(v):
 l=(v.lower()*Q).floor().unique_fmpz();u=(v.upper()*Q).ceil().unique_fmpz()
 assert l is not None and u is not None
 return [int(l),int(u)]
w=load('return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
raw=load('return059/RETURN059_raw_receipt.json','2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8')
mom=load('return066/RETURN066_moments.json','b9aa229e7147b8d47d56b0cbfec290255af557cfb99260624d1c4b576369d184')
pi=arb.pi();a=A(F(5,4));D=2*a;half=A(F(1,2));one=arb(1)
psi0=acb(a).digamma().real
shifts=[(arb(n).log(),arb(p).log()/arb(n).sqrt()) for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]]
bs=[A(F(4*m+1,2)) for m in range(1,49)]
es=[(-b*D).exp() for b in bs]
BB=A(F(197,2));decay=(-BB*D).exp()/(1-(-2*D).exp())
assert decay/(BB*BB-half*half)<arb('1e-101')
assert decay/(BB-half)<arb('1e-101')
pad=arb(0,arb('1e-100'))
ks0=sorted(j for j in set(w['even_exceptional']+w['unpaid_odd']) if j>0)
assert len(ks0)==42 and len(raw['ks'])==44
def coefficient(par,typ,k):
 if typ:
  assert par==1 and k==0
  ell=half*half;b=(a*half).sinh();v=half*(a*half).cosh()
  p0=acb(a+half/2).digamma().real;p1=acb(a-half/2).digamma().real
  g=(p0+p1)/2-psi0;h=(p0-p1)/(4*half)
  M=a.sinh()-a
  def val(t):return (half*t).sinh(),half*(half*t).cosh()
 elif k==0:
  assert par==0
  return [arb(0),one,arb(0),arb(0),4*(a*half).sinh()*(a*half).exp()]
 else:
  r=2*pi*k/5;ell=-r*r
  co=arb([1,0,-1,0][k%4]);si=arb([0,1,0,-1][k%4])
  b=si if par else co;v=r*co if par else -r*si
  z=acb(a,r/2).digamma();g=z.real-psi0;h=z.imag/(2*r)
  M=((half*b-v)*(half*a).exp()-(-1 if par else 1)*(half*b+v)*(-half*a).exp())/(half*half-ell)
  def val(t):
   s,c=(r*t).sin_cos()
   return (s,r*c) if par else (c,-r*s)
 e=sum(t/(beta*beta-ell) for t,beta in zip(es,bs))+pad
 f=sum(t*beta/(beta*beta-ell) for t,beta in zip(es,bs))+pad
 if par:
  c=b*(-g/2+f)+v*(-h+e)+sum(weight*val(a-s)[0] for s,weight in shifts)
  d=2*M*(a*half).exp()
 else:
  c=ell*b*(h-e)+v*(g/2-f)-sum(weight*val(a-s)[1] for s,weight in shifts)
  d=M*(a*half).exp()
 return [ell,b,v,c,d]
records=[];allrows=[];fail=[]
for par in (0,1):
 ids=([(0,0)]if par==0 else [])+[(0,j)for j in ks0]+([(1,0)]if par else [])+[(0,j)for j in raw['ks']]
 assert len(ids)==len(mom['blocks'][par]['raw_coeff'])==87
 rows=[];contained=0
 for i,(typ,k) in enumerate(ids):
  z=coefficient(par,typ,k);v=[enc(x) for x in z];ref=mom['blocks'][par]['raw_coeff'][i]
  assert len(v)==len(ref)==5
  for j,((l,u),(a0,b0)) in enumerate(zip(v,ref)):
   ok=a0<=l<=u<=b0;contained+=int(ok)
   if not ok:fail.append([par,typ,k,j,[l,u],[a0,b0]])
  rows.append(v)
 records.append({'parity':par,'raw_generators':87,'entries':435,'contained':contained})
 allrows.append(rows)
 print('Arb coefficient containment',records[-1],flush=True)
ok=all(r['contained']==r['entries'] for r in records)
(out/'receipt.json').write_text(json.dumps({'status':'PASS'if ok else 'NOT_CERTIFIED','scope':'independent Arb raw ell,b,v,c,d evaluation; physical endpoint formulas, 48-term directed decay bound, no historical evaluator','scale':Q,'records':records,'raw_coeff':allrows,'mismatches':fail[:20]},indent=2)+'\n')
assert ok,fail[:3]
print('RETURN101-B ALL 870 RAW COEFFICIENTS PASS',flush=True)
