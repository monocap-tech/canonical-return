"""RETURN101-C: Arb retained Gram attachment, independent numerical primitives.
Use RETURN098-A's physical exponential formulas, not historical evaluators.
Compare strict containment in the archived RETURN086 Hraw and RETURN093-B
alternate Hraw; no new terminal solve or proof promotion.
"""
import sys,json,hashlib,math
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache as cache
from flint import arb,acb,fmpq,ctx
sys.set_int_max_str_digits(0)
root=Path(sys.argv[1]);bank=Path(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=True)
ctx.prec=448;S=10**85
def A(v):
 v=F(v);return arb(fmpq(v.numerator,v.denominator))
def load(p,h=None):
 data=p.read_bytes()
 if h:assert hashlib.sha256(data).hexdigest()==h,p
 return json.loads(data)
raw=load(root/'return059/RETURN059_raw_receipt.json','2433554360bde233df3718c68323c5e0d007e4976eb1c5e4a8cf7b842648ded8')
ledger=load(root/'return047/RETURN047_payment_receipt.json','def4c4d30dbd4467d9eee44e18edfefa407c349928c9fd21de4851e87650cc4f')
a=A(F(5,4));D=2*a;half=A(F(1,2));pi=arb.pi();zero=(0,0);one=(1,0)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
@cache(None)
def z(x):return acb(A(F(x[0],2)),2*pi*x[1]/5)
@cache(None)
def exa(x):return (A(F(x[0],2))*a).exp()*[acb(1),acb(0,1),acb(-1),acb(0,-1)][x[1]%4]
@cache(None)
def I(x):return acb(D) if x==zero else (exa(x)-exa(neg(x)))/z(x)
@cache(None)
def Ip(x):return acb(0) if x==zero else (a*(exa(x)+exa(neg(x)))-I(x))/z(x)
@cache(None)
def psi(x):return (acb(a)+z(x)/2).digamma()
# Independent exact Bernoulli recurrence for the rare coincident-node
# trigamma; periodic-Bernoulli remainder bound from RETURN093-B.
B=[F(1)]
for n in range(1,81):B.append(-sum(F(math.comb(n+1,k))*B[k] for k in range(n))/F(n+1))
N=256;P=40
trerr=A(abs(B[2*P])/F(N**(2*P+1)))
def pad(v,e):return acb(v.real+arb(0,e.upper()),v.imag+arb(0,e.upper()))
@cache(None)
def tri(x):
 v=acb(a)+z(x)/2;w=v+N
 assert w.real>N
 return pad(1/w+1/(2*w*w)+sum(A(B[2*j])/w**(2*j+1) for j in range(1,P+1))+sum(1/(v+j)**2 for j in range(N)),trerr)
@cache(None)
def ss(x,y):return tri(x)/4 if x==y else (psi(y)-psi(x))/(2*z(sub(y,x)))
betas=[A(F(4*m+1,2)) for m in range(1,65)]
weights=[(-b*D).exp() for b in betas]
BB=A(F(261,2));decerr=(-BB*D).exp()/((BB-half)**2*(1-(-2*D).exp()))
@cache(None)
def den(x):return [1/(b+z(x)) for b in betas]
@cache(None)
def dec(x,y):return pad(sum(w*u*v for w,u,v in zip(weights,den(x),den(y))),decerr)
@cache(None)
def triangle(x,y):
 t=add(y,one)
 if t==zero:return Ip(sub(x,one))+a*I(sub(x,one))
 return (I(add(x,y))-exa(neg(t))*I(sub(x,one)))/z(t)
G0=acb(A(F(1,4))).digamma().real-pi.log()
@cache(None)
def gpart(x):return (psi(x)+psi(neg(x)))/2-psi(zero)
shifts=[(arb(n).log(),arb(p).log()/arb(n).sqrt()) for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]]
assert arb(11).log()<D<arb(13).log()
@cache(None)
def shiftI(x,j):
 s,w=shifts[j];v=z(x)*(a-s/2)
 return acb(D-s) if x==zero else (v.exp()-(-v).exp())/z(x)
@cache(None)
def pair(x,j):
 s,w=shifts[j];v=z(x)*s/2
 return w*(v.exp()+(-v).exp())
@cache(None)
def entry(x,y):
 u=add(x,y);v=sub(y,x)
 gamma=(G0+gpart(x)+4)*I(u)-triangle(x,y)-triangle(y,x)
 gamma+=exa(neg(u))*ss(x,neg(y))+exa(u)*ss(neg(x),y)
 gamma-=exa(v)*dec(x,neg(y))+exa(neg(v))*dec(neg(x),y)
 return gamma-sum(pair(v,j)*shiftI(u,j) for j in range(8))+I(add(x,one))*I(sub(y,one))+I(sub(x,one))*I(add(y,one))
def enc(v):
 l=(v.lower()*S).floor().unique_fmpz();u=(v.upper()*S).ceil().unique_fmpz()
 assert l is not None and u is not None
 return [int(l),int(u)]
def modes(k,p):
 return [((0,k),acb(half) if not p else acb(0,-half)),((0,-k),acb(half) if not p else acb(0,half))]
ks=sorted(j for j in set(ledger['even_exceptional']+ledger['unpaid_odd']) if j>0)
assert len(ks)==42 and len(raw['ks'])==44
rec=[];fail=[]
for p in (0,1):
 cols=([[(zero,acb(1))]] if not p else [])+[modes(k,p) for k in ks]
 if p:cols.append([((1,0),acb(half)),((-1,0),acb(-half))])
 cols += [modes(k,p) for k in raw['ks']];assert len(cols)==87
 H=[[None]*87 for i in range(87)]
 for i,c in enumerate(cols):
  for j,d in enumerate(cols[:i+1]):
   v=sum(w*t*entry(x,y) for x,w in c for y,t in d)
   assert v.imag.contains(0),(p,i,j)
   H[i][j]=H[j][i]=enc(v.real)
  if (i+1)%16==0:print('Arb retained rows',p,i+1,flush=True)
 for suffix in ('','_alt'):
  name=f'matrices_{p}{suffix}.json'
  for label,path in [('original',root/'return086'/name),('alternate',bank/('reconstructed_'+name))]:
   d=load(path);Q=d['HS'];assert Q==10**75 and d['C']==raw['blocks'][p]['C'] and d['CS']==raw['T']
   contained=0
   for i in range(87):
    for j in range(87):
     l,u=H[i][j];lo,hi=d['Hraw'][i][j]
     ok=lo*S<=l*Q<=u*Q<=hi*S
     contained+=int(ok)
     if not ok and len(fail)<20:fail.append([p,suffix,label,i,j,[l,u],[lo,hi]])
   rec.append({'parity':p,'variant':suffix or 'generator','reference':label,'entries':7569,'contained':contained})
   print('Arb retained comparison',rec[-1],flush=True)
 (out/f'H_{p}.json').write_text(json.dumps({'parity':p,'scale':S,'Hraw':H},separators=(',',':'))+'\n')
ok=all(r['contained']==r['entries'] for r in rec)
(out/'receipt.json').write_text(json.dumps({'status':'PASS'if ok else 'NOT_CERTIFIED','scope':'independent Arb physical retained Gram, strict containment in both original and alternate terminal inputs; no terminal rerun','records':rec,'mismatch_samples':fail},indent=2)+'\n')
assert ok,fail[:3]
print('RETURN101-C ALL RETAINED GRAM CONTAINMENT PASS',sum(r['entries'] for r in rec),flush=True)
